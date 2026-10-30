#!/usr/bin/env python3
"""准备有来源证据的生产图谱；源FSL只读，患者数据不参与。"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import inspect
import json
import os
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

import nibabel as nib
import nilearn
import numpy as np
import requests
import yaml
from nilearn.datasets import fetch_atlas_harvard_oxford
from nilearn.image import resample_to_img

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wmh_gmv.config import load_config


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else root / path).resolve()


def read_sources(path: Path) -> dict:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"version", "fsl_root", "fsl_version_file", "expected_fsl_version", "working_cache", "output_manifest",
                "allow_lateralized_download", "download_timeout_seconds", "allow_library_split_fallback",
                "overwrite_existing", "expected_spatial_unit", "expected_voxel_sizes", "expected_template",
                "maximum_opposite_hemisphere_fraction", "atlases"}
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError(f"Atlas-source config keys must exactly match: {sorted(required)}")
    if value["version"] != 1 or value["expected_spatial_unit"] != "mm":
        raise ValueError("Only source schema version 1 and explicit mm units are supported")
    if float(value["download_timeout_seconds"]) <= 0:
        raise ValueError("download_timeout_seconds must be positive")
    if not 0 <= float(value["maximum_opposite_hemisphere_fraction"]) < .5:
        raise ValueError("maximum_opposite_hemisphere_fraction must be in [0, 0.5)")
    base = {"target_name", "method", "source_image", "source_labels", "labels_index_offset", "expected_shape"}
    for atlas in value["atlases"]:
        extra = {"public_atlas_name", "fallback_atlas_name", "public_source_url"} if atlas["method"] == "nilearn_lateralized" else set()
        if set(atlas) != base | extra or atlas["method"] not in ("copy", "nilearn_lateralized"):
            raise ValueError(f"Invalid atlas source definition: {atlas.get('target_name')}")
    return value


def xml_labels(path: Path, offset: int, background: bool) -> dict[int, str]:
    output = {0: "Background"} if background else {}
    for element in ET.parse(path).findall(".//label"):
        index = int(element.attrib["index"]) + offset
        if index in output:
            raise ValueError(f"Repeated label index in {path}: {index}")
        output[index] = str(element.text).strip()
    return output


def write_labels(path: Path, labels: dict[int, str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["index", "name"])
        writer.writerows(sorted(labels.items()))


@contextlib.contextmanager
def network_policy(allowed: bool, timeout: float):
    original = requests.Session.send

    def send(session, request, **kwargs):
        if not allowed:
            raise RuntimeError("Network downloads disabled for local atlas preparation")
        kwargs["timeout"] = timeout
        return original(session, request, **kwargs)

    with patch.object(requests.Session, "send", send):
        yield


def validate_image(image, source: dict, options: dict, config: dict, target: dict) -> dict:
    if tuple(image.shape) != tuple(source["expected_shape"]):
        raise ValueError(f"Unexpected shape for {target['name']}: {image.shape}")
    if image.header.get_xyzt_units()[0] != options["expected_spatial_unit"]:
        raise ValueError(f"Atlas unit is not explicitly mm: {target['name']}; source headers will not be rewritten")
    affine = np.asarray(image.affine)
    if not np.isfinite(affine).all() or abs(np.linalg.det(affine[:3, :3])) <= np.finfo(float).eps:
        raise ValueError(f"Invalid affine: {target['name']}")
    if not np.allclose(nib.affines.voxel_sizes(affine), options["expected_voxel_sizes"], rtol=0,
                       atol=config["analysis"]["acute"]["grid_affine_tolerance"]):
        raise ValueError(f"Unexpected voxel sizes: {target['name']}")
    if target["template"] != options["expected_template"] or target["template"] != config["analysis"]["acute"]["template"]:
        raise ValueError(f"Template declarations disagree: {target['name']}")
    values = image.get_fdata(dtype=np.float32)
    if not np.isfinite(values).all():
        raise ValueError(f"Atlas contains nonfinite values: {target['name']}")
    if target["kind"] == "labels":
        if not np.allclose(values, np.rint(values), rtol=0, atol=config["analysis"]["acute"]["binary_tolerance"]):
            raise ValueError(f"Noninteger deterministic atlas: {target['name']}")
        present = [int(x) for x in np.unique(values)]
    else:
        scale = config["analysis"]["acute"]["probability_scale"]
        if values.min() < 0 or values.max() > scale:
            raise ValueError(f"Probability atlas exceeds configured scale: {target['name']}")
        present = list(range(values.shape[3]))
    return {"shape": list(image.shape), "voxel_sizes_mm": nib.affines.voxel_sizes(affine).tolist(),
            "spatial_unit": image.header.get_xyzt_units()[0], "affine": affine.tolist(), "present_indices": present,
            "coordinate_template": target["template"]}


def lateralized(source: dict, source_image: Path, source_xml: Path, options: dict, cache: Path):
    public_error = None
    if options["allow_lateralized_download"]:
        public_cache = cache / "public"
        try:
            with network_policy(True, float(options["download_timeout_seconds"])):
                atlas = fetch_atlas_harvard_oxford(source["public_atlas_name"], data_dir=str(public_cache),
                                                   symmetric_split=False, resume=True, verbose=0)
            source_path = Path(atlas.filename)
            image = nib.load(source_path)
            label_path = public_cache / "fsl/data/atlases/HarvardOxford-Cortical-Lateralized.xml"
            metadata = {"derivation": "official_pre_lateralized_nilearn_distribution", "source_type": "downloaded_FSL5_distribution",
                        "source_url": source["public_source_url"], "source_image_sha256": sha256(source_path),
                        "source_labels_sha256": sha256(label_path), "nilearn_version": nilearn.__version__,
                        "source_image": source["public_source_url"] + "#" + str(source_path.relative_to(public_cache / "fsl")),
                        "source_labels": source["public_source_url"] + "#data/atlases/HarvardOxford-Cortical-Lateralized.xml",
                        "source_fsl_version": "FSL5_distribution_provided_by_Nilearn", "_local_source_image": str(source_path)}
            return image, dict(enumerate(atlas.labels)), metadata
        except Exception as error:
            public_error = f"{type(error).__name__}: {error}"
            if not options["allow_library_split_fallback"]:
                raise
    if not options["allow_library_split_fallback"]:
        raise RuntimeError("No locally installed lateralized atlas; library fallback is disabled")
    # 拷贝两项源文件到隔离缓存，禁止Nilearn改写原FSL目录或联网。
    local_cache = cache / "local"
    local_image = local_cache / "fsl" / source["source_image"]
    local_labels = local_cache / "fsl" / source["source_labels"]
    local_image.parent.mkdir(parents=True, exist_ok=True)
    local_labels.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_image, local_image)
    shutil.copyfile(source_xml, local_labels)
    with network_policy(False, float(options["download_timeout_seconds"])):
        atlas = fetch_atlas_harvard_oxford(source["fallback_atlas_name"], data_dir=str(local_cache),
                                           symmetric_split=False, resume=False, verbose=0)
    image = atlas.maps if hasattr(atlas.maps, "affine") else nib.load(atlas.maps)
    # 此非公开辅助函数来自已锁定Nilearn版本，记录源代码哈希；不实现另一套左右划分算法。
    from nilearn.datasets.atlas import _compute_symmetric_split
    if sorted(np.unique(image.get_fdata()).astype(int).tolist()) != list(range(len(atlas.labels))):
        raise ValueError("Library split requires contiguous present labels; refusing ambiguous label-name pairing")
    values, names = _compute_symmetric_split("HarvardOxford", image, list(atlas.labels))
    source_header = nib.load(source_image).header.copy()
    output = nib.Nifti1Image(values.astype(np.int16), image.affine, header=source_header)
    metadata = {"derivation": "nilearn_library_symmetric_split", "source_type": "installed_FSL",
                "source_image_sha256": sha256(source_image), "source_labels_sha256": sha256(source_xml),
                "nilearn_version": nilearn.__version__, "library_function": "nilearn.datasets.atlas._compute_symmetric_split",
                "library_function_sha256": hashlib.sha256(inspect.getsource(_compute_symmetric_split).encode()).hexdigest(),
                "public_download_error": public_error, "header_policy": "copied_known_mm_source_header; no source mutation"}
    return output, dict(enumerate(names)), metadata


def validate_left_right(image, labels: dict[int, str], tolerance: float, maximum_opposite_fraction: float) -> dict:
    values = image.get_fdata()
    checked = 0
    details = []
    for index, name in labels.items():
        if not (name.startswith("Left ") or name.startswith("Right ")):
            continue
        voxels = np.argwhere(values == index)
        if not len(voxels):
            continue
        x = nib.affines.apply_affine(image.affine, voxels)[:, 0]
        opposite = x > tolerance if name.startswith("Left ") else x < -tolerance
        wrong_centroid = x.mean() >= 0 if name.startswith("Left ") else x.mean() <= 0
        if wrong_centroid or opposite.mean() > maximum_opposite_fraction:
            raise ValueError(f"Lateralized label violates world-coordinate hemisphere: {index} {name}; x_min={x.min()}, x_max={x.max()}, x_mean={x.mean()}, opposite_fraction={opposite.mean()}")
        details.append({"index": index, "name": name, "voxel_count": len(x), "x_min_mm": float(x.min()),
                        "x_max_mm": float(x.max()), "x_mean_mm": float(x.mean()),
                        "opposite_hemisphere_voxels": int(opposite.sum()), "opposite_hemisphere_fraction": float(opposite.mean())})
        checked += 1
    if checked == 0:
        raise ValueError("No left/right labels could be verified")
    return {"left_right_world_x_check": "passed_centroid_and_fraction", "lateralized_nonempty_labels_checked": checked,
            "maximum_allowed_opposite_fraction": maximum_opposite_fraction, "lateralization_qc": details}


def publish(source: Path, destination: Path, overwrite: bool) -> str:
    if destination.exists():
        if sha256(source) == sha256(destination):
            return "unchanged"
        if not overwrite:
            raise FileExistsError(f"Refusing to replace different existing atlas resource: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".partial")
    shutil.copyfile(source, temporary)
    os.replace(temporary, destination)
    return "created"


def prepare(config_path: Path, sources_path: Path) -> dict:
    config = load_config(config_path)
    options = read_sources(sources_path)
    root = Path(config["project"]["project_root"])
    fsl_root = resolve(root, options["fsl_root"])
    version_path = resolve(fsl_root, options["fsl_version_file"])
    version = version_path.read_text().strip()
    if version != options["expected_fsl_version"]:
        raise ValueError(f"FSL version mismatch: expected {options['expected_fsl_version']}, observed {version}")
    targets = {atlas["name"]: atlas for atlas in config["paths"]["atlases"]}
    if set(targets) != {atlas["target_name"] for atlas in options["atlases"]}:
        raise ValueError("Atlas source list must match all configured production atlas names")
    cache_root = resolve(root, options["working_cache"])
    cache_root.mkdir(parents=True, exist_ok=True)
    records = []
    with tempfile.TemporaryDirectory(prefix="prepare-", dir=cache_root) as temporary:
        work = Path(temporary)
        staged = []
        for source in options["atlases"]:
            target = targets[source["target_name"]]
            input_image = resolve(fsl_root, source["source_image"])
            input_xml = resolve(fsl_root, source["source_labels"])
            before = {"image": sha256(input_image), "labels": sha256(input_xml)}
            original = nib.load(input_image)
            validate_image(original, source, options, config, target)
            output_image = work / (target["name"] + ".nii.gz")
            output_labels = work / (target["name"] + ".tsv")
            if source["method"] == "copy":
                shutil.copyfile(input_image, output_image)
                labels = xml_labels(input_xml, int(source["labels_index_offset"]), target["kind"] == "labels")
                provenance = {"derivation": "byte_identical_copy", "source_type": "installed_FSL",
                              "source_image_sha256": before["image"], "source_labels_sha256": before["labels"]}
            else:
                image, labels, provenance = lateralized(source, input_image, input_xml, options, cache_root / "nilearn-cache")
                # 显式还原到配置源图谱的FSL采样网格；不更改模板，不执行配准。
                resampled = image.shape != original.shape or not np.allclose(image.affine, original.affine, rtol=0,
                                                atol=config["analysis"]["acute"]["grid_affine_tolerance"])
                if resampled:
                    image = resample_to_img(image, original, interpolation="nearest", force_resample=True, copy_header=True)
                local_source = provenance.pop("_local_source_image", None)
                if local_source and not resampled:
                    shutil.copyfile(local_source, output_image)
                else:
                    nib.save(image, output_image)
                provenance["grid_resampled"] = bool(resampled)
                provenance.update(validate_left_right(image, labels, config["analysis"]["acute"]["grid_affine_tolerance"],
                                                        float(options["maximum_opposite_hemisphere_fraction"])))
            write_labels(output_labels, labels)
            validation = validate_image(nib.load(output_image), source, options, config, target)
            if not set(validation["present_indices"]).issubset(labels):
                raise ValueError(f"Output dictionary is missing actual indices: {target['name']}")
            selected = set(target["include_labels"] or labels) - set(target["exclude_labels"])
            if not selected.issubset(labels):
                raise ValueError(f"Unknown selected label: {target['name']}")
            if before != {"image": sha256(input_image), "labels": sha256(input_xml)}:
                raise RuntimeError("Installed source files changed during preparation")
            if target.get("sha256") and target["sha256"].lower() != sha256(output_image):
                raise ValueError(f"Prepared image does not match configured hash: {target['name']}")
            record = {"name": target["name"], "path": target["path"], "labels_path": target["labels"],
                      "project_relative_path": str(Path(target["path"]).relative_to(root)),
                      "project_relative_labels_path": str(Path(target["labels"]).relative_to(root)),
                      "license": target["license"], "license_source": target["source"],
                      "source_image": str(input_image), "source_labels": str(input_xml), "source_fsl_version": version,
                      "source_xml_index_offset": source["labels_index_offset"], **provenance, **validation,
                      "installed_reference_image": str(input_image), "installed_reference_image_sha256": before["image"],
                      "dictionary_label_count": len(labels), "selected_label_count": len(selected),
                      "output_sha256": sha256(output_image), "labels_sha256": sha256(output_labels)}
            staged.append((output_image, Path(target["path"]), output_labels, Path(target["labels"]), record))
        # 所有图谱先通过检查，再发布；已存在不同内容默认拒绝覆盖。
        for image, destination, labels, label_destination, record in staged:
            for src, dst in ((image, destination), (labels, label_destination)):
                if dst.exists() and sha256(src) != sha256(dst) and not options["overwrite_existing"]:
                    raise FileExistsError(f"Refusing replacement: {dst}")
        for image, destination, labels, label_destination, record in staged:
            record["publication"] = publish(image, destination, options["overwrite_existing"])
            publish(labels, label_destination, options["overwrite_existing"])
            records.append(record)
    manifest = {"created_at_utc": datetime.now(timezone.utc).isoformat(), "sources_config_sha256": sha256(sources_path),
                "main_config": str(config_path.resolve()), "fsl_root": str(fsl_root), "fsl_version": version,
                "nilearn_version": nilearn.__version__, "nibabel_version": nib.__version__,
                "source_files_modified": False, "atlases": records}
    manifest_path = resolve(root, options["output_manifest"])
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_manifest = manifest_path.with_suffix(".json.partial")
    tmp_manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    os.replace(tmp_manifest, manifest_path)
    return {"manifest": str(manifest_path), "fsl_version": version,
            "atlases": [{"name": record["name"], "labels": record["selected_label_count"],
                         "derivation": record["derivation"], "sha256": record["output_sha256"]} for record in records]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--sources", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.config, args.sources), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
