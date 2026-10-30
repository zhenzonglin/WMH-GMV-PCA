"""既有急性病灶的图谱负荷提取；不修改WMH或GMV，也不执行配准。"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np
from nilearn.image import resample_to_img


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _labels(path: Path) -> dict[int, str]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        index_key = next((key for key in ("index", "label", "value") if key in fields), None)
        name_key = next((key for key in ("name", "label_name", "region") if key in fields), None)
        if index_key is None or name_key is None:
            raise ValueError(f"Atlas label table requires index/name columns: {path}")
        result: dict[int, str] = {}
        for row in reader:
            index = int(row[index_key])
            if index in result:
                raise ValueError(f"Repeated atlas label index {index}: {path}")
            result[index] = str(row[name_key])
    return result


def _geometry(image: nib.spatialimages.SpatialImage, name: str, expected_unit: str) -> None:
    if len(image.shape) not in (3, 4):
        raise ValueError(f"{name} must be a 3D label or 4D probability image")
    affine = np.asarray(image.affine)
    if not np.isfinite(affine).all() or abs(np.linalg.det(affine[:3, :3])) <= np.finfo(float).eps:
        raise ValueError(f"{name} has invalid/singular spatial affine")
    units = image.header.get_xyzt_units()[0]
    if units != expected_unit:
        raise ValueError(f"{name} spatial units must explicitly be {expected_unit}: {units}; no unit reinterpretation is performed")


def _on_mask_grid(image: nib.spatialimages.SpatialImage, target: nib.spatialimages.SpatialImage,
                  interpolation: str, allow: bool, affine_tolerance: float) -> np.ndarray:
    same_grid = image.shape[:3] == target.shape[:3] and np.allclose(image.affine, target.affine, rtol=0, atol=affine_tolerance)
    if same_grid:
        return image.get_fdata(dtype=np.float32)
    if not allow:
        raise ValueError("Atlas and lesion grids differ; grid resampling is disabled")
    # 调整的是同一模板下的采样网格，绝不能将此解释为配准。
    result = resample_to_img(image, target, interpolation=interpolation,
                             force_resample=True, copy_header=True)
    return result.get_fdata(dtype=np.float32)


def _write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if value is None else value for key, value in row.items()})


def extract_acute(subject: dict[str, Any], config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """返回任务回执；缺输入等待，坏几何/非二值数据抛错，交由worker记录失败。"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(subject.get(key) in (None, "") for key in ("subject_key", "participant_id", "cohort")):
        raise ValueError("Acute feature extraction requires nonempty subject identifiers and cohort")
    cfg = config["analysis"]["acute"]
    expected_unit = str(cfg["expected_spatial_unit"])
    if expected_unit != "mm":
        raise ValueError("Acute volume extraction supports only explicit millimetre (mm) spatial units")
    atlases = config["paths"]["atlases"]
    mask_value = subject.get("lesion_mask")
    if mask_value is None or str(mask_value).strip() in ("", "NA", "nan"):
        return {"status": "waiting_input", "reason": "Missing acute lesion mask", "artifacts": [], "n_subjects": 1}
    mask_path = Path(str(mask_value))
    missing = [str(mask_path)] if not mask_path.is_file() else []
    for atlas in atlases:
        for key in ("path", "labels"):
            path = Path(atlas[key])
            if not path.is_file():
                missing.append(str(path))
    if missing or not atlases:
        return {"status": "waiting_input", "reason": "Missing lesion/atlas inputs: " + ", ".join(missing),
                "artifacts": [], "n_subjects": 1}
    declared = subject.get("lesion_template")
    if not declared or str(declared) in ("NA", "nan"):
        return {"status": "waiting_input", "reason": "Missing lesion_template declaration; grid shape is not registration evidence",
                "artifacts": [], "n_subjects": 1}
    expected_template = str(cfg["template"])
    if str(declared) != expected_template:
        raise ValueError(f"Lesion template {declared!r} differs from configured template {expected_template!r}")
    if expected_template.lower().startswith("synthetic") and config["project"]["mode"] != "synthetic":
        raise ValueError("Synthetic template is only allowed in project.mode=synthetic")
    mask_image = nib.load(str(mask_path))
    _geometry(mask_image, "Lesion", expected_unit)
    if len(mask_image.shape) != 3:
        raise ValueError("Acute lesion mask must be 3D")
    mask_data = mask_image.get_fdata(dtype=np.float32)
    tolerance = float(cfg["binary_tolerance"])
    if not np.isfinite(mask_data).all() or not np.all((np.abs(mask_data) <= tolerance) | (np.abs(mask_data - 1) <= tolerance)):
        raise ValueError("Acute lesion mask must contain only finite binary 0/1 values")
    lesion = mask_data > 0.5
    lesion_n = int(lesion.sum())
    voxel_ml = float(abs(np.linalg.det(mask_image.affine[:3, :3])) / 1000)
    row: dict[str, Any] = {key: str(subject.get(key, "")) for key in ("subject_key", "participant_id", "cohort")}
    row["acute_volume_ml"] = lesion_n * voxel_ml
    row["acute_empty_mask"] = lesion_n == 0
    union_support = np.zeros(lesion.shape, dtype=bool)
    atlas_records: list[dict[str, Any]] = []
    qc: list[dict[str, Any]] = []
    used_columns: set[str] = set()
    for atlas in atlases:
        if atlas["template"] != declared:
            raise ValueError(f"Atlas {atlas['name']} uses a different declared template")
        atlas_path = Path(atlas["path"])
        digest = _sha256(atlas_path)
        if atlas.get("sha256") and str(atlas["sha256"]).lower() != digest:
            raise ValueError(f"Atlas SHA256 mismatch: {atlas['name']}")
        label_path = Path(atlas["labels"])
        labels = _labels(label_path)
        image = nib.load(str(atlas_path))
        _geometry(image, str(atlas["name"]), expected_unit)
        raw = image.get_fdata(dtype=np.float32)
        if not np.isfinite(raw).all():
            raise ValueError(f"Atlas contains nonfinite values: {atlas['name']}")
        include = {int(v) for v in atlas["include_labels"]}
        exclude = {int(v) for v in atlas["exclude_labels"]}
        selected = sorted((include if include else set(labels)) - exclude)
        if not set(selected).issubset(labels):
            raise ValueError(f"Configured labels absent from dictionary: {atlas['name']}")
        atlas_name = re.sub(r"[^A-Za-z0-9_]+", "_", str(atlas["name"]))
        if atlas["kind"] == "labels":
            if raw.ndim != 3 or not np.allclose(raw, np.rint(raw), rtol=0, atol=tolerance):
                raise ValueError(f"Deterministic atlas must have 3D integer labels: {atlas['name']}")
            actual = {int(v) for v in np.unique(raw)} - {0}
            if not actual.issubset(labels):
                raise ValueError(f"Label dictionary does not cover atlas values: {atlas['name']}")
            sampled = np.rint(_on_mask_grid(image, mask_image, "nearest", bool(cfg["allow_grid_resampling"]), float(cfg["grid_affine_tolerance"]))).astype(np.int32)
            weights = ((index, (sampled == index).astype(np.float32)) for index in selected)
        elif atlas["kind"] == "probability":
            scale = float(cfg["probability_scale"])
            if scale <= 0 or raw.ndim != 4 or np.min(raw) < -tolerance or np.max(raw) > scale + tolerance:
                raise ValueError(f"Probability atlas values/dimensions outside configured range: {atlas['name']}")
            if any(index < 0 or index >= raw.shape[3] for index in selected):
                raise ValueError("Probability atlas dictionary indices must be zero-based volume indices")
            sampled = np.clip(_on_mask_grid(image, mask_image, "linear", bool(cfg["allow_grid_resampling"]), float(cfg["grid_affine_tolerance"])) / scale, 0, 1)
            weights = ((index, sampled[..., index]) for index in selected)
        else:
            raise ValueError(f"Unknown atlas kind: {atlas['kind']}")
        for index, weight in weights:
            column = f"acute_{atlas_name}_{index}"
            if column in used_columns:
                raise ValueError(f"Repeated acute feature column: {column}")
            used_columns.add(column)
            denominator = float(weight.sum(dtype=np.float64))
            numerator = float(weight[lesion].sum(dtype=np.float64))
            row[column] = numerator / denominator if denominator > 0 else None
            union_support |= weight > 0
            qc.append({"feature": column, "atlas": atlas["name"], "index": index, "name": labels[index],
                       "atlas_weight_sum": denominator, "lesion_weight_sum": numerator,
                       "status": "valid" if denominator > 0 else "empty_region_on_grid"})
        atlas_records.append({"name": atlas["name"], "path": str(atlas_path), "sha256": digest,
                              "labels_sha256": _sha256(label_path), "template": atlas["template"],
                              "kind": atlas["kind"], "license": atlas["license"], "source": atlas["source"],
                              "selected_indices": selected})
    row["acute_uncovered_fraction"] = float(np.count_nonzero(lesion & ~union_support) / lesion_n) if lesion_n else None
    _write_tsv(output_dir / "acute.tsv", [row], list(row))
    qc_fields = ["feature", "atlas", "index", "name", "atlas_weight_sum", "lesion_weight_sum", "status"]
    _write_tsv(output_dir / "region_qc.tsv", qc, qc_fields)
    manifest = {"lesion_path": str(mask_path), "lesion_sha256": _sha256(mask_path),
                "lesion_template": declared, "voxel_volume_ml": voxel_ml,
                "spatial_units": mask_image.header.get_xyzt_units()[0], "atlases": atlas_records}
    (output_dir / "atlas_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return {"status": "completed", "reason": "Empty lesion mask; exclude from acute location modelling" if not lesion_n else "",
            "artifacts": ["acute.tsv", "region_qc.tsv", "atlas_manifest.json"], "n_subjects": 1,
            "metrics": {"acute_volume_ml": row["acute_volume_ml"], "region_count": len(used_columns),
                        "uncovered_fraction": row["acute_uncovered_fraction"]}}
