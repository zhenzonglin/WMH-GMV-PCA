#!/usr/bin/env python3
"""Generate explicitly synthetic, config-driven inputs for the real backends."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np
import yaml


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def make(config_path: Path) -> Path:
    spec = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if spec["scenario"] not in {"degenerate_v1", "spatial_contrast_v2"}:
        raise ValueError("Unknown explicitly configured synthetic scenario")
    root = (config_path.parent / spec["project_root"]).resolve()
    output = (root / spec["output_dir"]).resolve()
    config_dir = (root / spec["config_dir"]).resolve()
    if output.exists() and any(output.iterdir()) and not spec["overwrite"]:
        raise FileExistsError("Synthetic output exists; change output_dir or set overwrite in generation config")
    output.mkdir(parents=True, exist_ok=True); config_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(spec["seed"])
    shape = tuple(spec["atlas_shape"]); voxel = spec["voxel_size_mm"]
    atlas = np.zeros(shape, dtype=np.int16)
    for j, indices in enumerate(np.array_split(np.arange(shape[0]), 3), 1):
        atlas[indices, :, :] = j
    affine = np.diag([voxel, voxel, voxel, 1.0])
    atlas_path = output / "synthetic_atlas.nii.gz"
    atlas_image = nib.Nifti1Image(atlas, affine)
    atlas_image.header.set_xyzt_units("mm")
    nib.save(atlas_image, atlas_path)
    atlas_labels = output / "synthetic_atlas.tsv"
    write_tsv(atlas_labels, [{"index": i, "name": f"Synthetic_region_{i}"} for i in (1, 2, 3)])
    cohort_paths = {}
    for cohort, n in spec["cohort_sizes"].items():
        directory = output / cohort; directory.mkdir(exist_ok=True)
        masks = directory / "masks"; masks.mkdir(exist_ok=True)
        roster, metadata, wmh_rows, gmv_rows, clinical = [], [], [], [], []
        for index in range(n):
            identifier = f"{index + 1:06d}"
            age = rng.uniform(*spec["age_range"]); sex = int(rng.integers(2))
            icv = rng.normal(1450 + 70 * sex, 80)
            latent = rng.normal(size=2) + (0 if cohort == "LS" else 1.0)
            w = np.empty(spec["wmh_regions"])
            for j in range(len(w)):
                signal = latent[j % 2] * (0.18 if cohort == "LS" else 0.65)
                mean = np.exp(-0.1 + .035 * (age - 60) + .15 * sex + signal)
                w[j] = 0 if rng.uniform() < .20 else rng.gamma(5, mean / 5)
            g = np.empty(spec["gmv_regions"])
            weights = np.linspace(.8, 1.2, len(g)); weights /= weights.sum()
            # 两个空间对比与独立全局因子分开，避免总GM由前两个PC近乎确定。
            angles = 2 * np.pi * np.arange(len(g)) / len(g)
            patterns = np.column_stack([np.cos(angles), np.sin(angles)])
            patterns -= np.sum(patterns * weights[:, None], axis=0)
            pattern_sd = np.sqrt(np.sum(patterns ** 2 * weights[:, None], axis=0))
            patterns /= np.where(pattern_sd > 0, pattern_sd, 1)
            global_gmv = rng.normal(0, spec["gmv_global_sd"]) if spec["scenario"] == "spatial_contrast_v2" else 0
            for j in range(len(g)):
                spatial = float(patterns[j] @ latent) if spec["scenario"] == "spatial_contrast_v2" else latent[j % 2]
                shift = (0 if cohort == "LS" else spec["gmv_spatial_amplitude"] * spatial)
                g[j] = 620 * weights[j] * np.exp(-.005 * (age - 60) + .75 * np.log(icv / 1450) - shift + global_gmv + rng.normal(0, spec["gmv_regional_noise_sd"]))
            lesion_path = ""
            acute_ml = 0.0
            if cohort != "LS":
                lesion = np.zeros(shape, dtype=np.uint8)
                region = int(rng.integers(1, 4))
                selected = np.flatnonzero(atlas == region)
                count = int(rng.integers(5, max(6, len(selected) // 2)))
                lesion.flat[rng.choice(selected, count, replace=False)] = 1
                lesion_file = masks / f"{identifier}_lesion.nii.gz"
                lesion_image = nib.Nifti1Image(lesion, affine)
                lesion_image.header.set_xyzt_units("mm")
                nib.save(lesion_image, lesion_file)
                lesion_path = str(lesion_file)
                acute_ml = count * voxel ** 3 / 1000
            roster.append({"participant_id": identifier})
            metadata.append({"participant_id": identifier,"age":age,"sex":sex,"icv_ml":icv,
                "wmh_total_ml":float(w.sum()),"gmv_total_ml":float(g.sum()),
                "lesion_mask":lesion_path,"lesion_template":"SYNTHETIC"})
            wmh_rows.append({"participant_id":identifier,**{f"wmh_region_{j+1:02d}_ml":value for j,value in enumerate(w)}})
            gmv_rows.append({"participant_id":identifier,**{f"muse_{j+1:03d}_ml":value for j,value in enumerate(g)}})
            if cohort != "LS":
                base = {"participant_id":identifier,"prestroke_mrs":int(rng.integers(0,4)),
                    "nihss":int(rng.poisson(6)),"prior_stroke":int(rng.binomial(1,.2)),
                    "hypertension":int(rng.binomial(1,.6)),"diabetes":int(rng.binomial(1,.3)),
                    "af":int(rng.binomial(1,.2)),"ivt":int(rng.binomial(1,.3)),"evt":int(rng.binomial(1,.2)),
                    "index_date":"2025-01-01","death_date":""}
                person = rng.normal(0, spec["patient_random_intercept_sd"])
                dead = False
                for month, date, death_date in [(3,"2025-04-01","2025-03-01"),(6,"2025-07-01","2025-06-01"),(12,"2026-01-01","2025-12-01")]:
                    score = int(np.digitize(person + .35 * latent[0] + .2 * latent[1] + .4 * acute_ml + rng.logistic() - .02 * month,
                                             [-2,-1,-.3,.3,1,2]))
                    if dead or score == 6:
                        if not dead: base["death_date"] = death_date
                        score = 6; dead = True
                    base[f"mrs_m{month}"] = score
                    base[f"death_m{month}"] = int(dead)
                    base[f"visit_date_m{month}"] = date
                clinical.append(base)
        for filename, rows in [("roster.tsv",roster),("metadata.tsv",metadata),("wmh.tsv",wmh_rows),("gmv.tsv",gmv_rows)]:
            write_tsv(directory / filename, rows)
        if clinical: write_tsv(directory / "clinical.tsv", clinical)
        cohort_paths[cohort] = {"roster":str(directory / "roster.tsv"),"derivatives_root":str(directory),
            "wmh_table":str(directory / "wmh.tsv"),"gmv_table":str(directory / "gmv.tsv"),
            "metadata_table":str(directory / "metadata.tsv"),"clinical_table":str(directory / "clinical.tsv") if clinical else None,
            "raw_records":None}
    configs = {name:yaml.safe_load((root / "config" / f"{name}.yaml").read_text(encoding="utf-8"))
               for name in ("project","paths","analysis","variables","runtime")}
    configs["project"]["project"].update({"project_root":str(root),"mode":"synthetic"})
    configs["project"]["stages"]["clinical_b1"] = bool(spec["enable_clinical"])
    configs["project"]["stages"]["clinical_b2"] = bool(spec["enable_clinical"])
    paths = configs["paths"]
    paths.update({"cohorts":cohort_paths,"output_root":"runs/synthetic","state_root":"state/synthetic",
                  "log_root":"logs/synthetic","temp_root":"tmp/synthetic"})
    paths["atlases"] = [{"name":"Synthetic","kind":"labels","path":str(atlas_path),"labels":str(atlas_labels),
        "template":"SYNTHETIC","sha256":hashlib.sha256(atlas_path.read_bytes()).hexdigest(),
        "license":"CC0-synthetic","source":"synthetic-generator","include_labels":[],"exclude_labels":[0]}]
    analysis = configs["analysis"]
    analysis["seed"] = spec["seed"]
    analysis["normative"].update({"folds":spec["normative_folds"],"min_positive":spec["normative_min_positive"]})
    analysis["pca"].update({"max_components":spec["pca_max_components"],"cv_repeats":spec["pca_cv_repeats"],
        "bootstrap_repeats":spec["pca_bootstrap_repeats"]})
    analysis["acute"].update({"template":"SYNTHETIC","probability_scale":1})
    analysis["clinical"]["enable_internal_validation"] = bool(spec["enable_internal_validation"])
    variables = configs["variables"]
    variables["imaging"].update({"expected_wmh_regions":spec["wmh_regions"],"expected_gmv_regions":spec["gmv_regions"]})
    runtime = configs["runtime"]
    runtime["execution"].update({"cores":spec["cores"],"max_parallel_jobs":spec["parallel_jobs"],"memory_mb":spec["memory_mb"]})
    for name,value in configs.items():
        (config_dir / f"{name}.yaml").write_text(yaml.safe_dump(value,allow_unicode=True,sort_keys=False),encoding="utf-8")
    (output / "SYNTHETIC_ONLY.txt").write_text(f"Generated test fixtures: {spec['scenario']}. No real patient data. Reduced resampling counts are for smoke testing only.\n",encoding="utf-8")
    return config_dir / "project.yaml"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    print(make(args.config.resolve()))
