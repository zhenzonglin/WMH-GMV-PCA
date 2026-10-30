"""建立固定DAG：每个区域、CV重复、bootstrap和患者批次均有独立回执。"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .config import digest
from .ingest import collect_cohort, sha256_file


def build_manifest(config: dict, run_dir: str | Path) -> dict:
    run_dir = Path(run_dir).resolve()
    inputs_dir = run_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    frames, features, evidence, clinical_snapshots = {}, {}, {}, {}
    for cohort in ("LS", "C3", "GB"):
        frame, local_features, sources = collect_cohort(config, cohort)
        if features and local_features != features:
            raise ValueError(f"Region columns/order differ between LS and {cohort}")
        features = local_features
        frames[cohort] = frame
        frame.to_csv(inputs_dir / f"{cohort}.tsv", sep="\t", index=False)
        for source in sources:
            observed = frame.attrs.get("source_hashes", {}).get(str(source))
            current = sha256_file(source)
            if observed and observed != current:
                raise ValueError(f"Input changed during snapshot preparation: {source}")
            evidence[str(source)] = observed or current
        clinical_path = config["paths"]["cohorts"][cohort]["clinical_table"]
        clinical_snapshots[cohort] = None
        if clinical_path and Path(clinical_path).is_file():
            content = Path(clinical_path).read_bytes()
            digest_now = hashlib.sha256(content).hexdigest()
            if evidence.get(clinical_path) != digest_now:
                raise ValueError(f"Clinical input changed during snapshot preparation: {clinical_path}")
            destination = inputs_dir / f"{cohort}_clinical{Path(clinical_path).suffix}"
            destination.write_bytes(content)
            clinical_snapshots[cohort] = str(destination)
    for atlas in config["paths"]["atlases"]:
        for field in ("path", "labels"):
            path = Path(atlas[field])
            if path.is_file():
                evidence[str(path)] = sha256_file(path)
    root = Path(config["project"]["project_root"])
    code = {}
    for folder, pattern in (("src", "*.py"), ("R", "*.R"), ("workflow", "*")):
        for path in sorted((root / folder).rglob(pattern)):
            if path.is_file() and "__pycache__" not in str(path):
                code[str(path.relative_to(root))] = sha256_file(path)
    fingerprint = digest({"files": evidence, "implementation": code,
                          "declared_inputs": config["paths"]["cohorts"]})
    (inputs_dir / "source_manifest.json").write_text(json.dumps(
        {"input_fingerprint": fingerprint, "files": evidence, "implementation": code}, indent=2))
    (inputs_dir / "features.json").write_text(json.dumps(features, indent=2))
    tasks, by_id = [], {}

    def add(task_id, stage, action, cohort="ALL", deps=(), unit="cohort", subjects=None, **parameters):
        item = {"task_id": task_id, "stage": stage, "unit": unit, "cohort": cohort,
                "parameters": {"action": action, **parameters}, "inputs": {},
                "dependencies": list(dict.fromkeys(deps)),
                "output_dir": str(run_dir / "artifacts" / task_id),
                "n_subjects": len(frames[cohort]) if cohort in frames else None}
        if subjects is not None:
            item["subject_keys"] = subjects
            item["n_subjects"] = len(subjects)
        tasks.append(item)
        by_id[task_id] = item
        return task_id

    def output(task_id, filename):
        return str(Path(by_id[task_id]["output_dir"]) / filename)

    def patient_batches(cohort, stage, action, deps=(), **parameters):
        frame = frames[cohort]
        size = config["runtime"]["execution"]["patient_batch_size"]
        ids = []
        for start in range(0, len(frame), size):
            keys = frame.iloc[start:start + size]["subject_key"].tolist()
            local_parameters = dict(parameters)
            if "input_hashes" in local_parameters:
                relevant = set(frame.iloc[start:start + size]["lesion_mask"].tolist())
                relevant.update(a[field] for a in config["paths"]["atlases"] for field in ("path", "labels"))
                local_parameters["input_hashes"] = {k: v for k, v in local_parameters["input_hashes"].items() if k in relevant}
            tid = add(f"{cohort}.{stage}.{start // size:05d}", stage, action, cohort,
                      deps=deps, unit="patient_batch", subjects=keys,
                      source=str(inputs_dir / f"{cohort}.tsv"), features=features,
                      **local_parameters)
            ids.append(tid)
        return ids

    ingested = {}
    for cohort in frames:
        batch = patient_batches(cohort, "ingest", "ingest")
        ingested[cohort] = add(f"{cohort}.inputs", "input_matrix", "assemble_rows", cohort,
                              deps=batch, files=[output(t, "rows.tsv") for t in batch])
    if not config["stages"]["representation"]:
        return {"schema_version": 1, "input_fingerprint": fingerprint, "features": features,
                "cohorts": {c: {"n": len(f), "subject_keys": f.subject_key.tolist()} for c, f in frames.items()},
                "tasks": tasks, "target_task_ids": [t["task_id"] for t in tasks], "blocked_task_ids": []}
    matrices, score_tasks = {}, {}
    for modality in ("wmh", "gmv"):
        models = {}
        for index, feature in enumerate(features[modality]):
            tid = add(f"norm.{modality}.{index:03d}", f"norm_fit_{modality}", "r", "LS",
                      deps=[ingested["LS"]], unit="region", r_action="norm_fit", feature=feature,
                      modality=modality, r_inputs={"train": output(ingested["LS"], "data.tsv")})
            models[feature] = tid
        for cohort in frames:
            applied = []
            for index, feature in enumerate(features[modality]):
                model = models[feature]
                tid = add(f"{cohort}.norm.{modality}.{index:03d}", f"norm_apply_{modality}", "r", cohort,
                          deps=[ingested[cohort], model], unit="region", r_action="norm_apply",
                          feature=feature, modality=modality,
                          r_inputs={"data": output(ingested[cohort], "data.tsv"), "model": output(model, "model.rds")})
                applied.append(tid)
            matrix = add(f"{cohort}.{modality}.matrix", f"deviation_matrix_{modality}", "assemble_deviations", cohort,
                         deps=[ingested[cohort], *applied], base=output(ingested[cohort], "data.tsv"),
                         files=[output(t, "scores.tsv") for t in applied], modality=modality,
                         features=features[modality])
            matrices[cohort, modality] = matrix
            patient_batches(cohort, f"{modality}_deviation", "slice_patient", deps=[matrix],
                            table=output(matrix, "data.tsv"), required_suffix="_dev",
                            active_features=output(matrix, "active_features.json"))

        train = matrices["C3", modality]
        score_tasks[modality] = add_pca(add, output, config, modality, train,
                                      output(train, "data.tsv"), config["analysis"]["pca"]["max_components"])

    projected, joint_tasks = {}, {}
    for cohort in ("C3", "GB"):
        for modality in ("wmh", "gmv"):
            fitted = score_tasks[modality]
            matrix = matrices[cohort, modality]
            tid = add(f"{cohort}.{modality}.project", f"pca_project_{modality}", "r", cohort,
                      deps=[fitted, matrix], r_action="pca_project", modality=modality,
                      r_inputs={"data": output(matrix, "data.tsv"), "model": output(fitted, "model.rds")})
            projected[cohort, modality] = tid
            patient_batches(cohort, f"{modality}_spatial_scores", "slice_patient", deps=[tid],
                            table=output(tid, "scores.tsv"), required_prefix=f"{modality.upper()}_PC")
        deps = [projected[cohort, "wmh"], projected[cohort, "gmv"]]
        r_inputs = {"wmh_scores": output(deps[0], "scores.tsv"), "gmv_scores": output(deps[1], "scores.tsv")}
        if cohort == "GB":
            deps.append(joint_tasks["C3"])
            r_inputs["model"] = output(joint_tasks["C3"], "model.rds")
        joint = add(f"{cohort}.joint", "joint_phenotype", "r", cohort, deps=deps,
                    r_action="joint", r_inputs=r_inputs)
        joint_tasks[cohort] = joint
        patient_batches(cohort, "joint_scores", "slice_patient", deps=[joint],
                        table=output(joint, "scores.tsv"), required_prefix="JOINT_")

    acute_projected, acute_raw, acute_fit = {}, {}, None
    if config["stages"]["acute_location"]:
        for cohort in ("C3", "GB"):
            batch = patient_batches(cohort, "acute_features", "acute_extract", deps=[ingested[cohort]],
                                    input_hashes=evidence)
            acute_raw[cohort] = add(f"{cohort}.acute.matrix", "acute_matrix", "assemble_rows", cohort,
                                  deps=batch, files=[output(t, "rows.tsv") for t in batch])
        residual = add("C3.acute.residual", "acute_volume_adjustment", "r", "C3", deps=[acute_raw["C3"]],
                       r_action="acute_residualize", modality="acute",
                       r_inputs={"train": output(acute_raw["C3"], "data.tsv")})
        fitted = add_pca(add, output, config, "acute", residual, output(residual, "residuals.tsv"),
                         config["analysis"]["acute"]["max_components"])
        acute_fit = fitted
        for cohort in ("C3", "GB"):
            r = residual
            if cohort == "GB":
                r = add("GB.acute.residual", "acute_volume_adjustment", "r", "GB",
                        deps=[residual, acute_raw["GB"]], r_action="acute_residualize", modality="acute",
                        r_inputs={"data": output(acute_raw["GB"], "data.tsv"), "model": output(residual, "model.rds")})
            acute_projected[cohort] = add(f"{cohort}.acute.project", "acute_position", "r", cohort,
                deps=[r, fitted], r_action="pca_project", modality="acute",
                r_inputs={"data": output(r, "residuals.tsv"), "model": output(fitted, "model.rds")})

    clinical_models = {}
    for cohort in ("C3", "GB"):
        deps = [ingested[cohort], joint_tasks[cohort], projected[cohort, "wmh"], projected[cohort, "gmv"]]
        inputs = [output(ingested[cohort], "data.tsv"), output(projected[cohort, "wmh"], "scores.tsv"),
                  output(projected[cohort, "gmv"], "scores.tsv"), output(joint_tasks[cohort], "scores.tsv")]
        for modality in ("wmh", "gmv"):
            deps.append(matrices[cohort, modality])
            inputs.append(output(matrices[cohort, modality], "data.tsv"))
        for mapping in (acute_raw, acute_projected):
            if cohort in mapping:
                deps.append(mapping[cohort])
                inputs.append(output(mapping[cohort], "data.tsv" if mapping is acute_raw else "scores.tsv"))
        prep = add(f"{cohort}.clinical.prepare", "clinical_variables", "clinical_prepare", cohort,
                   deps=deps, files=inputs, clinical_table=clinical_snapshots[cohort])
        if cohort == "C3" and config["stages"]["clinical_b1"] and config["analysis"]["clinical"]["enable_internal_validation"]:
            add("C3.internal_validation", "internal_prediction_validation", "internal_validation", "C3",
                deps=[prep, matrices["C3", "wmh"], matrices["C3", "gmv"]], unit="model",
                inputs={"wmh_deviations": output(matrices["C3", "wmh"], "data.tsv"),
                        "gmv_deviations": output(matrices["C3", "gmv"], "data.tsv"),
                        "clinical_wide": output(prep, "clinical_wide.tsv"), "run_dir": str(run_dir)})
        for action, enable in (("clinical_b1_single", "clinical_b1"),
                               ("clinical_b1_longitudinal", "clinical_b1"), ("clinical_b2", "clinical_b2")):
            if not config["stages"][enable]:
                continue
            r_inputs = {"data": output(prep, "clinical_wide.tsv"), "long": output(prep, "clinical_long.tsv")}
            modeldeps = [prep]
            if action == "clinical_b2":
                for modality in ("wmh", "gmv"):
                    r_inputs[f"{modality}_component_metadata"] = output(score_tasks[modality], "component_metadata.tsv")
                    modeldeps.append(score_tasks[modality])
                if acute_fit:
                    r_inputs["acute_component_metadata"] = output(acute_fit, "component_metadata.tsv")
                    modeldeps.append(acute_fit)
            if cohort == "GB":
                reference = clinical_models["C3", action]
                modeldeps.append(reference)
                r_inputs["reference_models"] = output(reference, "models.rds")
            clinical_models[cohort, action] = add(f"{cohort}.{action}", action, "r", cohort,
                deps=modeldeps, unit="model", r_action=action, r_inputs=r_inputs)
        if config["stages"]["export_statistician"]:
            patient_batches(cohort, "handoff", "slice_patient", deps=[prep], table=output(prep, "clinical_wide.tsv"))
            add(f"{cohort}.statistician_package", "statistician_package", "handoff_package", cohort,
                deps=[prep], source_dir=by_id[prep]["output_dir"])
    return {"schema_version": 1, "input_fingerprint": fingerprint, "features": features,
            "cohorts": {c: {"n": len(f), "subject_keys": f.subject_key.tolist()} for c, f in frames.items()},
            "tasks": tasks, "target_task_ids": [t["task_id"] for t in tasks], "blocked_task_ids": []}


def add_pca(add, output, config, modality, input_task, matrix, maximum):
    reference = add(f"{modality}.pca.reference", f"pca_reference_{modality}", "r", "C3", deps=[input_task],
                    r_action="pca_reference", modality=modality, k="max", r_inputs={"train": matrix})
    cvs, bootstraps = [], []
    for rep in range(1, config["analysis"]["pca"]["cv_repeats"] + 1):
        cvs.append(add(f"{modality}.cv.{rep:04d}", f"pca_cv_{modality}", "r", "C3", deps=[input_task],
            unit="repeat", r_action="pca_cv", modality=modality, replicate=rep, r_inputs={"train": matrix}))
    for k in range(1, maximum + 1):
        for rep in range(1, config["analysis"]["pca"]["bootstrap_repeats"] + 1):
            bootstraps.append(add(f"{modality}.bootstrap.k{k}.{rep:04d}", f"bootstrap_{modality}", "r", "C3",
                deps=[reference], unit="repeat", r_action="pca_bootstrap", modality=modality, k=k, replicate=rep,
                r_inputs={"train": matrix, "reference": output(reference, "model.rds")}))
    select = add(f"{modality}.pca.select", f"pca_selection_{modality}", "r", "C3", deps=[*cvs, *bootstraps],
        r_action="pca_select", modality=modality,
        r_inputs={"cv": [output(t, "criterion.tsv") for t in cvs],
                  "bootstrap": [output(t, "bootstrap.tsv") for t in bootstraps]})
    return add(f"{modality}.pca.fit", f"pca_fit_{modality}", "r", "C3", deps=[input_task, select],
               r_action="pca_fit", modality=modality,
               r_inputs={"train": matrix, "selection": output(select, "selection.json")})
