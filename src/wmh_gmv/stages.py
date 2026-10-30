"""科学步骤适配层。核心统计计算委托固定的R包实现。"""
from __future__ import annotations

import json
import os
import subprocess
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd


def _write_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    os.replace(temporary, path)


def table(path: str | Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t", dtype={"subject_key": str, "participant_id": str, "cohort": str})
    if "subject_key" not in data or data.subject_key.isna().any() or data.subject_key.duplicated().any():
        raise ValueError(f"Missing/duplicated subject_key in {path}")
    return data


def join_tables(frames: list[pd.DataFrame]) -> pd.DataFrame:
    if not frames:
        return pd.DataFrame(columns=["subject_key", "participant_id", "cohort"])
    out = frames[0]
    for frame in frames[1:]:
        new = [c for c in frame if c not in out or c == "subject_key"]
        out = out.merge(frame[new], on="subject_key", how="left", validate="one_to_one")
    return out


def r_request(request: dict, config: dict, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    request = {**request, "config": config, "output_dir": str(output_dir)}
    action = request["action"]
    inputs = request.get("inputs", {})
    for name, value in inputs.items():
        values = value if isinstance(value, list) else [value]
        if name in {"cv", "bootstrap"}:
            request["inputs"][name] = [p for p in values if Path(p).is_file()]
            continue
        if value and not Path(value).is_file():
            upstream = Path(value).parent / "task_result.json"
            if upstream.is_file():
                previous_result = json.loads(upstream.read_text())
                if previous_result.get("status") == "not_estimable":
                    return {"status": "not_estimable", "reason":
                            f"Upstream scientific stage not estimable ({name}): {previous_result.get('reason', '')}",
                            "artifacts": []}
            return {"status": "waiting_input", "reason": f"Required upstream artifact unavailable: {name}", "artifacts": []}
    if action.startswith("pca_") and action not in {"pca_project", "pca_select"}:
        data = table(inputs["train"])
        suffix = "_resid" if request.get("modality") == "acute" else "_dev"
        request["features"] = [c for c in data if c.endswith(suffix) and pd.to_numeric(data[c], errors="coerce").notna().any()]
        if not request["features"] and action in {"pca_cv", "pca_bootstrap"}:
            # K0仍可拟合；没有区域时无需伪造遮蔽预测误差。
            if action == "pca_cv":
                pd.DataFrame([{"k": 0, "msep": 0, "replicate": request.get("replicate", 1),
                               "internal_repeats": 1, "masked_cells": 0, "zero_features": True}]).to_csv(
                    output_dir / "criterion.tsv", sep="\t", index=False)
                return {"status": "completed", "reason": "No estimable features: mean-only K0 comparator",
                        "artifacts": ["criterion.tsv"], "n_subjects": len(data)}
            return {"status": "not_estimable", "reason": "No estimable features", "artifacts": []}
    if action == "acute_residualize" and "train" in inputs:
        data = table(inputs["train"])
        age = pd.to_numeric(data.get("age", pd.Series(np.nan, index=data.index)), errors="coerce")
        lo, hi = config["analysis"]["age_range"]
        training = data.loc[age.between(lo, hi)]
        exclude = {"acute_volume_ml", "acute_empty_mask", "acute_uncovered_fraction"}
        minimum = max(config["analysis"]["acute"]["minimum_region_affected"],
                      int(np.ceil(len(training) * config["analysis"]["acute"]["minimum_region_fraction"])))
        request["features"] = [c for c in data if c.startswith("acute_") and c not in exclude
                               and (pd.to_numeric(training[c], errors="coerce") > 0).sum() >= minimum]
    _write_json(output_dir / "request.json", request)
    env = dict(os.environ)
    env["R_LIBS"] = config["paths"]["r_library"] + os.pathsep + env.get("R_LIBS", "")
    threads = str(config["runtime"]["execution"]["blas_threads_per_job"])
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[name] = threads
    with (output_dir / "R.log").open("w", encoding="utf-8") as log:
        process = subprocess.run([config["paths"]["rscript"],
            str(Path(config["project"]["project_root"]) / "R/run_stage.R"), str(output_dir / "request.json")],
            stdout=log, stderr=subprocess.STDOUT, env=env, check=False)
    if process.returncode:
        tail = (output_dir / "R.log").read_text(encoding="utf-8", errors="replace")[-3500:]
        raise RuntimeError(f"R stage {action} exited {process.returncode}:\n{tail}")
    result = json.loads((output_dir / "result.json").read_text())
    artifacts = result.get("artifacts", [])
    if isinstance(artifacts, str):
        artifacts = [artifacts]
    result["artifacts"] = [*artifacts, "request.json", "R.log", "result.json"]
    return result


def _patient_batch(task: dict, config: dict, output: Path) -> dict:
    from .ingest import normalize_subject
    params = task["parameters"]
    action = params["action"]
    keys = task["subject_keys"]
    states = {key: {"status": "pending", "reason": ""} for key in keys}
    previous = Path(task["output_dir"])
    retained_rows = []
    retry_keys = set(keys)
    if config["runtime"]["execution"]["retry_failed_on_resume"] and (previous / "subjects.json").is_file():
        old_states = json.loads((previous / "subjects.json").read_text())
        if set(old_states) == set(keys) and any(s["status"] == "failed" for s in old_states.values()):
            states = old_states
            retry_keys = {key for key in keys if states[key]["status"] == "failed"}
            for child in previous.iterdir():
                if child.name == "task_result.json":
                    continue
                if child.is_dir():
                    shutil.copytree(child, output / child.name, dirs_exist_ok=True)
                else:
                    shutil.copy2(child, output / child.name)
            old_rows = table(previous / "rows.tsv")
            retained_rows = old_rows.loc[~old_rows.subject_key.isin(retry_keys)].to_dict(orient="records")
    _write_json(output / "subjects.json", states)
    rows, audit = retained_rows, []
    def verify_sources():
        from .ingest import sha256_file
        for filename, expected in params.get("input_hashes", {}).items():
            if not Path(filename).is_file() or sha256_file(filename) != expected:
                raise ValueError(f"input_changed_since_run_snapshot: {filename}")
    verify_sources()
    source_path = params.get("table") or params.get("source")
    source = table(source_path) if source_path and Path(source_path).is_file() else None
    if source is not None:
        source = source.set_index("subject_key", drop=False)
    for index, key in enumerate(keys):
        if key not in retry_keys:
            continue
        states[key] = {"status": "running", "reason": ""}
        _write_json(output / "subjects.json", states)
        try:
            if source is None or key not in source.index:
                states[key] = {"status": "waiting_input", "reason": "Upstream subject row unavailable"}
                continue
            row = source.loc[key].where(pd.notna(source.loc[key]), None).to_dict()
            if action == "ingest":
                normalized, issue = normalize_subject(row, params["features"], config)
                rows.append(normalized)
                audit.append(issue)
                states[key] = {"status": "completed", "reason": "; ".join(issue["issues"])}
            elif action == "acute_extract":
                from .acute import extract_acute
                directory = output / f"subject_{index:05d}"
                directory.mkdir(exist_ok=True)
                result = extract_acute(row, config, directory)
                if result["status"] == "completed":
                    acute_row = table(directory / "acute.tsv").iloc[0].to_dict()
                    acute_row["age"] = row.get("age")
                    rows.append(acute_row)
                states[key] = {"status": result["status"], "reason": result.get("reason", "")}
            else:
                required = [c for c in row if
                            (params.get("required_suffix") and c.endswith(params["required_suffix"])) or
                            (params.get("required_prefix") and c.startswith(params["required_prefix"]))]
                if params.get("active_features"):
                    required = json.loads(Path(params["active_features"]).read_text())["features"]
                if (params.get("required_prefix") or params.get("required_suffix")) and not required:
                    states[key] = {"status": "not_estimable", "reason": "No retained estimable component/region"}
                elif required and not all(pd.notna(row[c]) for c in required):
                    states[key] = {"status": "not_estimable", "reason": "Unavailable subject projection"}
                else:
                    states[key] = {"status": "completed", "reason": ""}
                rows.append(row)
        except (ValueError, OSError, KeyError) as exc:
            states[key] = {"status": "failed", "reason": str(exc)}
        finally:
            # 人数只在对应行产物实际落盘后推进。
            pd.DataFrame(rows if rows else [], columns=None if rows else ["subject_key", "participant_id", "cohort"]).to_csv(
                output / "rows.tsv", sep="\t", index=False)
            _write_json(output / "subjects.json", states)
    _write_json(output / "audit.json", audit)
    verify_sources()
    artifacts = [str(p.relative_to(output)) for p in output.rglob("*") if p.is_file()]
    return {"status": "completed", "artifacts": artifacts, "n_subjects": len(keys),
            "metrics": {"failed_subjects": sum(s["status"] == "failed" for s in states.values())}}


def execute_task(task: dict, config: dict, staging_dir: str | Path) -> dict:
    output = Path(staging_dir)
    output.mkdir(parents=True, exist_ok=True)
    params = task["parameters"]
    action = params["action"]
    if task["unit"] == "patient_batch":
        return _patient_batch(task, config, output)
    if action == "assemble_rows":
        frames = [table(path) for path in params["files"] if Path(path).is_file()]
        data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["subject_key", "participant_id", "cohort"])
        if data.subject_key.duplicated().any():
            raise ValueError("Duplicate patient across batches")
        data.to_csv(output / "data.tsv", sep="\t", index=False)
        return {"status": "completed", "artifacts": ["data.tsv"], "n_subjects": len(data)}
    if action == "assemble_deviations":
        available = [table(p) for p in params["files"] if Path(p).is_file()]
        active = [c for frame in available for c in frame if c.endswith("_dev")]
        data = join_tables([table(params["base"]), *available])
        for feature in params["features"]:
            if feature + "_dev" not in data:
                data[feature + "_dev"] = np.nan
        data.to_csv(output / "data.tsv", sep="\t", index=False)
        _write_json(output / "active_features.json", {"features": active})
        return {"status": "completed", "artifacts": ["data.tsv", "active_features.json"], "n_subjects": len(data)}
    if action == "r":
        request = {k: v for k, v in params.items() if k not in {"action", "r_action", "r_inputs"}}
        request.update(action=params["r_action"], inputs=params["r_inputs"], cohort=task["cohort"])
        return r_request(request, config, output)
    if action == "clinical_prepare":
        from .clinical import prepare_clinical
        data = join_tables([table(p) for p in params["files"] if Path(p).is_file()])
        data.to_csv(output / "imaging.tsv", sep="\t", index=False)
        clinical = params.get("clinical_table")
        result = prepare_clinical(output / "imaging.tsv", Path(clinical) if clinical else None,
                                  config, output, task["cohort"])
        result["artifacts"] = [*result.get("artifacts", []), "imaging.tsv"]
        if clinical and Path(clinical).is_file():
            filename = "clinical_source" + Path(clinical).suffix
            shutil.copy2(clinical, output / filename)
            result["artifacts"].append(filename)
        return result
    if action == "internal_validation":
        from .validation import run_internal_validation
        return run_internal_validation(config, params["inputs"], output)
    if action == "handoff_package":
        from .handoff import package_handoff
        return package_handoff(Path(params["source_dir"]), output, task["cohort"], config)
    if action == "report":
        from .report import create_report
        return create_report(Path(params["run_dir"]), output, config)
    if action == "test_sleep" and config["project"]["mode"] == "synthetic":
        time.sleep(params["seconds"])
        (output / "result.txt").write_text("synthetic lifecycle test\n")
        return {"status": "completed", "artifacts": ["result.txt"]}
    raise ValueError(f"Unsupported task action: {action}")
