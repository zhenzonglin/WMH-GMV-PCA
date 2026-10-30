from __future__ import annotations

import copy
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from wmh_gmv.config import load_config
from wmh_gmv.stages import execute_task, r_request

ROOT = Path(__file__).resolve().parents[1]


def test_unestimable_upstream_is_not_reported_as_missing_patient_data(tmp_path):
    cfg = load_config(ROOT / "config/project.yaml")
    parent = tmp_path / "model"
    parent.mkdir()
    (parent / "task_result.json").write_text(json.dumps({"status": "not_estimable", "reason": "No stable acute component"}))
    result = r_request({"action": "pca_project", "inputs": {"model": str(parent / "model.rds")}}, cfg, tmp_path / "projection")
    assert result["status"] == "not_estimable"
    assert "No stable acute component" in result["reason"]


def test_unknown_config_key_rejected(tmp_path):
    for path in (ROOT / "config").glob("*.yaml"):
        (tmp_path / path.name).write_bytes(path.read_bytes())
    runtime = yaml.safe_load((tmp_path / "runtime.yaml").read_text())
    runtime["execution"]["misspelled_cores"] = 999
    (tmp_path / "runtime.yaml").write_text(yaml.safe_dump(runtime))
    with pytest.raises(ValueError, match="misspelled_cores"):
        load_config(tmp_path / "project.yaml")


def test_stage_dependencies_are_explicit(tmp_path):
    for path in (ROOT / "config").glob("*.yaml"):
        (tmp_path / path.name).write_bytes(path.read_bytes())
    project = yaml.safe_load((tmp_path / "project.yaml").read_text())
    project["stages"]["representation"] = False
    (tmp_path / "project.yaml").write_text(yaml.safe_dump(project))
    with pytest.raises(ValueError, match="Dependent analysis"):
        load_config(tmp_path / "project.yaml")


def base_task(tmp_path, action):
    source = tmp_path / "source.tsv"
    pd.DataFrame([{"subject_key": "C3:0001", "participant_id": "0001", "cohort": "C3",
                   "good_dev": 1.0, "unestimable_dev": None}]).to_csv(source, sep="\t", index=False)
    return {"task_id": "test.batch", "unit": "patient_batch", "cohort": "C3", "subject_keys": ["C3:0001"],
            "output_dir": str(tmp_path / "published"), "parameters": {"action": action, "source": str(source)}}


def test_unestimable_roi_does_not_invalidate_every_patient(tmp_path):
    cfg = load_config(ROOT / "config/project.yaml")
    task = base_task(tmp_path, "slice_patient")
    active = tmp_path / "active.json"
    active.write_text(json.dumps({"features": ["good_dev"]}))
    task["parameters"].update(required_suffix="_dev", active_features=str(active))
    result = execute_task(task, cfg, tmp_path / "stage")
    state = json.loads((tmp_path / "stage/subjects.json").read_text())
    assert result["status"] == "completed"
    assert state["C3:0001"]["status"] == "completed"


def test_source_change_is_not_silently_consumed(tmp_path):
    cfg = load_config(ROOT / "config/project.yaml")
    task = base_task(tmp_path, "acute_extract")
    task["parameters"]["input_hashes"] = {task["parameters"]["source"]: "wrong-checksum"}
    with pytest.raises(ValueError, match="input_changed_since_run_snapshot"):
        execute_task(task, cfg, tmp_path / "stage")


def test_batch_retry_keeps_successful_rows(tmp_path, monkeypatch):
    from wmh_gmv import acute
    cfg = load_config(ROOT / "config/project.yaml")
    cfg["runtime"]["execution"]["retry_failed_on_resume"] = True
    task = base_task(tmp_path, "acute_extract")
    source = Path(task["parameters"]["source"])
    rows = pd.read_csv(source, sep="\t", dtype=str)
    rows = pd.concat([rows, rows.assign(subject_key="C3:0002", participant_id="0002")])
    rows.to_csv(source, sep="\t", index=False)
    task["subject_keys"].append("C3:0002")
    previous = Path(task["output_dir"])
    previous.mkdir()
    (previous / "subjects.json").write_text(json.dumps({"C3:0001": {"status": "completed"},
                                                       "C3:0002": {"status": "failed"}}))
    pd.DataFrame([{"subject_key": "C3:0001", "participant_id": "0001", "cohort": "C3",
                   "acute_volume_ml": 123.0}]).to_csv(previous / "rows.tsv", sep="\t", index=False)
    called = []
    def fake_extract(row, config, out):
        called.append(row["subject_key"])
        pd.DataFrame([{**row, "acute_volume_ml": 456.0}]).to_csv(out / "acute.tsv", sep="\t", index=False)
        return {"status": "completed", "artifacts": ["acute.tsv"]}
    monkeypatch.setattr(acute, "extract_acute", fake_extract)
    execute_task(task, cfg, tmp_path / "stage")
    assert called == ["C3:0002"]
    actual = pd.read_csv(tmp_path / "stage/rows.tsv", sep="\t")
    assert actual.loc[actual.subject_key == "C3:0001", "acute_volume_ml"].item() == 123.0
