"""Lifecycle contracts and an actual short Snakemake drain/resume integration."""
from __future__ import annotations

import importlib.metadata
import os
import shutil
import sys
import time
import types
from pathlib import Path

import pytest

from wmh_gmv import lifecycle, monitor
from wmh_gmv.worker import run_task
from wmh_gmv.state import (artifact_receipts, atomic_json, check_process,
                           process_identity, read_json, task_dir, verify_artifacts)


def minimal_config(tmp_path: Path) -> dict:
    return {
        "project": {"project_root": str(tmp_path / "project"), "mode": "synthetic"},
        "paths": {"output_root": str(tmp_path / "outputs"), "state_root": str(tmp_path / "state"),
                  "log_root": str(tmp_path / "logs"), "temp_root": str(tmp_path / "temp"),
                  "python": sys.executable, "rscript": shutil.which("Rscript") or sys.executable},
        "runtime": {
            "execution": {"cores": 1, "max_parallel_jobs": 1, "blas_threads_per_job": 1,
                          "memory_mb": 32, "background": True, "keep_going": True,
                          "automatic_retries": 0, "retry_failed_on_resume": False, "startup_timeout_seconds": 10},
            "monitor": {"heartbeat_seconds": 0.1, "stale_heartbeat_seconds": 2, "refresh_seconds": 0.1,
                        "outputs": ["json", "tsv", "html"]},
            "shutdown": {"mode": "drain", "force_kill": False, "wait_for_exit_seconds": 0.2,
                         "drain_timeout_seconds": None},
            "resume": {"reuse_verified_outputs": True, "require_matching_scientific_config": True},
            "environment": {"scheduler_version": "9.27.1"}},
        "_config_path": str(tmp_path / "config.yaml"), "_config_hash": "config", "_scientific_hash": "science",
    }


def make_run(config: dict) -> Path:
    root = Path(config["paths"]["output_root"])
    run = root / "test-run"
    run.mkdir(parents=True)
    atomic_json(root / ".active-run.json", {"run_dir": str(run)})
    atomic_json(run / "run.json", {"run_id": run.name, "status": "completed",
                                  "config_path": config["_config_path"], "scientific_hash": config["_scientific_hash"],
                                  "project_root": config["project"]["project_root"]})
    return run


def test_process_identity_rejects_reused_pid() -> None:
    identity = process_identity(os.getpid())
    assert check_process(identity)[0] == "alive"
    altered = {**identity, "start_ticks": "wrong"}
    assert check_process(altered)[0] == "mismatch"


def test_artifact_hash_detects_same_size_mutation(tmp_path: Path) -> None:
    file = tmp_path / "data.tsv"
    file.write_text("one")
    receipt = {"artifacts": artifact_receipts(tmp_path, ["data.tsv"])}
    assert verify_artifacts(receipt)[0]
    file.write_text("two")
    assert not verify_artifacts(receipt)[0]


def test_missing_manifest_is_unknown_not_zero(tmp_path: Path) -> None:
    config = minimal_config(tmp_path)
    make_run(config)
    value = monitor.snapshot(config)
    assert value["status"] == "unknown"
    assert value["total_people"] is None


def test_monitor_current_attempt_excludes_old_failure(tmp_path: Path) -> None:
    config = minimal_config(tmp_path)
    run = make_run(config)
    task = {"task_id": "person1", "stage": "projection", "cohort": "C3", "unit": "patient", "subject_key": "C3:001"}
    atomic_json(run / "manifest.json", {"tasks": [task]})
    location = task_dir(run, "person1")
    atomic_json(location / "state.json", {"status": "running", "attempt_id": "second", "process": process_identity(os.getpid()), "history_failures": 1})
    from wmh_gmv.state import utc_now
    atomic_json(location / "heartbeat.json", {"time": utc_now(), "attempt_id": "second"})
    row = monitor.snapshot(config)["patient_steps"][0]
    assert row["running"] == 1 and row["failed"] == 0
    assert row["historical_failures"] == 1 and row["count_invariant_ok"]


def test_patient_batch_progress_is_not_model_count(tmp_path: Path) -> None:
    config = minimal_config(tmp_path)
    run = make_run(config)
    staging = tmp_path / "staging"
    staging.mkdir()
    task = {"task_id": "batch", "stage": "norm_apply", "cohort": "C3", "unit": "patient_batch",
            "subject_keys": ["C3:001", "C3:002", "C3:003"], "output_dir": str(run / "artifacts" / "batch")}
    atomic_json(run / "manifest.json", {"tasks": [task]})
    atomic_json(task_dir(run, "batch") / "state.json", {"status": "running", "process": process_identity(os.getpid()),
                                                        "staging_dir": str(staging), "attempt_id": "1"})
    atomic_json(staging / "subjects.json", {"C3:001": {"status": "completed"}, "C3:002": {"status": "failed"}})
    value = monitor.snapshot(config)
    row = value["patient_steps"][0]
    assert (row["completed"], row["failed"], row["pending"], row["total"]) == (1, 1, 1, 3)
    assert value["group_steps"] == []


def test_stop_rejects_identity_mismatch_without_request(tmp_path: Path) -> None:
    config = minimal_config(tmp_path)
    run = make_run(config)
    meta = read_json(run / "run.json")
    meta.update(status="running", scheduler={**process_identity(os.getpid()), "start_ticks": "wrong"})
    atomic_json(run / "run.json", meta)
    with pytest.raises(RuntimeError, match="identity mismatch"):
        lifecycle.stop_project(config)
    assert not (run / "stop_request.json").exists()


def test_resume_blocks_historical_failure_and_downstream(tmp_path: Path) -> None:
    config = minimal_config(tmp_path)
    run = make_run(config)
    tasks = [{"task_id": name, "output_dir": str(run / "artifacts" / name), "dependencies": deps}
             for name, deps in [("a", []), ("b", ["a"]), ("c", [])]]
    atomic_json(task_dir(run, "a") / "state.json", {"status": "failed", "reason": "known failure"})
    result = lifecycle._prepare_resume(run, {"tasks": tasks}, config)
    assert result["blocked_task_ids"] == ["a", "b"]
    assert result["target_task_ids"] == ["c"]


def _completed_task(run: Path, task: dict, failed_subjects: int = 0) -> None:
    output = Path(task["output_dir"])
    output.mkdir(parents=True)
    atomic_json(output / "subjects.json", {"C3:001": {"status": "completed"}})
    (output / "rows.tsv").write_text("subject_key\nC3:001\n")
    atomic_json(task_dir(run, task["task_id"]) / "state.json", {
        "status": "completed", "scientific_hash": "science", "attempt_id": "previous",
        "metrics": {"failed_subjects": failed_subjects},
        "artifacts": artifact_receipts(output, ["subjects.json", "rows.tsv"]),
    })
    (task_dir(run, task["task_id"]) / "finished.flag").write_text("completed\n")


def test_partial_batch_retry_invalidates_downstream_only_when_enabled(tmp_path: Path) -> None:
    config = minimal_config(tmp_path)
    run = make_run(config)
    tasks = [
        {"task_id": "batch", "unit": "patient_batch", "dependencies": [], "output_dir": str(run / "artifacts" / "batch")},
        {"task_id": "model", "unit": "model", "dependencies": ["batch"], "output_dir": str(run / "artifacts" / "model")},
    ]
    _completed_task(run, tasks[0], failed_subjects=1)
    _completed_task(run, tasks[1])
    lifecycle._prepare_resume(run, {"tasks": tasks}, config)
    assert all((task_dir(run, task["task_id"]) / "finished.flag").exists() for task in tasks)
    assert lifecycle.completion_status(run, tasks, 0) == "completed_with_failures"
    config["runtime"]["execution"]["retry_failed_on_resume"] = True
    result = lifecycle._prepare_resume(run, {"tasks": tasks}, config)
    assert result["blocked_task_ids"] == []
    assert all(not (task_dir(run, task["task_id"]) / "finished.flag").exists() for task in tasks)
    assert all(read_json(task_dir(run, task["task_id"]) / "state.json")["status"] == "pending" for task in tasks)
    assert all(Path(task["output_dir"]).exists() for task in tasks)


@pytest.mark.parametrize("retry", [False, True])
def test_worker_partial_batch_cache_respects_retry_setting(tmp_path: Path, monkeypatch, retry: bool) -> None:
    config = minimal_config(tmp_path)
    config["runtime"]["execution"]["retry_failed_on_resume"] = retry
    run = make_run(config)
    task = {"task_id": "batch", "unit": "patient_batch", "dependencies": [], "output_dir": str(run / "artifacts" / "batch")}
    _completed_task(run, task, failed_subjects=1)
    atomic_json(run / "manifest.json", {"tasks": [task]})
    atomic_json(run / "resolved_config.json", config)
    calls = []

    def execute(_task, _config, staging):
        calls.append(True)
        # 科学层在旧发布目录读取成功患者，worker必须保留该目录直至新结果校验完成。
        assert (Path(task["output_dir"]) / "rows.tsv").exists()
        atomic_json(staging / "subjects.json", {"C3:001": {"status": "completed"}})
        (staging / "rows.tsv").write_text("subject_key\nC3:001\n")
        return {"status": "completed", "metrics": {"failed_subjects": 0}, "artifacts": ["subjects.json", "rows.tsv"]}

    monkeypatch.setitem(sys.modules, "wmh_gmv.stages", types.SimpleNamespace(execute_task=execute))
    assert run_task(run, "batch") == 0
    assert bool(calls) == retry
    current = read_json(task_dir(run, "batch") / "state.json")
    assert (current["attempt_id"] != "previous") == retry
    assert lifecycle.completion_status(run, [task], 0) == ("completed" if retry else "completed_with_failures")


def test_failed_patient_recovery_recomputes_descendant_but_preserves_unrelated_cache(tmp_path: Path, monkeypatch) -> None:
    config = minimal_config(tmp_path)
    config["runtime"]["execution"]["retry_failed_on_resume"] = True
    run = make_run(config)
    tasks = [
        {"task_id": "source", "unit": "patient_batch", "dependencies": [], "output_dir": str(run / "artifacts" / "source")},
        {"task_id": "aggregate", "unit": "model", "dependencies": ["source"], "output_dir": str(run / "artifacts" / "aggregate")},
        {"task_id": "unrelated", "unit": "model", "dependencies": [], "output_dir": str(run / "artifacts" / "unrelated")},
    ]
    for task in tasks:
        _completed_task(run, task)
    source_state = read_json(task_dir(run, "source") / "state.json")
    atomic_json(task_dir(run, "source") / "state.json", {**source_state, "status": "failed"})
    unrelated_time = (task_dir(run, "unrelated") / "finished.flag").stat().st_mtime_ns
    manifest = lifecycle._prepare_resume(run, {"tasks": tasks}, config)
    atomic_json(run / "manifest.json", manifest)
    atomic_json(run / "resolved_config.json", config)
    calls = []

    def execute(task, _config, staging):
        calls.append(task["task_id"])
        if task["task_id"] == "source":
            (staging / "rows.tsv").write_text("subject_key\nC3:001\nC3:002\n")
            atomic_json(staging / "subjects.json", {"C3:001": {"status": "completed"}, "C3:002": {"status": "completed"}})
            return {"status": "completed", "n_subjects": 2, "metrics": {"failed_subjects": 0}, "artifacts": ["rows.tsv", "subjects.json"]}
        count = len((Path(tasks[0]["output_dir"]) / "rows.tsv").read_text().splitlines()) - 1
        (staging / "model_count.txt").write_text(str(count))
        return {"status": "completed", "n_subjects": count, "artifacts": ["model_count.txt"]}

    monkeypatch.setitem(sys.modules, "wmh_gmv.stages", types.SimpleNamespace(execute_task=execute))
    assert run_task(run, "source") == 0
    assert run_task(run, "aggregate") == 0
    assert run_task(run, "unrelated") == 0
    assert calls == ["source", "aggregate"]
    assert read_json(task_dir(run, "aggregate") / "state.json")["n_subjects"] == 2
    assert read_json(task_dir(run, "aggregate") / "state.json")["attempt_id"] != "previous"
    assert read_json(task_dir(run, "unrelated") / "state.json")["attempt_id"] == "previous"
    assert (task_dir(run, "unrelated") / "finished.flag").stat().st_mtime_ns == unrelated_time


def test_failed_run_still_generates_report_after_terminal_status(tmp_path: Path, monkeypatch) -> None:
    config = minimal_config(tmp_path)
    run = make_run(config)
    meta = read_json(run / "run.json")
    meta["status"] = "completed_with_failures"
    atomic_json(run / "run.json", meta)
    calls = []

    def create_report(run_dir, output, _config):
        calls.append(read_json(run_dir / "run.json")["status"])
        assert output.is_dir()
        (output / "report.html").write_text("<p>failed and completed tasks are both reported</p>")
        return {"status": "completed", "artifacts": ["report.html"]}

    monkeypatch.setitem(sys.modules, "wmh_gmv.report", types.SimpleNamespace(create_report=create_report))
    result = lifecycle.finalize_report(run, config)
    assert calls == ["completed_with_failures"]
    assert result["status"] == "completed_with_failures" and result["report_status"] == "completed"
    assert Path(result["report_path"]).is_file()


def test_report_error_is_recorded_without_rewriting_scientific_status(tmp_path: Path, monkeypatch) -> None:
    config = minimal_config(tmp_path)
    run = make_run(config)

    def create_report(*_args):
        raise RuntimeError("synthetic report error")

    monkeypatch.setitem(sys.modules, "wmh_gmv.report", types.SimpleNamespace(create_report=create_report))
    result = lifecycle.finalize_report(run, config)
    assert result["status"] == "completed" and result["report_status"] == "failed"
    assert "synthetic report error" in result["report_error"]
    assert Path(result["report_error_log"]).is_file()


def test_html_escapes_error_content() -> None:
    value = {"status": "unknown", "updated_at": "now", "patient_steps": [], "group_steps": [],
             "errors": ["<script>alert(1)</script>"]}
    rendered = monitor.html_report(value)
    assert "<script>" not in rendered
    assert "&lt;script&gt;" in rendered


def test_cooperative_checkpoint_reuses_verified_staging(tmp_path: Path, monkeypatch) -> None:
    config = minimal_config(tmp_path)
    run = make_run(config)
    task = {"task_id": "cv", "output_dir": str(run / "artifacts" / "cv"), "stage": "cv",
            "cohort": "C3", "unit": "model", "dependencies": []}
    atomic_json(run / "manifest.json", {"tasks": [task]})
    atomic_json(run / "resolved_config.json", config)

    def interrupted(_task, _config, staging):
        (staging / "checkpoint.txt").write_text("completed fold 1")
        atomic_json(run / "stop_request.json", {"requested_at": "now"})
        return {"status": "interrupted", "reason": "Stopped after fold 1", "artifacts": ["checkpoint.txt"]}

    monkeypatch.setitem(sys.modules, "wmh_gmv.stages", types.SimpleNamespace(execute_task=interrupted))
    assert run_task(run, "cv") == 75
    first = read_json(task_dir(run, "cv") / "state.json")
    assert first["status"] == "pending" and first["resume_available"]
    assert not (task_dir(run, "cv") / "finished.flag").exists()
    (run / "stop_request.json").unlink()

    def resumed(_task, _config, staging):
        assert str(staging) == first["staging_dir"]
        assert (staging / "checkpoint.txt").read_text() == "completed fold 1"
        (staging / "done.txt").write_text("fold 2 complete")
        return {"status": "completed", "artifacts": ["checkpoint.txt", "done.txt"]}

    monkeypatch.setitem(sys.modules, "wmh_gmv.stages", types.SimpleNamespace(execute_task=resumed))
    assert run_task(run, "cv") == 0
    assert read_json(task_dir(run, "cv") / "state.json")["checkpoint_reused"]


def _wait_until(predicate, timeout: float = 30) -> None:
    start = time.monotonic()
    while time.monotonic() - start < timeout:
        if predicate():
            return
        time.sleep(0.1)
    raise AssertionError("Timed out waiting for lifecycle condition")


def test_real_scheduler_drain_and_resume(tmp_path: Path, monkeypatch) -> None:
    try:
        version = importlib.metadata.version("snakemake")
    except importlib.metadata.PackageNotFoundError:
        pytest.skip("Real scheduler integration requires pinned Snakemake")
    if version != "9.27.1":
        pytest.skip("Only verified scheduler version is tested")
    config = minimal_config(tmp_path)
    root = Path(config["project"]["project_root"])
    package = root / "src" / "wmh_gmv"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    source = Path(lifecycle.__file__).parent
    for name in ("lifecycle.py", "state.py", "worker.py"):
        shutil.copyfile(source / name, package / name)
    (package / "stages.py").write_text(
        "import time\n"
        "def execute_task(task, config, staging_dir):\n"
        "    time.sleep(task['parameters']['seconds'])\n"
        "    (staging_dir / 'done.txt').write_text(task['task_id'])\n"
        "    return {'status':'completed','artifacts':['done.txt']}\n"
    )
    workflow = root / "workflow"
    workflow.mkdir()
    (workflow / "Snakefile").write_text('''import json
import sys
from pathlib import Path
run = config["run_dir"]
manifest = json.load(open(config["manifest"]))
tasks = {t["task_id"]: t for t in manifest["tasks"]}
python = sys.executable
rule all:
    input: [f"{run}/tasks/{x}/finished.flag" for x in manifest["target_task_ids"]]
rule work:
    input: lambda w: [f"{run}/tasks/{x}/finished.flag" for x in tasks[w.taskid]["dependencies"]]
    output: f"{run}/tasks/{{taskid}}/finished.flag"
    resources: mem_mb=1
    shell: "{python:q} -m wmh_gmv.worker --run-dir {run:q} --task-id {wildcards.taskid:q}"
''')

    def build_manifest(_config, run):
        return {"input_fingerprint": "unchanged", "tasks": [
            {"task_id": "a", "unit": "patient", "subject_key": "C3:001", "stage": "a", "cohort": "C3", "dependencies": [],
             "output_dir": str(run / "artifacts" / "a"), "parameters": {"seconds": 4}},
            {"task_id": "b", "unit": "patient", "subject_key": "C3:001", "stage": "b", "cohort": "C3", "dependencies": ["a"],
             "output_dir": str(run / "artifacts" / "b"), "parameters": {"seconds": 0.1}}]}
    monkeypatch.setitem(sys.modules, "wmh_gmv.pipeline", types.SimpleNamespace(build_manifest=build_manifest))
    first = lifecycle.start_project(config)
    run = Path(first["run_dir"])
    try:
        _wait_until(lambda: (task_dir(run, "a") / "state.json").exists() and read_json(task_dir(run, "a") / "state.json")["status"] == "running")
        assert lifecycle.start_project(config)["status"] == "already_running"
        stopped = lifecycle.stop_project(config)
        assert stopped["status"] in ("draining", "stopped")
        _wait_until(lambda: read_json(run / "run.json")["status"] == "stopped" and not lifecycle._is_live(run))
        assert (task_dir(run, "a") / "finished.flag").exists()
        assert not (task_dir(run, "b") / "finished.flag").exists()
        receipt = read_json(task_dir(run, "a") / "state.json")
        flag_time = (task_dir(run, "a") / "finished.flag").stat().st_mtime_ns
        resumed = lifecycle.start_project(config)
        assert resumed["resumed"] and resumed["run_dir"] == str(run)
        _wait_until(lambda: read_json(run / "run.json")["status"] == "completed" and not lifecycle._is_live(run))
        assert read_json(task_dir(run, "a") / "state.json")["attempt_id"] == receipt["attempt_id"]
        assert (task_dir(run, "a") / "finished.flag").stat().st_mtime_ns == flag_time
        assert verify_artifacts(receipt)[0]
        assert (task_dir(run, "b") / "finished.flag").exists()
    finally:
        if lifecycle._is_live(run):
            lifecycle.stop_project(config)
            _wait_until(lambda: not lifecycle._is_live(run))
