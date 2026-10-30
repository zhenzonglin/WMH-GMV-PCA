"""单个可恢复任务：调用科学模块，验证产物后原子发布。"""
from __future__ import annotations

import argparse
import os
import threading
import traceback
import uuid
from pathlib import Path
from typing import Any

from .state import (append_event, artifact_receipts, atomic_json, process_identity,
                    read_json, task_dir, utc_now, verify_artifacts, within, retry_partial_batch)


def run_task(run_dir: Path, task_id: str) -> int:
    run_dir = run_dir.resolve()
    config = read_json(run_dir / "resolved_config.json")
    manifest = read_json(run_dir / "manifest.json")
    matches = [item for item in manifest["tasks"] if item["task_id"] == task_id]
    if len(matches) != 1:
        raise ValueError(f"Task not found or duplicated: {task_id}")
    task = matches[0]
    location = task_dir(run_dir, task_id)
    location.mkdir(parents=True, exist_ok=True)
    old: dict[str, Any] = {}
    if (location / "state.json").exists():
        old = read_json(location / "state.json")
        if old.get("status") in ("completed", "not_estimable", "waiting_input"):
            valid, _ = verify_artifacts(old)
            if valid and old.get("scientific_hash") == config["_scientific_hash"] and not retry_partial_batch(task, old, config):
                if not (location / "finished.flag").exists():
                    (location / "finished.flag").write_text(old.get("status", "completed") + "\n")
                return 0
        if (location / "finished.flag").exists():
            (location / "finished.flag").unlink()
    if (run_dir / "stop_request.json").exists():
        # 调度器与启动边界上的竞争也不应启动新的科学计算。
        atomic_json(location / "state.json", {"status": "pending", "reason": "Stop requested", "task_id": task_id})
        return 75
    attempt_id = uuid.uuid4().hex
    attempt = location / "attempts" / attempt_id
    attempt.mkdir(parents=True)
    output = within(task["output_dir"], run_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging_root = Path(config["paths"]["temp_root"]) / run_dir.name
    staging_root.mkdir(parents=True, exist_ok=True)
    if staging_root.stat().st_dev != output.parent.stat().st_dev:
        raise ValueError("temp_root and output_root must share a filesystem for atomic publication")
    staging = staging_root / f"{task_id}.{attempt_id}"
    checkpoint_reused = False
    if old.get("resume_available") and old.get("scientific_hash") == config["_scientific_hash"]:
        previous = within(old["staging_dir"], staging_root)
        valid, _ = verify_artifacts({"artifacts": old.get("checkpoint_artifacts", [])})
        if previous.is_dir() and valid:
            staging = previous
            checkpoint_reused = True
    if not checkpoint_reused:
        staging.mkdir()
    base = {"task_id": task_id, "attempt_id": attempt_id, "started_at": utc_now(),
            "scientific_hash": config["_scientific_hash"], "process": process_identity(os.getpid()),
            "status": "running", "history_failures": old.get("history_failures", 0),
            "staging_dir": str(staging), "checkpoint_reused": checkpoint_reused}
    atomic_json(location / "state.json", base)
    append_event(location / "events.jsonl", {"event": "started", "attempt_id": attempt_id})
    done = threading.Event()

    def heartbeat() -> None:
        interval = config["runtime"]["monitor"]["heartbeat_seconds"]
        while not done.is_set():
            atomic_json(location / "heartbeat.json", {"time": utc_now(), "attempt_id": attempt_id})
            done.wait(interval)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        from .stages import execute_task
        result = execute_task(task, config, staging)
        status = result.get("status")
        if status == "interrupted":
            if not (run_dir / "stop_request.json").exists():
                raise ValueError("Cooperative interruption requires a recorded stop request")
            checkpoints = artifact_receipts(staging, result.get("artifacts", []))
            final = {**base, "status": "pending", "reason": result.get("reason", "Stopped at checkpoint"),
                     "finished_at": utc_now(), "resume_available": bool(checkpoints), "checkpoint_artifacts": checkpoints}
            atomic_json(location / "state.json", final)
            atomic_json(attempt / "receipt.json", final)
            append_event(location / "events.jsonl", {"event": "checkpoint_stopped", "attempt_id": attempt_id})
            return 75
        if status not in ("completed", "not_estimable", "waiting_input"):
            raise ValueError(f"Unsupported scientific result status: {status}")
        names = result.get("artifacts", [])
        if status == "completed" and not names:
            raise ValueError("Completed task declares no validated artifacts")
        atomic_json(staging / "task_result.json", result)
        names = sorted(set([*names, "task_result.json"]))
        artifact_receipts(staging, names)  # 发布之前校验路径、类型及所有文件。
        if output.exists():
            archive = within(output.parent / f".{output.name}.superseded.{attempt_id}", run_dir)
            os.replace(output, archive)
        os.replace(staging, output)
        receipts = artifact_receipts(output, names)
        final = {**base, "status": status, "reason": result.get("reason", ""),
                 "finished_at": utc_now(), "n_subjects": result.get("n_subjects", task.get("n_subjects")),
                 "artifacts": receipts, "metrics": result.get("metrics", {})}
        atomic_json(location / "state.json", final)
        atomic_json(attempt / "receipt.json", final)
        append_event(location / "events.jsonl", {"event": status, "attempt_id": attempt_id})
        (location / "finished.flag").write_text(status + "\n", encoding="utf-8")
        return 0
    except Exception as exc:
        error = traceback.format_exc()
        (attempt / "error.log").write_text(error, encoding="utf-8")
        final = {**base, "status": "failed", "reason": str(exc), "finished_at": utc_now(),
                 "history_failures": base["history_failures"] + 1, "error_log": str(attempt / "error.log")}
        atomic_json(location / "state.json", final)
        atomic_json(attempt / "receipt.json", final)
        append_event(location / "events.jsonl", {"event": "failed", "attempt_id": attempt_id, "reason": str(exc)})
        flag = location / "finished.flag"
        if flag.exists():
            flag.unlink()
        print(error, flush=True)
        return 1
    finally:
        done.set()
        thread.join(timeout=2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--task-id", required=True)
    args = parser.parse_args()
    try:
        return run_task(args.run_dir, args.task_id)
    except Exception as exc:
        # 配置、临时目录等启动期错误也要成为明确失败，不能永远显示待完成。
        location = task_dir(args.run_dir, args.task_id)
        attempt = location / "attempts" / uuid.uuid4().hex
        attempt.mkdir(parents=True, exist_ok=True)
        error = traceback.format_exc()
        (attempt / "startup_error.log").write_text(error, encoding="utf-8")
        if (location / "state.json").exists():
            (attempt / "previous_state.json").write_bytes((location / "state.json").read_bytes())
        atomic_json(location / "state.json", {"task_id": args.task_id, "status": "failed", "reason": str(exc),
                                              "finished_at": utc_now(), "error_log": str(attempt / "startup_error.log"),
                                              "history_failures": 1})
        print(error, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
