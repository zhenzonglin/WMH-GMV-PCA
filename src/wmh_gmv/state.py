"""原子状态、产物校验及 Linux 进程身份。监测不改变任务状态。"""
from __future__ import annotations

import hashlib
import json
import os
import re
import socket
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: str | Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def read_json(path: str | Path) -> Any:
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_task_id(task_id: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", task_id) or task_id in (".", ".."):
        raise ValueError(f"Invalid task_id: {task_id!r}")
    return task_id


def within(path: str | Path, root: str | Path) -> Path:
    result, base = Path(path).resolve(), Path(root).resolve()
    if result == base or base not in result.parents:
        raise ValueError(f"Path must be a child of {base}: {result}")
    return result


def task_dir(run_dir: str | Path, task_id: str) -> Path:
    return Path(run_dir) / "tasks" / safe_task_id(task_id)


def process_identity(pid: int) -> dict[str, Any]:
    """不能只用 PID 判断身份；start_ticks 和 boot_id 防止 PID 复用。"""
    base = Path("/proc") / str(int(pid))
    stat = (base / "stat").read_text()
    fields = stat[stat.rfind(")") + 2:].split()
    if fields[0] == "Z":
        raise ProcessLookupError(f"PID {pid} is a zombie")
    command = (base / "cmdline").read_bytes().decode("utf-8", "replace").rstrip("\0").split("\0")
    return {
        "pid": int(pid), "pgid": int(fields[2]), "start_ticks": fields[19],
        "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
        "hostname": socket.gethostname(), "cwd": os.readlink(base / "cwd"),
        "argv": command, "uid": base.stat().st_uid,
    }


def check_process(record: dict[str, Any] | None) -> tuple[str, str]:
    if not record:
        return "unknown", "Process identity is not registered"
    try:
        live = process_identity(int(record["pid"]))
    except (FileNotFoundError, ProcessLookupError):
        return "gone", "Registered process has exited"
    except (OSError, KeyError, ValueError) as exc:
        return "unknown", str(exc)
    for key in ("pid", "pgid", "start_ticks", "boot_id", "hostname", "cwd", "argv", "uid"):
        if record.get(key) != live.get(key):
            return "mismatch", f"Process identity mismatch: {key}"
    return "alive", ""


def handles_sigterm(pid: int) -> bool:
    try:
        for line in (Path("/proc") / str(pid) / "status").read_text().splitlines():
            if line.startswith("SigCgt:"):
                return bool(int(line.split()[1], 16) & (1 << 14))
    except OSError:
        return False
    return False


def artifact_receipts(output_dir: Path, names: list[str]) -> list[dict[str, Any]]:
    receipts = []
    for name in names:
        if (output_dir / name).is_symlink():
            raise ValueError(f"Declared artifact cannot be a symlink: {name}")
        path = within(output_dir / name, output_dir)
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Declared artifact is not a regular file: {path}")
        info = path.stat()
        receipts.append({"path": str(path), "sha256": sha256_file(path),
                         "size": info.st_size, "mtime_ns": info.st_mtime_ns})
    return receipts


def verify_artifacts(state: dict[str, Any], full_hash: bool = True) -> tuple[bool, str]:
    if not state.get("artifacts"):
        return False, "No artifact receipts"
    try:
        for item in state["artifacts"]:
            path = Path(item["path"])
            info = path.stat()
            if not path.is_file() or path.is_symlink() or info.st_size != item["size"]:
                return False, f"Artifact missing or size/type changed: {path}"
            if full_hash or info.st_mtime_ns != item.get("mtime_ns"):
                if sha256_file(path) != item["sha256"]:
                    return False, f"Artifact checksum mismatch: {path}"
    except (OSError, KeyError, TypeError) as exc:
        return False, str(exc)
    return True, ""


def has_failed_patients(task: dict[str, Any], state: dict[str, Any]) -> bool:
    """批次程序成功退出不代表批次内所有患者成功。"""
    return (task.get("unit") == "patient_batch"
            and int(state.get("metrics", {}).get("failed_subjects", 0) or 0) > 0)


def retry_partial_batch(task: dict[str, Any], state: dict[str, Any], config: dict[str, Any]) -> bool:
    return (state.get("status") == "completed" and has_failed_patients(task, state)
            and config["runtime"]["execution"]["retry_failed_on_resume"])


def append_event(path: Path, event: dict[str, Any]) -> None:
    """每个任务由一个进程写入独立事件文件，避免共享追加竞争。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({"time": utc_now(), **event}, ensure_ascii=False, allow_nan=False) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def inspect_task(run_dir: Path, task: dict[str, Any], full_hash: bool = False) -> dict[str, Any]:
    location = task_dir(run_dir, task["task_id"])
    state_path = location / "state.json"
    if not state_path.exists():
        status = "unknown" if (location / "events.jsonl").exists() else "pending"
        return {"status": status, "reason": "State missing" if status == "unknown" else "Not started"}
    try:
        state = read_json(state_path)
        if not isinstance(state, dict):
            raise ValueError("Task state is not an object")
    except (OSError, ValueError) as exc:
        return {"status": "unknown", "reason": str(exc)}
    status = state.get("status", "unknown")
    if status == "running":
        alive, reason = check_process(state.get("process"))
        if alive == "gone":
            return {**state, "status": "pending", "reason": "Interrupted; resume required"}
        if alive != "alive":
            return {**state, "status": "unknown", "reason": reason}
    if status in ("completed", "not_estimable", "waiting_input"):
        valid, reason = verify_artifacts(state, full_hash)
        if not valid:
            return {**state, "status": "unknown", "reason": reason}
    return state
