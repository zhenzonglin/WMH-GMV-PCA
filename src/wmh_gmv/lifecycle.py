"""配置驱动的启动、验证后续跑和安全排空；不接管其他项目进程。"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import signal
import shlex
import subprocess
import sys
import time
import traceback
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .state import (append_event, artifact_receipts, atomic_json, check_process, handles_sigterm, has_failed_patients,
                    inspect_task, process_identity, read_json, safe_task_id,
                    task_dir, utc_now, verify_artifacts, within, retry_partial_batch)


def _root(config: dict[str, Any]) -> Path:
    return Path(config["project"]["project_root"])


def _pointer(config: dict[str, Any]) -> Path:
    return Path(config["paths"]["output_root"]) / ".active-run.json"


def _active_run(config: dict[str, Any]) -> Path | None:
    pointer = _pointer(config)
    if not pointer.exists():
        return None
    value = read_json(pointer)
    return within(value["run_dir"], config["paths"]["output_root"])


def _environment(config: dict[str, Any]) -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(_root(config) / "src") + os.pathsep + env.get("PYTHONPATH", "")
    threads = str(config["runtime"]["execution"]["blas_threads_per_job"])
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = threads
    # 实际R入口将项目库前置；保留R默认库，以免切断已安装的传递依赖。
    if config["paths"].get("r_library"):
        env["R_LIBS"] = config["paths"]["r_library"] + os.pathsep + env.get("R_LIBS", "")
    env["PYTHONUNBUFFERED"] = "1"
    return env


def check_resources(config: dict[str, Any]) -> None:
    execution = config["runtime"]["execution"]
    available_cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    cpu_limit = Path("/sys/fs/cgroup/cpu.max")
    if cpu_limit.exists():
        quota, period = cpu_limit.read_text().split()
        if quota != "max":
            available_cpus = min(available_cpus or 1, max(1, int(quota) // int(period)))
    if execution["cores"] > (available_cpus or 1):
        raise RuntimeError(f"Configured cores={execution['cores']} exceeds available={available_cpus}")
    available_kb = None
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            available_kb = int(line.split()[1])
    cgroup_max, cgroup_current = Path("/sys/fs/cgroup/memory.max"), Path("/sys/fs/cgroup/memory.current")
    if cgroup_max.exists() and cgroup_current.exists() and cgroup_max.read_text().strip() != "max":
        free_kb = max(0, int(cgroup_max.read_text()) - int(cgroup_current.read_text())) // 1024
        available_kb = min(available_kb or free_kb, free_kb)
    if available_kb is not None and execution["memory_mb"] * 1024 > available_kb:
        raise RuntimeError(f"Configured memory_mb={execution['memory_mb']} exceeds available={available_kb // 1024}")


def _check_dependencies(config: dict[str, Any]) -> str:
    executable = Path(config["paths"]["python"])
    if not executable.is_file():
        raise FileNotFoundError(f"Configured Python executable unavailable: {executable}")
    check = subprocess.run([str(executable), "-c", "import importlib.metadata; print(importlib.metadata.version('snakemake'))"],
                           capture_output=True, text=True, check=True, env=_environment(config))
    actual = check.stdout.strip()
    expected = config["runtime"]["environment"]["scheduler_version"]
    if actual != expected:
        raise RuntimeError(f"Snakemake version {actual} differs from configured verified version {expected}")
    if not Path(config["paths"]["rscript"]).is_file():
        raise FileNotFoundError("Configured Rscript executable unavailable")
    return actual


def _code_fingerprint(config: dict[str, Any]) -> str:
    digest = hashlib.sha256()
    root = _root(config)
    files = []
    for directory in ("src", "R", "workflow", "scripts"):
        files.extend(path for path in (root / directory).rglob("*")
                     if path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc")
    files.extend(root / name for name in ("pyproject.toml", "uv.lock", "renv.lock") if (root / name).is_file())
    for path in sorted(files):
        digest.update(str(path.relative_to(root)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _software_inventory(config: dict[str, Any]) -> dict[str, Any]:
    """实际环境而非要求文件：库升级后旧计算不能悄然被当作同版本缓存。"""
    code = ("import json, sys, importlib.metadata as m; "
            "print(json.dumps({'python':sys.version, 'packages':"
            "{d.metadata['Name']:d.version for d in m.distributions() if d.metadata['Name']}}))")
    python = subprocess.run([config["paths"]["python"], "-c", code], check=True,
                            capture_output=True, text=True, env=_environment(config))
    result = json.loads(python.stdout)
    packages = config["runtime"]["environment"].get("r_packages", [])
    rcode = ("a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); p <- a[-1]; "
             "v <- vapply(p, function(x) if(requireNamespace(x,quietly=TRUE)) as.character(packageVersion(x)) else NA_character_, character(1)); "
             "d <- loadedNamespaces(); dv <- setNames(lapply(d,function(x) as.character(packageVersion(x))),d); "
             "cat(jsonlite::toJSON(list(R=R.version.string,packages=as.list(v),dependencies=dv),auto_unbox=TRUE,na='null'))")
    if packages:
        r = subprocess.run([config["paths"]["rscript"], "-e", rcode, config["paths"]["r_library"], *packages], check=True,
                           capture_output=True, text=True, env=_environment(config))
        result["r"] = json.loads(r.stdout)
        missing = [name for name, value in result["r"]["packages"].items() if value is None]
        if missing:
            raise RuntimeError("Required R packages unavailable: " + ", ".join(missing))
    return result


def _is_live(run_dir: Path) -> bool:
    try:
        meta = read_json(run_dir / "run.json")
    except (OSError, ValueError):
        return False
    return any(check_process(meta.get(key))[0] == "alive" for key in ("supervisor", "scheduler"))


def _prepare_resume(run_dir: Path, manifest: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    failed = set()
    invalidated = set()
    tasks = manifest["tasks"]
    for task in tasks:
        safe_task_id(task["task_id"])
        within(task["output_dir"], run_dir)
        location = task_dir(run_dir, task["task_id"])
        location.mkdir(parents=True, exist_ok=True)
        state_path = location / "state.json"
        flag = location / "finished.flag"
        if state_path.exists():
            # 不可读回执不能被自动当作待完成处理。
            old = read_json(state_path)
            if not isinstance(old, dict):
                raise ValueError(f"Invalid task receipt requires reconciliation: {state_path}")
            valid = old.get("scientific_hash") == config["_scientific_hash"] and verify_artifacts(old)[0]
            partial_retry = retry_partial_batch(task, old, config)
            if old.get("status") in ("completed", "not_estimable", "waiting_input") and valid and not partial_retry:
                if not flag.exists():
                    flag.write_text(old["status"] + "\n", encoding="utf-8")
                continue
            if old.get("status") == "failed" and not config["runtime"]["execution"]["retry_failed_on_resume"]:
                failed.add(task["task_id"])
            if old.get("status") == "running":
                live, why = check_process(old.get("process"))
                if live == "alive":
                    raise RuntimeError(f"Cannot resume while task {task['task_id']} is still active")
                if live in ("unknown", "mismatch"):
                    raise RuntimeError(f"Task identity requires reconciliation: {task['task_id']}: {why}")
            if flag.exists():
                flag.unlink()
            invalidated.add(task["task_id"])
            if task["task_id"] not in failed:
                reason = "Retry failed patients in completed batch" if partial_retry else "Resume required; cached receipt not valid"
                atomic_json(state_path, {**old, "status": "pending", "reason": reason})
        else:
            if flag.exists():
                flag.unlink()
            atomic_json(state_path, {"status": "pending", "reason": "Not started", "task_id": task["task_id"]})
            invalidated.add(task["task_id"])
    # 新结果必须向下游传播；不能让worker用旧completed回执短路新依赖。
    changed = True
    while changed:
        changed = False
        for task in tasks:
            if task["task_id"] not in invalidated and set(task.get("dependencies", [])) & invalidated:
                invalidated.add(task["task_id"])
                location = task_dir(run_dir, task["task_id"])
                old = read_json(location / "state.json")
                atomic_json(location / "state.json", {**old, "status": "pending", "resume_available": False,
                                                       "reason": "Upstream result will be recomputed; cached result invalidated"})
                flag = location / "finished.flag"
                if flag.exists():
                    flag.unlink()
                append_event(location / "events.jsonl", {"event": "dependency_cache_invalidated"})
                changed = True
    for task in tasks:
        if set(task.get("dependencies", [])) & invalidated:
            location = task_dir(run_dir, task["task_id"])
            old = read_json(location / "state.json")
            if old.get("resume_available"):
                atomic_json(location / "state.json", {**old, "resume_available": False,
                                                       "reason": "Upstream changed; previous checkpoint is invalid"})
    blocked = set(failed)
    changed = True
    while changed:
        changed = False
        for task in tasks:
            if task["task_id"] not in blocked and set(task.get("dependencies", [])) & blocked:
                blocked.add(task["task_id"])
                changed = True
    manifest["blocked_task_ids"] = sorted(blocked)
    manifest["target_task_ids"] = [task["task_id"] for task in tasks if task["task_id"] not in blocked]
    return manifest


def start_project(config: dict[str, Any]) -> dict[str, Any]:
    from .pipeline import build_manifest
    output_root = Path(config["paths"]["output_root"])
    for key in ("output_root", "state_root", "log_root", "temp_root"):
        Path(config["paths"][key]).mkdir(parents=True, exist_ok=True)
    lock_path = output_root / ".project.lock"
    lock_fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o600)
    inherited = False
    try:
        try:
            fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            active = _active_run(config)
            return {"status": "already_running", "run_dir": str(active) if active else None,
                    "reason": "Output directory is locked by an active controller"}
        prior = _active_run(config)
        if prior and _is_live(prior):
            return {"status": "already_running", "run_dir": str(prior), "reason": "Verified controller or scheduler still active"}
        check_resources(config)
        scheduler_version = _check_dependencies(config)
        code_fingerprint = _code_fingerprint(config)
        software = _software_inventory(config)
        software_fingerprint = hashlib.sha256(json.dumps(software, sort_keys=True).encode()).hexdigest()
        run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
        run_dir = output_root / run_id
        run_dir.mkdir()
        manifest = build_manifest(config, run_dir)
        resumed = False
        if prior and config["runtime"]["resume"]["reuse_verified_outputs"]:
            previous = read_json(prior / "run.json")
            same_inputs = manifest.get("input_fingerprint") and manifest["input_fingerprint"] == previous.get("input_fingerprint")
            same_science = config["_scientific_hash"] == previous.get("scientific_hash")
            same_code = code_fingerprint == previous.get("code_fingerprint") and scheduler_version == previous.get("scheduler_version")
            same_software = software_fingerprint == previous.get("software_fingerprint")
            if same_inputs and same_science and same_code and same_software:
                # 留下本次预检记录，避免覆盖上一轮输入身份。
                atomic_json(run_dir / "preflight.json", {"status": "resumed_existing", "resume_run_dir": str(prior), "time": utc_now()})
                run_dir, run_id = prior, prior.name
                manifest = build_manifest(config, run_dir)
                resumed = True
        manifest = _prepare_resume(run_dir, manifest, config)
        atomic_json(run_dir / "manifest.json", manifest)
        history = run_dir / "config_history"
        history.mkdir(exist_ok=True)
        atomic_json(history / f"{uuid.uuid4().hex}.json", config)
        atomic_json(run_dir / "resolved_config.json", config)
        if (run_dir / "stop_request.json").exists():
            history = run_dir / "control_history"
            history.mkdir(exist_ok=True)
            os.replace(run_dir / "stop_request.json", history / f"stop_{uuid.uuid4().hex}.json")
        meta = {"run_id": run_id, "run_dir": str(run_dir), "status": "starting", "started_at": utc_now(),
                "config_path": config["_config_path"], "config_hash": config["_config_hash"],
                "scientific_hash": config["_scientific_hash"], "input_fingerprint": manifest.get("input_fingerprint"),
                "project_root": str(_root(config)), "scheduler_version": scheduler_version, "resumed": resumed}
        meta["code_fingerprint"] = code_fingerprint
        meta["software_fingerprint"] = software_fingerprint
        atomic_json(run_dir / "software_versions.json", software)
        try:
            meta["git_commit"] = subprocess.check_output(["git", "-C", str(_root(config)), "rev-parse", "HEAD"],
                                                        stderr=subprocess.DEVNULL, text=True).strip()
        except (OSError, subprocess.CalledProcessError):
            meta["git_commit"] = None
        log_dir = Path(config["paths"]["log_root"]) / run_id
        log_dir.mkdir(parents=True, exist_ok=True)
        meta["controller_log"] = str(log_dir / "controller.log")
        meta["scheduler_log"] = str(log_dir / "snakemake.log")
        atomic_json(run_dir / "run.json", meta)
        pointer = {"run_id": run_id, "run_dir": str(run_dir), "updated_at": utc_now()}
        atomic_json(_pointer(config), pointer)
        atomic_json(Path(config["paths"]["state_root"]) / "latest.json", pointer)
        append_event(run_dir / "control_events.jsonl", {"event": "resume" if resumed else "start", "run_id": run_id})
        command = [config["paths"]["python"], "-m", "wmh_gmv.lifecycle", "_supervise",
                   "--run-dir", str(run_dir), "--lock-fd", str(lock_fd)]
        with Path(meta["controller_log"]).open("a", encoding="utf-8") as log:
            process = subprocess.Popen(command, cwd=_root(config), env=_environment(config),
                                       stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True, pass_fds=(lock_fd,))
        inherited = True
        os.close(lock_fd)
        lock_fd = -1
        deadline = time.monotonic() + config["runtime"]["execution"]["startup_timeout_seconds"]
        while time.monotonic() < deadline:
            meta = read_json(run_dir / "run.json")
            if meta.get("supervisor") or process.poll() is not None:
                break
            time.sleep(0.1)
        if process.poll() is not None and not meta.get("supervisor"):
            meta.update(status="controller_failed", reason="Controller exited before registering its identity", finished_at=utc_now())
            atomic_json(run_dir / "run.json", meta)
        if not config["runtime"]["execution"]["background"]:
            process.wait()
            meta = read_json(run_dir / "run.json")
        return {"status": meta["status"], "run_id": run_id, "run_dir": str(run_dir),
                "resumed": resumed, "scheduler_log": meta["scheduler_log"], "controller_log": meta["controller_log"],
                "monitor_command": shlex.join([config["paths"]["python"], str(_root(config) / "scripts" / "monitor_project.py"), "--config", config["_config_path"]]),
                "stop_command": shlex.join(["bash", str(_root(config) / "scripts" / "stop_project.sh"), "--config", config["_config_path"]])}
    finally:
        if lock_fd >= 0:
            if not inherited:
                fcntl.flock(lock_fd, fcntl.LOCK_UN)
            os.close(lock_fd)


def completion_status(run_dir: Path, tasks: list[dict[str, Any]], return_code: int) -> str:
    states = [inspect_task(run_dir, item) for item in tasks]
    statuses = [item["status"] for item in states]
    for task, state in zip(tasks, states):
        if has_failed_patients(task, state):
            statuses.append("failed")
        if task.get("unit") == "patient_batch" and state["status"] == "completed":
            try:
                progress = read_json(Path(task["output_dir"]) / "subjects.json")
                statuses.extend(item.get("status", "unknown") for item in progress.values())
            except (OSError, ValueError):
                statuses.append("unknown")
    if (run_dir / "stop_request.json").exists():
        return "stopped"
    if return_code != 0 or "failed" in statuses or "unknown" in statuses:
        return "completed_with_failures"
    if "waiting_input" in statuses or "pending" in statuses:
        return "waiting_input"
    return "completed"


def finalize_report(run_dir: Path, config: dict[str, Any]) -> dict[str, Any]:
    """终态后只汇总已有产物；报告失败不改写科学运行状态。"""
    meta = read_json(run_dir / "run.json")
    original_status = meta["status"]
    if original_status not in ("completed", "completed_with_failures", "stopped", "waiting_input", "controller_failed"):
        return meta
    try:
        from .report import create_report
    except ModuleNotFoundError as exc:
        if exc.name != "wmh_gmv.report":
            # 模块存在但其依赖缺失属于报告失败，不能假装模块未提供。
            meta.update(report_status="failed", report_error=str(exc), report_finished_at=utc_now())
        else:
            meta.update(report_status="skipped_module_unavailable", report_path=None, report_finished_at=utc_now())
        atomic_json(run_dir / "run.json", meta)
        return meta
    except Exception as exc:
        error_path = run_dir / f"report_import_error_{uuid.uuid4().hex}.log"
        error_path.write_text(traceback.format_exc(), encoding="utf-8")
        meta.update(report_status="failed", report_error=str(exc), report_error_log=str(error_path), report_finished_at=utc_now())
        atomic_json(run_dir / "run.json", meta)
        return meta
    report_dir = within(run_dir / "report", run_dir)
    meta.update(report_status="running", report_path=str(report_dir / "report.html"), report_started_at=utc_now())
    atomic_json(run_dir / "run.json", meta)
    try:
        if report_dir.exists():
            history = run_dir / "report_history"
            history.mkdir(exist_ok=True)
            os.replace(report_dir, history / uuid.uuid4().hex)
        report_dir.mkdir()
        result = create_report(run_dir, report_dir, config)
        if result.get("status") != "completed":
            raise RuntimeError(f"Report did not complete: {result.get('status')}: {result.get('reason', '')}")
        names = result.get("artifacts", [])
        if not names or not (report_dir / "report.html").is_file():
            raise RuntimeError("Report did not produce its index and declared artifacts")
        receipts = artifact_receipts(report_dir, names)
        meta.update(report_status="completed", report_artifacts=receipts, report_finished_at=utc_now())
        meta.pop("report_error", None)
        meta.pop("report_error_log", None)
        append_event(run_dir / "control_events.jsonl", {"event": "report_completed", "report_path": meta["report_path"]})
    except Exception as exc:
        error_path = run_dir / f"report_error_{uuid.uuid4().hex}.log"
        error_path.write_text(traceback.format_exc(), encoding="utf-8")
        meta.update(report_status="failed", report_error=str(exc), report_error_log=str(error_path), report_finished_at=utc_now())
        append_event(run_dir / "control_events.jsonl", {"event": "report_failed", "reason": str(exc)})
    meta["status"] = original_status
    atomic_json(run_dir / "run.json", meta)
    return meta


def supervise(run_dir: Path, lock_fd: int) -> int:
    config = read_json(run_dir / "resolved_config.json")
    meta = read_json(run_dir / "run.json")
    meta["supervisor"] = process_identity(os.getpid())
    meta["heartbeat_at"] = utc_now()
    atomic_json(run_dir / "run.json", meta)
    execution = config["runtime"]["execution"]
    command = [config["paths"]["python"], "-m", "snakemake", "--snakefile", str(_root(config) / "workflow" / "Snakefile"),
               "--directory", str(_root(config)), "--cores", str(execution["cores"]),
               "--jobs", str(execution["max_parallel_jobs"]), "--resources", f"mem_mb={execution['memory_mb']}",
               f"project_job={execution['max_parallel_jobs']}",
               "--rerun-incomplete", "--retries", str(execution["automatic_retries"]),
               "--config", f"manifest={run_dir / 'manifest.json'}", f"run_dir={run_dir}",
               f"resolved_config={run_dir / 'resolved_config.json'}"]
    if execution["keep_going"]:
        command.append("--keep-going")
    sent = False
    try:
        if (run_dir / "stop_request.json").exists():
            meta.update(status="stopped", finished_at=utc_now(), reason="Stopped before scheduler launch")
            atomic_json(run_dir / "run.json", meta)
            return 0
        with Path(meta["scheduler_log"]).open("a", encoding="utf-8") as log:
            process = subprocess.Popen(command, cwd=_root(config), env=_environment(config), stdin=subprocess.DEVNULL,
                                       stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            meta["scheduler"] = process_identity(process.pid)
            meta["status"] = "running"
            atomic_json(run_dir / "run.json", meta)
            while process.poll() is None:
                meta["scheduler_ready"] = handles_sigterm(process.pid)
                if (run_dir / "stop_request.json").exists():
                    meta["status"] = "draining"
                    request = read_json(run_dir / "stop_request.json")
                    timeout = config["runtime"]["shutdown"]["drain_timeout_seconds"]
                    if timeout is not None:
                        elapsed = (datetime.now(timezone.utc) - datetime.fromisoformat(request["requested_at"])).total_seconds()
                        meta["drain_timeout_exceeded"] = elapsed > timeout
                    if not sent and meta["scheduler_ready"]:
                        alive, reason = check_process(meta["scheduler"])
                        if alive != "alive":
                            raise RuntimeError(f"Refusing signal: {reason}")
                        os.kill(process.pid, signal.SIGTERM)  # 仅调度器；不发送到进程组。
                        sent = True
                        append_event(run_dir / "control_events.jsonl", {"event": "drain_signal_sent", "pid": process.pid})
                meta["heartbeat_at"] = utc_now()
                atomic_json(run_dir / "run.json", meta)
                time.sleep(min(1.0, config["runtime"]["monitor"]["heartbeat_seconds"]))
            return_code = process.wait()
        tasks = read_json(run_dir / "manifest.json")["tasks"]
        status = completion_status(run_dir, tasks, return_code)
        meta.update(status=status, finished_at=utc_now(), scheduler_exit_code=return_code)
        atomic_json(run_dir / "run.json", meta)
        append_event(run_dir / "control_events.jsonl", {"event": status, "scheduler_exit_code": return_code})
        return 0 if status in ("completed", "stopped", "waiting_input") else 1
    except Exception as exc:
        meta.update(status="controller_failed", reason=str(exc), finished_at=utc_now())
        atomic_json(run_dir / "run.json", meta)
        raise
    finally:
        try:
            finalize_report(run_dir, config)
        except Exception as exc:
            # 即使文件系统故障无法写报告回执，也保留此前已落盘的科学终态。
            print(json.dumps({"event": "report_finalization_failed", "reason": str(exc)}), flush=True)
        finally:
            os.close(lock_fd)


def stop_project(config: dict[str, Any]) -> dict[str, Any]:
    run_dir = _active_run(config)
    if run_dir is None:
        return {"status": "not_started", "reason": "No active run"}
    meta = read_json(run_dir / "run.json")
    if Path(meta["config_path"]).resolve() != Path(config["_config_path"]).resolve():
        raise RuntimeError("Refusing stop: source config identity differs")
    if Path(meta["project_root"]).resolve() != _root(config).resolve():
        raise RuntimeError("Refusing stop: project identity differs")
    if meta["scientific_hash"] != config["_scientific_hash"]:
        raise RuntimeError("Refusing stop: scientific configuration changed; restore the running configuration")
    supervisor_state, supervisor_reason = check_process(meta.get("supervisor"))
    scheduler_state, scheduler_reason = check_process(meta.get("scheduler"))
    if "mismatch" in (supervisor_state, scheduler_state):
        raise RuntimeError(f"Refusing stop: process identity mismatch ({supervisor_reason}; {scheduler_reason})")
    if supervisor_state != "alive" and scheduler_state != "alive":
        if meta["status"] in ("completed", "stopped", "waiting_input", "completed_with_failures"):
            return {"status": meta["status"], "run_dir": str(run_dir), "reason": "No active verified process"}
        raise RuntimeError(f"Process identity requires reconciliation: {supervisor_reason}; {scheduler_reason}")
    if not (run_dir / "stop_request.json").exists():
        atomic_json(run_dir / "stop_request.json", {"requested_at": utc_now(), "run_id": meta["run_id"],
                                                   "config_path": config["_config_path"], "mode": "drain"})
    if supervisor_state != "alive" and scheduler_state == "alive":
        if not handles_sigterm(meta["scheduler"]["pid"]):
            raise RuntimeError("Scheduler has not registered its drain handler; stop request retained")
        os.kill(meta["scheduler"]["pid"], signal.SIGTERM)
    deadline = time.monotonic() + config["runtime"]["shutdown"]["wait_for_exit_seconds"]
    while time.monotonic() < deadline:
        meta = read_json(run_dir / "run.json")
        if meta["status"] in ("stopped", "completed", "waiting_input", "completed_with_failures"):
            return {"status": meta["status"], "run_dir": str(run_dir)}
        if supervisor_state != "alive" and check_process(meta.get("scheduler"))[0] == "gone":
            remaining = [inspect_task(run_dir, task) for task in read_json(run_dir / "manifest.json")["tasks"]]
            if not any(item["status"] in ("running", "unknown") for item in remaining):
                meta.update(status="stopped", finished_at=utc_now(), reason="Orphan scheduler drained; no active task remains")
                atomic_json(run_dir / "run.json", meta)
                return {"status": "stopped", "run_dir": str(run_dir)}
        time.sleep(0.25)
    return {"status": "draining", "run_dir": str(run_dir),
            "reason": "Current jobs are still finishing; no force kill was sent"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["start", "stop", "_supervise"])
    parser.add_argument("--config", type=Path)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--lock-fd", type=int)
    args = parser.parse_args()
    try:
        if args.action == "_supervise":
            return supervise(args.run_dir, args.lock_fd)
        if args.config is None:
            parser.error("--config is required")
        from .config import load_config
        config = load_config(args.config)
        result = start_project(config) if args.action == "start" else stop_project(config)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 2 if result["status"] in ("controller_failed", "error") else 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(json.dumps({"status": "error", "reason": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
