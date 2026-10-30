"""只读监测运行状态；人数与区域/模型/重复任务进度分别报告。"""
from __future__ import annotations

import argparse
import csv
import html
import io
import json
import os
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .state import atomic_json, check_process, inspect_task, read_json, task_dir, utc_now, within


COUNTS = ("completed", "failed", "running", "pending", "unknown", "not_applicable")
LABELS = {"completed": "已完成", "failed": "失败", "running": "进行中", "pending": "待完成",
          "unknown": "待核实", "not_applicable": "不适用", "total": "总人数/任务数"}


def _seconds_since(value: str) -> float:
    stamp = datetime.fromisoformat(value)
    if stamp.tzinfo is None:
        stamp = stamp.replace(tzinfo=timezone.utc)
    return max(0.0, (datetime.now(timezone.utc) - stamp).total_seconds())


def snapshot(config: dict[str, Any]) -> dict[str, Any]:
    pointer = Path(config["paths"]["output_root"]) / ".active-run.json"
    if not pointer.exists():
        return {"status": "not_started", "updated_at": utc_now(), "patient_steps": [], "group_steps": [],
                "errors": [], "reason": "尚未启动；实际名单人数尚未确定", "total_people": None}
    try:
        run_dir = within(read_json(pointer)["run_dir"], config["paths"]["output_root"])
        meta = read_json(run_dir / "run.json")
        manifest = read_json(run_dir / "manifest.json")
        tasks = manifest["tasks"]
        if not isinstance(tasks, list) or any(not isinstance(task, dict) for task in tasks):
            raise ValueError("Invalid task manifest")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return {"status": "unknown", "updated_at": utc_now(), "patient_steps": [], "group_steps": [],
                "errors": [str(exc)], "total_people": None, "reason": "运行记录或固定名单不可读取"}
    rows: dict[tuple[str, str, str], dict[str, Any]] = {}
    units: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    errors, warnings, waiting = [], [], []
    states = {task["task_id"]: inspect_task(run_dir, task) for task in tasks}
    progress_tasks = []
    for task in tasks:
        if task.get("unit") != "patient_batch":
            progress_tasks.append(task)
            continue
        batch_state = states[task["task_id"]]
        source_root = Path(task["output_dir"]) if batch_state["status"] in ("completed", "not_estimable", "waiting_input") else Path(batch_state.get("staging_dir", task["output_dir"]))
        progress_path = source_root / "subjects.json"
        progress = {}
        try:
            if progress_path.exists():
                progress = read_json(progress_path)
                if not isinstance(progress, dict) or any(not isinstance(item, dict) for item in progress.values()):
                    raise ValueError("Invalid per-patient progress object")
            elif batch_state["status"] in ("completed", "not_estimable", "waiting_input", "failed", "unknown"):
                progress = {key: {"status": "unknown", "reason": "Per-patient progress unavailable"} for key in task["subject_keys"]}
        except (OSError, ValueError) as exc:
            progress = {key: {"status": "unknown", "reason": str(exc)} for key in task["subject_keys"]}
        for subject in task["subject_keys"]:
            virtual_id = task["task_id"] + "::" + subject
            item = progress.get(subject, {"status": "pending", "reason": "Not started in batch"})
            item = dict(item)
            if item.get("status") == "running" and batch_state["status"] != "running":
                item["status"] = "unknown" if batch_state["status"] == "unknown" else "pending"
                item["reason"] = "Batch not active; patient checkpoint requires resume"
            if item.get("status") == "running":
                item["process"] = batch_state.get("process")
                item["attempt_id"] = batch_state.get("attempt_id")
            states[virtual_id] = item
            progress_tasks.append({**task, "task_id": virtual_id, "unit": "patient", "subject_key": subject,
                                   "heartbeat_task_id": task["task_id"], "n_subjects": 1})
    for task in progress_tasks:
        unit = task.get("unit", "model")
        cohort, stage = task.get("cohort", "ALL"), task["stage"]
        key = (cohort, stage, unit)
        if key not in rows:
            rows[key] = {"cohort": cohort, "stage": stage, "unit": unit, **dict.fromkeys(COUNTS, 0),
                         "total": 0, "n_subjects": None, "historical_failures": 0}
        row = rows[key]
        entity = task.get("subject_key", task["task_id"])
        if entity in units[key]:
            errors.append(f"Duplicate progress unit in manifest: {cohort}/{stage}/{entity}")
            continue
        units[key].add(entity)
        row["total"] += 1
        state = states[task["task_id"]]
        status = state["status"]
        classified = {"waiting_input": "pending", "not_estimable": "not_applicable", "stopped": "pending"}.get(status, status)
        if classified not in COUNTS:
            classified = "unknown"
        row[classified] += 1
        row["historical_failures"] += state.get("history_failures", 0)
        number = state.get("n_subjects", task.get("n_subjects"))
        if number is not None:
            row["n_subjects"] = max(row["n_subjects"] or 0, number)
        if classified in ("failed", "unknown"):
            errors.append({"task_id": task["task_id"], "status": classified, "reason": state.get("reason", ""),
                           "error_log": state.get("error_log")})
        if status == "waiting_input":
            waiting.append({"task_id": task["task_id"], "reason": state.get("reason", "")})
        if classified == "pending":
            blockers = [dep for dep in task.get("dependencies", [])
                        if states.get(dep, {}).get("status") in ("failed", "waiting_input", "not_estimable", "unknown")]
            if blockers:
                waiting.append({"task_id": task["task_id"], "reason": "等待上游", "dependencies": blockers})
        if classified == "running":
            heartbeat = task_dir(run_dir, task.get("heartbeat_task_id", task["task_id"])) / "heartbeat.json"
            try:
                value = read_json(heartbeat)
                if value["attempt_id"] != state["attempt_id"]:
                    warnings.append(f"{task['task_id']}: heartbeat belongs to another attempt")
                elif _seconds_since(value["time"]) > config["runtime"]["monitor"]["stale_heartbeat_seconds"]:
                    warnings.append(f"{task['task_id']}: heartbeat stale; process alive, algorithm progress unknown")
            except (OSError, ValueError, KeyError) as exc:
                warnings.append(f"{task['task_id']}: heartbeat unavailable: {exc}")
    for row in rows.values():
        # 冻结分母不缩小；不适用和待核实分别显示，纳入扩展计数恒等式。
        row["planned_total"] = row["total"]
        row["count_invariant_ok"] = row["total"] == sum(row[name] for name in COUNTS)
        if row["unit"] == "patient":
            row["n_subjects"] = row["total"]
    current_status = meta["status"]
    if current_status in ("running", "starting", "draining"):
        live_states = [check_process(meta.get(name))[0] for name in ("supervisor", "scheduler")]
        if "alive" not in live_states:
            current_status = "interrupted" if all(item in ("gone", "unknown") for item in live_states) else "unknown"
            warnings.append("No verified live controller/scheduler; recorded status is stale")
    people = {item["subject_key"] for item in progress_tasks if item.get("subject_key")}
    if meta.get("drain_timeout_exceeded"):
        warnings.append("超过配置的排空时限；当前任务仍在收尾，未强制终止")
    return {
        "run_id": meta["run_id"], "run_dir": str(run_dir), "status": current_status,
        "data_mode": config["project"].get("mode", "unknown"),
        "updated_at": utc_now(), "started_at": meta.get("started_at"), "git_commit": meta.get("git_commit"),
        "scientific_hash": meta.get("scientific_hash"), "scheduler_version": meta.get("scheduler_version"),
        "total_people": len(people) if people else manifest.get("n_subjects"),
        "patient_steps": [row for row in rows.values() if row["unit"] == "patient"],
        "group_steps": [row for row in rows.values() if row["unit"] != "patient"],
        "errors": errors, "warnings": warnings, "waiting_input": waiting,
        "logs": {"controller": meta.get("controller_log"), "scheduler": meta.get("scheduler_log")},
        "draining": current_status == "draining" or (run_dir / "stop_request.json").exists() and current_status != "stopped",
    }


def terminal_report(value: dict[str, Any]) -> str:
    total_label = "未知" if value.get("total_people") is None else str(value["total_people"])
    lines = [f"WMH-GMV-PCA | {value.get('data_mode', 'unknown')} | {value.get('run_id', '未启动')} | {value['status']}",
             f"更新时间: {value['updated_at']} | 人数: {total_label}"]
    if value.get("reason"):
        lines.append(value["reason"])
    for name, title in (("patient_steps", "患者步骤"), ("group_steps", "群体计算（计数单位不是患者）")):
        lines += ["", title, "队列  步骤                      单位       完成  失败  进行  待完  待核  不适用  总数  涉及人数"]
        for row in value[name]:
            nums = " ".join(f"{row[key]:5d}" for key in COUNTS)
            participants = "未知" if row["n_subjects"] is None else str(row["n_subjects"])
            lines.append(f"{row['cohort']:5} {row['stage'][:24]:24} {row['unit']:8} {nums} {row['total']:5} {participants:>8}")
    if value.get("warnings"):
        lines += ["", "提醒: " + json.dumps(value["warnings"], ensure_ascii=False)]
    if value.get("errors"):
        lines += ["", "错误: " + json.dumps(value["errors"][-5:], ensure_ascii=False)]
    if value.get("waiting_input"):
        lines += ["", f"等待输入/依赖任务: {len(value['waiting_input'])}，详见 JSON/HTML"]
    if value.get("logs"):
        lines.append("日志: " + str(value["logs"]["scheduler"]))
    return "\n".join(lines)


def html_report(value: dict[str, Any]) -> str:
    escaped = html.escape(terminal_report(value))
    detail = html.escape(json.dumps({"errors": value.get("errors"), "waiting_input": value.get("waiting_input")},
                                  ensure_ascii=False, indent=2))
    return """<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>WMH-GMV-PCA 进度</title>
<style>body{font:15px system-ui;margin:28px;color:#153348;background:#f6f8fa}pre{background:white;padding:20px;border:1px solid #d8e1e8;overflow:auto;line-height:1.7}h1{font-size:24px}summary{cursor:pointer}</style>
<h1>WMH-GMV-PCA 项目进度</h1><p>离线快照；心跳不代表算法迭代进度。</p><pre>""" + escaped + "</pre><details><summary>错误及缺项明细</summary><pre>" + detail + "</pre></details></html>"


def write_reports(value: dict[str, Any], config: dict[str, Any]) -> None:
    """只写派生监测快照，不修改任务回执、运行状态或结果。"""
    if not value.get("run_dir"):
        return
    output = Path(value["run_dir"]) / "monitor"
    output.mkdir(exist_ok=True)
    selected = config["runtime"]["monitor"]["outputs"]
    if "json" in selected:
        atomic_json(output / "progress.json", value)
    for extension, content in (("html", html_report(value)), ("tsv", _tsv(value))):
        if extension in selected:
            target = output / f"progress.{extension}"
            temporary = output / f".progress.{extension}.{os.getpid()}"
            temporary.write_text(content, encoding="utf-8")
            os.replace(temporary, target)


def _tsv(value: dict[str, Any]) -> str:
    handle = io.StringIO()
    names = ["cohort", "stage", "unit", *COUNTS, "total", "planned_total", "n_subjects", "historical_failures", "count_invariant_ok"]
    writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t")
    writer.writeheader()
    writer.writerows(value["patient_steps"] + value["group_steps"])
    return handle.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--once", action="store_true", help="只显示一次快照")
    args = parser.parse_args()
    from .config import load_config
    config = load_config(args.config)
    try:
        while True:
            value = snapshot(config)
            if "terminal" in config["runtime"]["monitor"]["outputs"]:
                if sys.stdout.isatty() and not args.once:
                    print("\033[2J\033[H", end="")
                print(terminal_report(value), flush=True)
            write_reports(value, config)
            if args.once:
                return 0 if value["status"] != "unknown" else 2
            time.sleep(config["runtime"]["monitor"]["refresh_seconds"])
    except KeyboardInterrupt:
        return 0  # 关闭监测不会关闭分析。


if __name__ == "__main__":
    raise SystemExit(main())
