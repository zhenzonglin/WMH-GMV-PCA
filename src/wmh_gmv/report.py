"""生成离线研究结果索引；状态、统计结果和未执行原因同时保留。"""
from __future__ import annotations

import html
import json
import shutil
from pathlib import Path

import pandas as pd
from .state import verify_artifacts


def create_report(run_dir: Path, output: Path, config: dict) -> dict:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    statuses, artifacts, verified_tasks = [], [], set()
    for task in manifest["tasks"]:
        path = run_dir / "tasks" / task["task_id"] / "state.json"
        state = json.loads(path.read_text()) if path.is_file() else {"status": "pending"}
        valid_outputs = state["status"] in {"completed", "not_estimable", "waiting_input"} and verify_artifacts(state)[0]
        if valid_outputs:
            verified_tasks.add(task["task_id"])
        statuses.append({"task_id": task["task_id"], "stage": task["stage"], "cohort": task["cohort"],
                         "status": state["status"], "reason": state.get("reason", ""),
                         "current_artifacts_verified": valid_outputs,
                         "n_subjects": state.get("n_subjects", task.get("n_subjects"))})
        directory = Path(task["output_dir"])
        if valid_outputs and directory.exists():
            for filename in ("metrics.json", "selection.tsv", "effects.tsv", "block_tests.tsv",
                             "diagnostics.tsv", "variable_audit.tsv", "eligibility.tsv", "coverage.tsv"):
                candidate = directory / filename
                if candidate.is_file():
                    artifacts.append({"task_id": task["task_id"], "file": filename,
                                      "path": str(candidate.relative_to(run_dir))})
    pd.DataFrame(statuses).to_csv(output / "stage_results.tsv", sep="\t", index=False)
    pd.DataFrame(artifacts).to_csv(output / "artifact_index.tsv", sep="\t", index=False)
    content = ["<!doctype html><html lang='zh-CN'><meta charset='utf-8'>",
               "<title>WMH-GMV-PCA 分析记录</title><style>body{font:16px system-ui;margin:32px;line-height:1.6;color:#183044}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border:1px solid #ccd6dd;padding:6px;text-align:left}img{max-width:100%}.tag{background:#e7eff6;padding:12px}</style>",
               "<h1>WMH-GMV-PCA 分析记录</h1>",
               f"<p class='tag'>运行 {html.escape(run_dir.name)} · 数据模式 {html.escape(config['project']['mode'])}</p>",
               "<p>合成数据结果只用于软件验证。不可估计、等待输入与已检验但未显著必须分别解释。</p>",
               "<h2>各阶段记录</h2>", pd.DataFrame(statuses).to_html(index=False, escape=True),
               "<h2>统计结果索引</h2>", pd.DataFrame(artifacts).to_html(index=False, escape=True)]
    images = output / "figures"
    images.mkdir(exist_ok=True)
    for task in manifest["tasks"]:
        if task["task_id"] not in verified_tasks:
            continue
        directory = Path(task["output_dir"])
        for image in sorted(directory.glob("*.png")):
            name = task["task_id"] + "_" + image.name
            shutil.copy2(image, images / name)
            content.append(f"<h3>{html.escape(task['task_id'])}</h3><img src='figures/{html.escape(name, quote=True)}'>")
    content.append("<p>GM119标准空间脑图未配置时，仅提供区域载荷表和热图；不把热图称为体素脑图。</p></html>")
    (output / "report.html").write_text("\n".join(content), encoding="utf-8")
    summary = {"run_id": run_dir.name, "mode": config["project"]["mode"], "task_count": len(statuses),
               "status_counts": pd.Series([s["status"] for s in statuses]).value_counts().to_dict(),
               "scientific_hash": config["_scientific_hash"]}
    (output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    return {"status": "completed", "artifacts": [str(p.relative_to(output)) for p in output.rglob("*") if p.is_file()]}
