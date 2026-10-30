"""交接数据只做整理与缺项说明，不拟合出血模型。"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd

from .ingest import sha256_file


def package_handoff(source: Path, output: Path, cohort: str, config: dict) -> dict:
    files = ["clinical_wide.tsv", "clinical_long.tsv", "variable_audit.tsv", "eligibility.tsv",
             "conflicts.tsv", "events.tsv", "medication_snapshots.tsv", "imaging.tsv", "clinical_source.tsv", "clinical_source.csv"]
    copied = []
    for name in files:
        if (source / name).is_file():
            shutil.copy2(source / name, output / name)
            copied.append(name)
    if "clinical_wide.tsv" not in copied:
        return {"status": "waiting_input", "reason": "Patient feature table unavailable", "artifacts": copied}
    data = pd.read_csv(output / "clinical_wide.tsv", sep="\t", dtype=str, keep_default_na=False)
    dictionary = []
    for column in data:
        definition, unit = "临床字段；编码和来源见变量配置及审计表", "source"
        if column.endswith("_ml"):
            definition, unit = "原始结构体积；ICV使用真实颅内容积", "mL"
        if column.endswith("_dev"):
            definition, unit = "LS年龄/性别/ICV健康参考偏离；正值为WMH增多或GMV减少", "reference deviation"
        if "_PC" in column:
            definition, unit = "固定C3空间模式投射；方向任意，需结合载荷解释", "PCA score"
        if column.endswith("_C3SD"):
            definition, unit = "相对于C3训练集的标准差单位", "C3 SD"
        if column.startswith("JOINT_"):
            definition, unit = "按各模态C3总方差平衡后的联合坐标", "block-scaled score"
        dictionary.append({"variable": column, "definition": definition, "unit": unit,
                           "missing": "NA/empty means unknown or unavailable, never no-event"})
    pd.DataFrame(dictionary).to_csv(output / "feature_dictionary.tsv", sep="\t", index=False)
    checklists = {
        "B3_EVT_bleeding": ["EVT及操作时间", "出血首次出现时间或区间", "复查影像时间", "出血类型",
                            "症状性判定标准及NIHSS变化", "再通等级", "围术期用药"],
        "B4_antithrombotic": ["药物种类/剂量/适应证", "起停换药时间", "脑出血类型及日期",
                              "缺血复发日期", "死亡和末次随访", "CMB/cSS可判读状态", "肾功能和既往出血"]}
    (output / "branch_requirements.json").write_text(json.dumps(checklists, ensure_ascii=False, indent=2))
    (output / "variables_snapshot.json").write_text(json.dumps(config["variables"], ensure_ascii=False, indent=2))
    (output / "README_统计师交接.md").write_text(
        f"# {cohort}统计师交接包\n\n本包不含B3/B4拟合结果。clinical_wide为一人一行，clinical_long为随访长表。"
        "events仅保留记录，不将无记录视为无事件。medication_snapshots是用药快照，不能直接作为连续用药区间。\n\n"
        "偏离、PCA和联合坐标必须使用本轮固定定义。WMH与GMV区域不能逐编号配对。"
        "缺项以variable_audit、conflicts及branch_requirements为准；未记录、不可判读和否定记录须区分。\n",
        encoding="utf-8")
    manifest = {p.name: sha256_file(p) for p in output.iterdir() if p.is_file()}
    (output / "SHA256.json").write_text(json.dumps(manifest, indent=2))
    return {"status": "completed", "artifacts": [p.name for p in output.iterdir() if p.is_file()],
            "n_subjects": len(data)}
