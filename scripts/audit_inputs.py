#!/usr/bin/env python3
"""只读取三队列输入，生成缺项检查；不拟合模型。"""
import argparse
import json
from pathlib import Path

from wmh_gmv.config import load_config
from wmh_gmv.ingest import collect_cohort
from wmh_gmv.ingest import normalize_subject
from wmh_gmv.clinical import prepare_clinical
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
args = parser.parse_args()
cfg = load_config(args.config)
report = {}
for cohort in ("LS", "C3", "GB"):
    try:
        frame, features, sources = collect_cohort(cfg, cohort)
        report[cohort] = {"status": "readable", "n": len(frame),
                          "features": {k: len(v) for k, v in features.items()},
                          "clinical_table_configured": cfg["paths"]["cohorts"][cohort]["clinical_table"] is not None,
                          "source_files": len(sources)}
        normalized = [normalize_subject(row,features,cfg)[0] for row in frame.to_dict(orient="records")]
        directory = Path(cfg["paths"]["state_root"]) / "audit" / cohort
        directory.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(normalized).to_csv(directory / "imaging_input.tsv",sep="\t",index=False)
        report[cohort]["reference_ready"] = {modality:sum(row[modality+"_ready"] for row in normalized)
                                               for modality in ("wmh","gmv")}
        if cohort != "LS":
            clinical = cfg["paths"]["cohorts"][cohort]["clinical_table"]
            clinical_path = Path(clinical) if clinical and Path(clinical).is_file() else None
            result = prepare_clinical(directory / "imaging_input.tsv",clinical_path,cfg,directory,cohort)
            report[cohort]["clinical_audit"] = str(directory / "variable_audit.tsv")
            report[cohort]["clinical_file_available"] = clinical_path is not None
    except (OSError, ValueError, KeyError) as exc:
        report[cohort] = {"status": "waiting_input", "n": None, "reason": str(exc)}
report["atlas_inputs"] = [{"name":a["name"],"image_exists":Path(a["path"]).is_file(),
                           "labels_exist":Path(a["labels"]).is_file(),"path":a["path"],"template":a["template"]}
                          for a in cfg["paths"]["atlases"]]
output = Path(cfg["paths"]["state_root"]) / "input_audit.json"
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"report": str(output), "cohorts": report}, ensure_ascii=False, indent=2))
