"""从统一配置安装R依赖，安装日志写入配置指定目录。"""
import argparse
import json
import os
import subprocess
from pathlib import Path
from wmh_gmv.config import load_config

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
args = parser.parse_args()
cfg = load_config(args.config)
temp = Path(cfg["paths"]["temp_root"])
logs = Path(cfg["paths"]["log_root"])
temp.mkdir(parents=True, exist_ok=True)
logs.mkdir(parents=True, exist_ok=True)
request = temp / "setup-config.json"
request.write_text(json.dumps(cfg, ensure_ascii=False, indent=2))
env = dict(os.environ)
env["R_LIBS"] = cfg["paths"]["r_library"] + os.pathsep + env.get("R_LIBS", "")
env["PATH"] = str(Path(cfg["paths"]["python"]).parent) + os.pathsep + env["PATH"]
log_path = logs / "install_r.log"
with log_path.open("w") as log:
    bootstrap = (
        "a<-commandArgs(TRUE);dir.create(a[1],recursive=TRUE,showWarnings=FALSE);"
        ".libPaths(c(a[1],.libPaths()));options(repos=c(CRAN=a[2]));"
        "p<-c('jsonlite','renv');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];"
        "if(length(m))install.packages(m,lib=a[1])"
    )
    subprocess.run([cfg["paths"]["rscript"], "-e", bootstrap, cfg["paths"]["r_library"],
                    cfg["runtime"]["environment"]["cran_repository"]], env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    result = subprocess.run([cfg["paths"]["rscript"], str(Path(cfg["project"]["project_root"]) / "scripts/install_r.R"), str(request)],
                            env=env, stdout=log, stderr=subprocess.STDOUT)
print(json.dumps({"returncode": result.returncode, "log": str(log_path)}))
raise SystemExit(result.returncode)
