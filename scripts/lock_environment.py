"""在测试通过后锁定实际Python与R依赖；不启动科学分析。"""
import argparse
import importlib.metadata
import json
import os
import subprocess
from pathlib import Path
from wmh_gmv.config import load_config

parser = argparse.ArgumentParser()
parser.add_argument("--config",required=True)
args = parser.parse_args()
cfg = load_config(args.config)
root = Path(cfg["project"]["project_root"])
packages = {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()
            if d.metadata["Name"] and d.metadata["Name"].lower() != "wmh-gmv-pca"}
(root / "requirements.lock").write_text("# Verified environment; regenerate only after tests.\n" +
    "\n".join(f"{name}=={version}" for name,version in sorted(packages.items(),key=lambda x:x[0].lower())) + "\n")
path = Path(cfg["paths"]["temp_root"]) / "lock-config.json"
path.parent.mkdir(parents=True,exist_ok=True)
path.write_text(json.dumps(cfg))
env = dict(os.environ)
env["R_LIBS"] = cfg["paths"]["r_library"] + os.pathsep + env.get("R_LIBS", "")
subprocess.run([cfg["paths"]["rscript"],str(root / "scripts/lock_environment.R"),str(path)],env=env,check=True)
print(json.dumps({"python_packages":len(packages),"python_lock":str(root / "requirements.lock"),"r_lock":str(root / "renv.lock")}))
