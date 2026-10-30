"""仅引导依赖环境；数据和科学计算仍由正式入口启动。"""
import argparse
import os
import shutil
import subprocess
from pathlib import Path
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
args = parser.parse_args()
source = Path(args.config).resolve()
master = yaml.safe_load(source.read_text(encoding="utf-8-sig"))
root = Path(master["project"]["project_root"])
root = (root if root.is_absolute() else source.parent / root).resolve()
paths = yaml.safe_load((source.parent / master["includes"]["paths"]).read_text())
runtime = yaml.safe_load((source.parent / master["includes"]["runtime"]).read_text())
executable = Path(os.path.abspath(root / paths["python"]))
uv = shutil.which("uv")
if uv is None:
    raise SystemExit("uv is not on PATH")
if not executable.exists():
    subprocess.run([uv, "venv", "--python", runtime["environment"]["minimum_python"], str(executable.parent.parent)], check=True)
isolated = subprocess.check_output([str(executable), "-c", "import sys;print(int(sys.prefix!=sys.base_prefix))"], text=True).strip()
if isolated != "1":
    raise SystemExit("Configured Python must be an isolated virtual environment")
lock = root / "requirements.lock"
requirements = ["-r", str(lock)] if lock.is_file() else [str(root) + "[test,build]"]
subprocess.run([uv, "pip", "install", "--python", str(executable), *requirements], check=True)
subprocess.run([uv, "pip", "install", "--python", str(executable), "--no-deps", "-e", str(root)], check=True)
(root / ".wmh-gmv-python").write_text(str(executable) + "\n", encoding="utf-8")
subprocess.run([str(executable), str(root / "scripts/setup_r.py"), "--config", str(source)], check=True)
