#!/usr/bin/env python
"""从项目内或任何工作目录调用监测入口。"""
import os
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
bootstrap = project_root / ".venv" / "bin" / "python"
pointer = project_root / ".wmh-gmv-python"
if pointer.is_file():
    bootstrap = Path(pointer.read_text(encoding="utf-8").strip())
if Path(sys.prefix).resolve() != bootstrap.parent.parent.resolve():
    if not bootstrap.is_file():
        raise SystemExit("Project bootstrap environment missing; run bash scripts/setup_project.sh first.")
    os.execv(str(bootstrap), [str(bootstrap), str(Path(__file__).resolve()), *sys.argv[1:]])
sys.path.insert(0, str(project_root / "src"))
from wmh_gmv.launcher import main

if __name__ == "__main__":
    raise SystemExit(main(["monitor", *sys.argv[1:]]))
