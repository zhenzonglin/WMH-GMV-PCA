"""引导环境只读配置，实际入口始终交给配置指定的Python执行。"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from .config import load_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["start", "stop", "monitor"])
    parser.add_argument("--config", required=True)
    args, remaining = parser.parse_known_args(argv)
    config = load_config(args.config)
    executable = config["paths"]["python"]
    module = "wmh_gmv.monitor" if args.action == "monitor" else "wmh_gmv.lifecycle"
    command = [executable, "-m", module]
    if args.action != "monitor":
        command.append(args.action)
    command += ["--config", config["_config_path"], *remaining]
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(Path(config["project"]["project_root"]) / "src") + os.pathsep + environment.get("PYTHONPATH", "")
    os.execve(executable, command, environment)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
