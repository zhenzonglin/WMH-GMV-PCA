"""统一配置读取：路径相对项目根目录，科学身份与运行参数分开。"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any

import jsonschema
import yaml


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                    allow_nan=False).encode()).hexdigest()


def _yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Configuration must be a mapping: {path}")
    return value


def load_config(path: str | Path) -> dict:
    source = Path(path).expanduser().resolve(strict=True)
    master = _yaml(source)
    extra = set(master) - {"project", "includes", "stages"}
    if extra:
        raise ValueError(f"Unknown main configuration keys: {sorted(extra)}")
    if set(master.get("includes", {})) != {"paths", "analysis", "variables", "runtime"}:
        raise ValueError("includes must define paths, analysis, variables and runtime")
    cfg = {"project": master["project"], "stages": master["stages"]}
    for section, filename in master["includes"].items():
        cfg[section] = _yaml((source.parent / filename).resolve())
    schema_path = Path(__file__).resolve().parents[2] / "config" / "schema.json"
    if schema_path.exists():
        validator = jsonschema.Draft202012Validator(json.loads(schema_path.read_text()))
        errors = sorted(validator.iter_errors(cfg), key=lambda e: str(list(e.path)))
        if errors:
            raise ValueError("Invalid configuration:\n" + "\n".join(
                f"{'.'.join(map(str, e.path))}: {e.message}" for e in errors))
    else:
        raise ValueError(f"Required configuration schema missing: {schema_path}")
    root = Path(cfg["project"]["project_root"]).expanduser()
    if not root.is_absolute():
        root = source.parent / root
    root = root.resolve()
    cfg["project"]["project_root"] = str(root)

    def resolve(value: str | None) -> str | None:
        if value is None:
            return None
        item = Path(value).expanduser()
        return str((item if item.is_absolute() else root / item).resolve())

    for name in ("output_root", "state_root", "log_root", "temp_root", "r_library"):
        cfg["paths"][name] = resolve(cfg["paths"][name])
    for name in ("python", "rscript"):
        executable = Path(cfg["paths"][name]).expanduser()
        cfg["paths"][name] = os.path.abspath(executable if executable.is_absolute() else root / executable)
    for cohort in cfg["paths"]["cohorts"].values():
        for name, value in cohort.items():
            cohort[name] = resolve(value)
    for atlas in cfg["paths"]["atlases"]:
        atlas["path"] = resolve(atlas["path"])
        atlas["labels"] = resolve(atlas["labels"])
    if cfg["analysis"]["age_range"][0] >= cfg["analysis"]["age_range"][1]:
        raise ValueError("age_range must be increasing")
    if cfg["runtime"]["shutdown"]["force_kill"]:
        raise ValueError("Automatic force kill is not supported by the safe shutdown contract")
    if not cfg["runtime"]["resume"]["require_matching_scientific_config"]:
        raise ValueError("Scientific identity matching must remain enabled")
    if not cfg["stages"]["representation"] and any(cfg["stages"][key] for key in
            ("acute_location", "clinical_b1", "clinical_b2", "export_statistician")):
        raise ValueError("Dependent analysis stages require representation=true; disable all for audit-only mode")
    if cfg["stages"]["clinical_b2"] and not cfg["stages"]["acute_location"]:
        raise ValueError("clinical_b2 requires acute_location=true")
    scientific = {key: cfg[key] for key in ("analysis", "variables", "stages")}
    scientific["mode"] = cfg["project"]["mode"]
    scientific["atlas_definitions"] = [
        {k: v for k, v in a.items() if k not in {"path", "labels"}}
        for a in cfg["paths"]["atlases"]
    ]
    cfg["_scientific_hash"] = digest(scientific)
    cfg["_config_path"] = str(source)
    cfg["_config_hash"] = digest(cfg)
    return cfg


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output")
    args = parser.parse_args()
    cfg = load_config(args.config)
    serialized = json.dumps(cfg, ensure_ascii=False, indent=2, allow_nan=False)
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized)


if __name__ == "__main__":
    main()
