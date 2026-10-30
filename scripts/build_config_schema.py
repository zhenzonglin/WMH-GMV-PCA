"""开发工具：从已审定模板生成严格结构；发布时使用提交的 schema.json。"""
import json
from pathlib import Path
import yaml

root = Path(__file__).resolve().parents[1]
master = yaml.safe_load((root / "config/project.yaml").read_text())
value = {"project": master["project"], "stages": master["stages"]}
for key, name in master["includes"].items():
    value[key] = yaml.safe_load((root / "config" / name).read_text())


def schema(item, path=()):
    if path in {("variables", "sex_map"), ("variables", "clinical", "fields"),
                ("variables", "clinical", "value_maps"), ("variables", "clinical", "source_value_maps"),
                ("variables", "clinical", "field_missing_codes")}:
        return {"type": "object"}
    if isinstance(item, dict):
        return {"type": "object", "additionalProperties": False,
                "required": list(item),
                "properties": {k: schema(v, (*path, k)) for k, v in item.items()}}
    if isinstance(item, list):
        if path == ("variables", "clinical", "missing_codes"):
            return {"type": "array"}
        return {"type": "array", "items": schema(item[0], (*path, "[]")) if item else {}}
    if item is None:
        return {"type": ["string", "null"]} if path[-1] != "drain_timeout_seconds" else {"type": ["number", "null"], "minimum": 0}
    if isinstance(item, bool):
        return {"type": "boolean"}
    if isinstance(item, int):
        return {"type": "integer", "minimum": 0}
    if isinstance(item, float):
        return {"type": "number", "minimum": 0}
    return {"type": "string"}


s = schema(value)
s["$schema"] = "https://json-schema.org/draft/2020-12/schema"
s["properties"]["project"]["properties"]["mode"]["enum"] = ["production", "synthetic"]
for key in ["cores", "max_parallel_jobs", "blas_threads_per_job", "memory_mb", "startup_timeout_seconds"]:
    s["properties"]["runtime"]["properties"]["execution"]["properties"][key]["minimum"] = 1
for key in ["refresh_seconds", "heartbeat_seconds", "stale_heartbeat_seconds"]:
    s["properties"]["runtime"]["properties"]["monitor"]["properties"][key]["minimum"] = 1
for key in ["cdf_clip"]:
    s["properties"]["analysis"]["properties"]["normative"]["properties"][key].update({"exclusiveMinimum": 0, "exclusiveMaximum": .5})
s["properties"]["analysis"]["properties"]["pca"]["properties"]["mask_fraction"].update({"exclusiveMinimum": 0, "exclusiveMaximum": 1})
(root / "config/schema.json").write_text(json.dumps(s, indent=2), encoding="utf-8")
