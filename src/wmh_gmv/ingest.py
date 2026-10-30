"""读取既有体积结果，保留字符串ID和原始数据；不重做影像分割。"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path: str | Path, id_column: str = "participant_id") -> pd.DataFrame:
    path = Path(path)
    content = path.read_bytes()
    frame = pd.read_csv(io.BytesIO(content), sep="\t" if path.suffix.lower() != ".csv" else ",",
                        dtype=str, keep_default_na=False, encoding="utf-8-sig")
    if id_column not in frame:
        raise ValueError(f"Missing ID column {id_column}: {path}")
    if frame[id_column].eq("").any() or frame[id_column].duplicated().any():
        raise ValueError(f"Empty/duplicate IDs in {path}")
    frame = frame.rename(columns={id_column: "participant_id"})
    frame.attrs["source_sha256"] = hashlib.sha256(content).hexdigest()
    return frame


def _number(value):
    try:
        val = float(value)
        return val if math.isfinite(val) else None
    except (ValueError, TypeError):
        return None


def _nested(data: dict, key: str):
    for part in key.split("."):
        if not isinstance(data, dict):
            return None
        data = data.get(part)
    return data


def _discover(spec: dict, modality: str, cfg: dict) -> Path | None:
    value = spec[f"{modality}_table"]
    if value:
        path = Path(value)
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    root = spec.get("derivatives_root")
    candidates = [] if root is None else [Path(root) / p for p in
                    cfg["variables"]["imaging"]["table_candidates"][modality]]
    found = list(dict.fromkeys(p.resolve() for p in candidates if p.is_file()))
    if len(found) > 1:
        raise ValueError(f"Ambiguous {modality} tables; set an explicit path: {found}")
    return found[0] if found else None


def collect_cohort(config: dict, cohort: str) -> tuple[pd.DataFrame, dict, list[Path]]:
    """只读预检并构建输入快照；返回来源清单用于内容指纹。"""
    spec = config["paths"]["cohorts"][cohort]
    variables = config["variables"]
    roster_path = Path(spec["roster"])
    frame = read_table(roster_path, variables["id_column"]).astype(object).set_index("participant_id", drop=False)
    if frame.empty:
        raise ValueError(f"Empty cohort roster: {cohort}; cannot start a reference/validation analysis")
    source_hashes = {str(roster_path): frame.attrs["source_sha256"]}
    sources = [roster_path]
    if spec["metadata_table"]:
        metadata = read_table(spec["metadata_table"], variables["id_column"]).set_index("participant_id")
        source_hashes[str(spec["metadata_table"])] = metadata.attrs["source_sha256"]
        for col in metadata:
            if col in frame:
                aligned = metadata[col].reindex(frame.index).fillna("")
                conflict = frame[col].ne("") & aligned.ne("") & frame[col].ne(aligned)
                if conflict.any():
                    raise ValueError(f"Conflicting roster/metadata values: {cohort}/{col}")
                frame[col] = frame[col].where(frame[col].ne(""), aligned)
            else:
                frame[col] = metadata[col].reindex(frame.index).fillna("")
        sources.append(Path(spec["metadata_table"]))
    for canonical, source in variables["metadata"].items():
        if source in frame and canonical != source:
            frame[canonical] = frame[source]
    snapshots = {}
    if spec["raw_records"]:
        raw_path = Path(spec["raw_records"])
        raw_bytes = raw_path.read_bytes()
        source_hashes[str(raw_path)] = hashlib.sha256(raw_bytes).hexdigest()
        with io.StringIO((gzip.decompress(raw_bytes) if raw_path.suffix == ".gz" else raw_bytes).decode("utf-8")) as handle:
            for line in handle:
                record = json.loads(line)
                pid = str(record["participant_id"])
                if pid in snapshots:
                    raise ValueError("Duplicate participant in raw JSON snapshot")
                snapshots[pid] = record
        sources.append(raw_path)
    for modality in ("wmh", "gmv"):
        table_path = _discover(spec, modality, config)
        if table_path:
            feature_table = read_table(table_path, variables["id_column"]).set_index("participant_id")
            source_hashes[str(table_path)] = feature_table.attrs["source_sha256"]
            prefix = variables["imaging"][f"{modality}_prefix"]
            suffix = variables["imaging"]["volume_suffix"]
            columns = [c for c in feature_table if c.startswith(prefix) and c.endswith(suffix)
                       and c not in {"wmh_total_ml", "gmv_total_ml"}]
            for col in columns:
                frame[col] = feature_table[col].reindex(frame.index).fillna("")
            sources.append(table_path)
    image_cfg = variables["imaging"]
    for pid in frame.index:
        t1, wmh = {}, {}
        snap = snapshots.get(pid, {})
        if snap:
            t1 = snap.get("features", {}).get("t1", {})
            wmh = snap.get("features", {}).get("wmh", {})
            for field in ("age", "sex"):
                if field not in frame or not str(frame.at[pid, field]):
                    frame.at[pid, field] = snap.get("participant_metadata", {}).get(field, "")
        elif spec["derivatives_root"]:
            for key, pattern in (("t1", "feature_json_pattern"), ("wmh", "wmh_json_pattern")):
                path = Path(spec["derivatives_root"]) / image_cfg[pattern].format(participant_id=pid)
                if path.is_file():
                    raw_bytes = path.read_bytes()
                    source_hashes[str(path)] = hashlib.sha256(raw_bytes).hexdigest()
                    data = json.loads(raw_bytes)
                    sources.append(path)
                    if key == "t1":
                        t1 = data
                    else:
                        wmh = data
        for name, value in wmh.get("raw_ml", {}).items():
            if name not in frame or str(frame.at[pid, name]) in {"", "nan"}:
                frame.at[pid, name] = value
        for label, value in t1.get(image_cfg["native_gm_json_key"], {}).items():
            name = f"muse_{label}_ml"
            if name not in frame or str(frame.at[pid, name]) in {"", "nan"}:
                frame.at[pid, name] = value
        if "icv_ml" not in frame or _number(frame.at[pid, "icv_ml"]) is None:
            frame.at[pid, "icv_ml"] = _nested(t1, image_cfg["icv_json_key"])
        if "wmh_total_ml" not in frame or _number(frame.at[pid, "wmh_total_ml"]) is None:
            # 只统计已有校正掩膜，不修改掩膜或以分区总和替代全脑总量。
            if spec["derivatives_root"]:
                mask = Path(spec["derivatives_root"]) / image_cfg["corrected_wmh_pattern"].format(participant_id=pid)
                if mask.is_file():
                    import nibabel as nib
                    before = sha256_file(mask)
                    image = nib.load(mask)
                    if image.header.get_xyzt_units()[0] != config["analysis"]["acute"]["expected_spatial_unit"]:
                        frame.at[pid, "wmh_total_issue"] = "Corrected WMH mask units unavailable/not mm; total not inferred"
                        continue
                    data = np.asarray(image.dataobj)
                    if not np.isin(data, [0, 1]).all():
                        raise ValueError(f"Corrected WMH mask is not binary: {mask}")
                    frame.at[pid, "wmh_total_ml"] = float(data.sum() * abs(np.linalg.det(image.affine[:3, :3])) / 1000)
                    if sha256_file(mask) != before:
                        raise ValueError(f"Input changed during mask read: {mask}")
                    source_hashes[str(mask)] = before
                    sources.append(mask)
    features = {}
    for modality in ("wmh", "gmv"):
        columns = sorted(c for c in frame if c.startswith(image_cfg[f"{modality}_prefix"])
                         and c.endswith(image_cfg["volume_suffix"])
                         and c not in {"wmh_total_ml", "gmv_total_ml"})
        expected = image_cfg[f"expected_{modality}_regions"]
        if len(columns) != expected:
            raise ValueError(f"{cohort}: expected {expected} {modality} regions; observed {len(columns)}")
        features[modality] = columns
    for name in ("age", "sex", "icv_ml", "wmh_total_ml", "gmv_total_ml", "lesion_mask", "lesion_template"):
        if name not in frame:
            frame[name] = ""
    frame["cohort"] = cohort
    frame["subject_key"] = cohort + ":" + frame["participant_id"]
    frame = frame.reset_index(drop=True).fillna("")
    project_root = Path(config["project"]["project_root"])
    frame["lesion_mask"] = frame["lesion_mask"].map(
        lambda value: str((Path(value) if Path(value).is_absolute() else project_root / value).resolve()) if value else "")
    if spec["clinical_table"] and Path(spec["clinical_table"]).is_file():
        sources.append(Path(spec["clinical_table"]))
    for value in frame["lesion_mask"]:
        if value and Path(value).is_file():
            sources.append(Path(value))
    frame.attrs["source_hashes"] = source_hashes
    return frame, features, list(dict.fromkeys(sources))


def normalize_subject(row: dict, features: dict, config: dict) -> tuple[dict, dict]:
    result = dict(row)
    sex_raw = str(row.get("sex", "")).strip()
    result["sex"] = config["variables"]["sex_map"].get(sex_raw)
    issues = []
    numeric = ["age", "icv_ml", "wmh_total_ml", "gmv_total_ml"] + features["wmh"] + features["gmv"]
    for col in numeric:
        result[col] = _number(row.get(col))
    if all(result[f] is not None and result[f] > 0 for f in features["gmv"]):
        result["gmv_total_ml"] = sum(result[f] for f in features["gmv"])
    lo, hi = config["analysis"]["age_range"]
    age_ok = result["age"] is not None and lo <= result["age"] <= hi
    cov_ok = age_ok and result["sex"] in (0, 1) and result["icv_ml"] is not None and result["icv_ml"] > 0
    for modality in ("wmh", "gmv"):
        valid = all(result[f] is not None and
                    (result[f] >= 0 if modality == "wmh" else result[f] > 0) for f in features[modality])
        result[f"{modality}_ready"] = bool(cov_ok and valid)
        if not valid:
            issues.append(f"{modality}_missing_or_invalid")
    if not age_ok:
        issues.append("outside_reference_age_or_missing_age")
    if not cov_ok:
        issues.append("reference_covariates_unavailable")
    result["reference_eligible"] = bool(cov_ok)
    return result, {"subject_key": row["subject_key"], "issues": issues,
                    "wmh_ready": result["wmh_ready"], "gmv_ready": result["gmv_ready"]}
