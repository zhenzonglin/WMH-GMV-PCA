"""Clinical field audit and lossless handoff; no outcome/exposure imputation."""
from __future__ import annotations

import calendar
import csv
import datetime as dt
import json
import math
from pathlib import Path
from typing import Any


BASELINE = ("prestroke_mrs", "nihss", "prior_stroke", "hypertension", "diabetes", "af", "ivt", "evt")
BINARY = {"prior_stroke", "hypertension", "diabetes", "af", "ivt", "evt", "sex",
          "any_ht", "symptomatic_ht", "recanalized", "cmb", "css"}


def _read(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="," if path.suffix.lower() == ".csv" else "\t")
        return list(reader), list(reader.fieldnames or [])


def _write(path: Path, rows: list[dict[str, Any]], columns: list[str] | None = None) -> None:
    columns = columns or list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _number(value: Any) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def _date(value: Any) -> dt.date | None:
    if not value:
        return None
    try:
        return dt.date.fromisoformat(str(value).strip()[:10])
    except ValueError:
        return None


def _month_date(origin: dt.date, months: int) -> dt.date:
    month_index = origin.year * 12 + origin.month - 1 + months
    year, month = divmod(month_index, 12)
    return dt.date(year, month + 1, min(origin.day, calendar.monthrange(year, month + 1)[1]))


def _valid(canonical: str, number: float | None) -> bool:
    if number is None:
        return False
    if canonical.startswith("mrs_m"):
        return number.is_integer() and 0 <= number <= 6
    if canonical == "prestroke_mrs":
        return number.is_integer() and 0 <= number <= 5
    if canonical == "nihss":
        return number.is_integer() and 0 <= number <= 42
    if canonical in BINARY or canonical.startswith("death_m"):
        return number in (0, 1)
    return True


def prepare_clinical(imaging_table: Path, clinical_table: Path | None, config: dict,
                     output_dir: Path, cohort: str) -> dict:
    """Join explicit IDs, audit every configured field and preserve missingness."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    image_rows, image_columns = _read(Path(imaging_table))
    settings = config["analysis"]["clinical"]
    mapping = config["variables"].get("clinical", {})
    fields = mapping.get("fields", {})
    missing_codes = {str(v) for v in mapping.get("missing_codes", [])} | {""}
    visits = sorted(int(v) for v in settings["visits"])
    configured_covariates = settings.get("covariates", ["age", "sex", *BASELINE])
    clinical_covariates = [name for name in configured_covariates if name not in {"age", "sex", "icv_ml", "wmh_total_ml", "gmv_total_ml"}]
    source_rows: list[dict[str, str]] = []
    source_columns: list[str] = []
    if clinical_table and Path(clinical_table).is_file():
        source_rows, source_columns = _read(Path(clinical_table))
    # 别名只用于查找列；多列同时存在且值不一致时必须记录冲突。
    def sources(name: str) -> list[str]:
        value = fields.get(name, name)
        values = value if isinstance(value, list) else [value]
        return [v for v in values if v in source_columns]
    id_columns = sources("participant_id")
    cohort_columns = sources("cohort")
    source_index: dict[str, dict[str, str]] = {}
    duplicates: set[str] = set()
    for row in source_rows:
        if cohort_columns and row.get(cohort_columns[0], "") != cohort:
            continue
        identifier = row.get(id_columns[0], "") if id_columns else ""
        if not identifier:
            continue
        if identifier in source_index:
            duplicates.add(identifier)
        source_index[identifier] = row

    canonical_fields = list(dict.fromkeys([*BASELINE, *clinical_covariates, "index_date", "death_date",
        *(f"{prefix}_m{v}" for v in visits for prefix in ("mrs", "death", "visit_date")),
        *(name for name in fields if name not in {"participant_id", "cohort", "age", "sex", "icv_ml", "wmh_total_ml", "gmv_total_ml"})]))
    audit = {name: {"variable": name, "source_columns": "|".join(sources(name)),
                   "column_present": bool(sources(name)), "available": 0,
                   "missing": 0, "invalid": 0, "total": len(image_rows)} for name in canonical_fields}
    conflicts: list[dict[str, Any]] = []
    wide: list[dict[str, Any]] = []
    long: list[dict[str, Any]] = []
    events: list[dict[str, Any]] = []
    medications: list[dict[str, Any]] = []
    eligibility: list[dict[str, Any]] = []
    for image in image_rows:
        identifier = image["participant_id"]
        key = image.get("subject_key", f"{cohort}:{identifier}")
        source = source_index.get(identifier, {}) if identifier not in duplicates else {}
        row: dict[str, Any] = dict(image)
        row["clinical_joined"] = int(bool(source))
        row["clinical_conflict"] = 0
        def conflict(variable: str, raw: Any, reason: str) -> None:
            row["clinical_conflict"] = 1
            conflicts.append({"subject_key": key, "variable": variable, "raw_value": raw, "reason": reason})
        if identifier in duplicates:
            conflict("participant_id", identifier, "duplicate_clinical_id")
        for name in canonical_fields:
            local_missing = {str(v) for v in mapping.get("field_missing_codes", {}).get(name, missing_codes)} | {""}
            all_values = [(c, source.get(c, "").strip()) for c in sources(name)]
            nonmissing = [(c,v) for c,v in all_values if v not in local_missing]
            def translate(column, value):
                rules = mapping.get("source_value_maps", {}).get(column,
                    mapping.get("value_maps", {}).get(name, {}))
                return rules.get(value, value)
            comparable = []
            for column, value in nonmissing:
                translated = translate(column, value)
                number = _number(translated)
                comparable.append(str(number) if number is not None else str(translated))
            if len(set(comparable)) > 1:
                conflict(name, "|".join(f"{c}={v}" for c, v in nonmissing), "conflicting_alias_columns_after_code_mapping")
            source_column, raw = nonmissing[0] if nonmissing else next(((c,v) for c,v in all_values if v), ("", ""))
            row[f"source_{name}"] = raw
            row[f"source_column_{name}"] = source_column
            row[name] = ""
            if not raw or raw in local_missing:
                audit[name]["missing"] += 1
                continue
            converted = translate(source_column, raw)
            if converted is None or str(converted) in local_missing:
                audit[name]["missing"] += 1
                continue
            if "date" in name:
                value = _date(converted)
                valid = value is not None
                normalized: Any = value.isoformat() if value else ""
            else:
                value_num = _number(converted)
                numeric_field = name in BASELINE or name.startswith(("mrs_m", "death_m")) or name in BINARY
                valid = _valid(name, value_num) if numeric_field else True
                normalized = value_num if numeric_field else converted
            if valid:
                row[name] = normalized
                audit[name]["available"] += 1
            else:
                audit[name]["invalid"] += 1
        row["cmb_assessment_status"] = ("not_assessable" if row.get("source_column_cmb") == "IMG_SVD_CMB" and row.get("source_cmb") == "99"
            else "available" if _number(row.get("cmb")) in (0,1) else "not_recorded_or_invalid")
        death_date = _date(row.get("death_date"))
        origin = _date(row.get("index_date"))
        if death_date and origin and death_date < origin:
            conflict("death_date", death_date, "death_before_index_date")
            death_date = None
        dead_by_previous = False
        for month in visits:
            mrs_name, death_name = f"mrs_m{month}", f"death_m{month}"
            visit_date = _date(row.get(f"visit_date_m{month}"))
            if visit_date is None and origin and settings.get("use_nominal_visit_dates_for_death", False):
                visit_date = _month_date(origin, month)
            flag = _number(row.get(death_name))
            mrs = _number(row.get(mrs_name))
            confirmed = dead_by_previous or flag == 1 or bool(death_date and visit_date and death_date <= visit_date)
            if confirmed:
                if flag == 0 or (mrs is not None and mrs < 6):
                    conflict(mrs_name, row.get(f"source_{mrs_name}"), "alive_score_or_flag_after_confirmed_death")
                row[mrs_name] = 6
                row[f"mrs_origin_m{month}"] = "confirmed_death"
                dead_by_previous = True
            elif mrs == 6:
                # mRS=6 本身为死亡记录；其来源和与显式生存记录的矛盾保留。
                if flag == 0 or (death_date and visit_date and death_date > visit_date):
                    conflict(mrs_name, mrs, "mrs_death_conflicts_with_alive_evidence")
                row[f"mrs_origin_m{month}"] = "reported_mrs_death"
                dead_by_previous = True
            else:
                row[f"mrs_origin_m{month}"] = "observed" if mrs is not None else "missing"
        required = list(dict.fromkeys(["age", "sex", "icv_ml", "wmh_total_ml", "gmv_total_ml", *configured_covariates]))
        missing_baseline = [name for name in required if _number(row.get(name)) is None]
        age = _number(row.get("age"))
        age_min, age_max = config["analysis"]["age_range"]
        row["clinical_age_eligible"] = int(age is not None and age_min <= age <= age_max)
        conflict_excluded = bool(settings.get("exclude_conflicts", True) and row["clinical_conflict"])
        reasons = missing_baseline + (["clinical_conflict"] if conflict_excluded else []) + ([] if row["clinical_age_eligible"] else ["age_outside_reference_range"])
        last_visit = max(visits)
        b1_single = not reasons and _number(row.get(f"mrs_m{last_visit}")) is not None
        b1_long = not reasons and any(_number(row.get(f"mrs_m{v}")) is not None for v in visits)
        b2 = b1_single and _number(row.get("acute_volume_ml")) is not None and any(c.startswith("ACUTE_PC") and _number(row.get(c)) is not None for c in row)
        row["eligible_b1_single"] = int(b1_single)
        row["eligible_b1_longitudinal"] = int(b1_long)
        row["eligible_b2"] = int(b2)
        eligibility.append({"subject_key": key, "b1_single": int(b1_single), "b1_longitudinal": int(b1_long),
                            "b2": int(b2), "missing_baseline": "|".join(reasons),
                            "missing_mrs12": int(_number(row.get(f"mrs_m{last_visit}")) is None),
                            "clinical_joined": row["clinical_joined"]})
        for month in visits:
            item = dict(row)
            item.update({"visit": month, "mrs": row.get(f"mrs_m{month}", ""),
                         "visit_date": row.get(f"visit_date_m{month}", ""),
                         "mrs_origin": row[f"mrs_origin_m{month}"]})
            long.append(item)
        for spec in mapping.get("events", []):
            confirmed = source.get(spec.get("confirmed", ""), "")
            date = source.get(spec.get("date", ""), "")
            unit = spec.get("time_unit", "calendar")
            relative_days = _number(date) if unit == "days" else None
            rules = mapping.get("source_value_maps", {}).get(spec.get("confirmed", ""), {})
            event_flag = _number(rules.get(confirmed, confirmed))
            # 保留记录缺失和原始代码，不将空白转换为“无事件”。
            events.append({"subject_key": key, "event_type": spec["type"], "confirmed_raw": confirmed,
                           "event_confirmed": event_flag if event_flag in (0,1) else "",
                           "event_time_raw": date, "time_unit": unit, "time_origin": spec.get("time_origin", "calendar"),
                           "relative_days": relative_days if relative_days is not None and relative_days >= 0 else "",
                           "time_valid": int(relative_days is not None and relative_days >= 0) if unit == "days" else int(_date(date) is not None),
                           "event_date_raw": date if unit == "calendar" else "", "date_valid": int(unit == "calendar" and _date(date) is not None),
                           "source_confirmed": spec.get("confirmed", ""), "source_date": spec.get("date", "")})
        for spec in mapping.get("medications", []):
            medications.append({"subject_key": key, "visit_month": spec["visit_month"], "drug": spec["drug"],
                                "status_raw": source.get(spec.get("status", ""), ""),
                                "dose_raw": source.get(spec.get("dose", ""), ""), "exposure_type": "snapshot"})
        wide.append(row)
    _write(output_dir / "clinical_wide.tsv", wide)
    _write(output_dir / "clinical_long.tsv", long)
    _write(output_dir / "variable_audit.tsv", list(audit.values()))
    _write(output_dir / "eligibility.tsv", eligibility)
    _write(output_dir / "conflicts.tsv", conflicts, ["subject_key", "variable", "raw_value", "reason"])
    _write(output_dir / "events.tsv", events, ["subject_key", "event_type", "confirmed_raw", "event_confirmed",
        "event_time_raw", "time_unit", "time_origin", "relative_days", "time_valid", "event_date_raw", "date_valid", "source_confirmed", "source_date"])
    _write(output_dir / "medication_snapshots.tsv", medications, ["subject_key", "visit_month", "drug", "status_raw", "dose_raw", "exposure_type"])
    missing_required = [name for name in clinical_covariates if not audit[name]["column_present"]]
    if not source_rows or not id_columns:
        missing_required.append("clinical_table_or_identifier")
    result = {"status": "completed", "reason": "Clinical audit completed; model readiness is reported separately.",
              "n_subjects": len(wide), "model_ready": not missing_required,
              "missing_required_columns": missing_required,
              "eligible_b1_single": sum(r["b1_single"] for r in eligibility),
              "eligible_b1_longitudinal": sum(r["b1_longitudinal"] for r in eligibility),
              "eligible_b2": sum(r["b2"] for r in eligibility),
              "artifacts": ["clinical_wide.tsv", "clinical_long.tsv", "variable_audit.tsv", "eligibility.tsv", "conflicts.tsv", "events.tsv", "medication_snapshots.tsv"]}
    (output_dir / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return result
