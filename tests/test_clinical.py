import csv
import tempfile
import unittest
from pathlib import Path

from wmh_gmv.clinical import prepare_clinical


def write_table(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def read_table(path):
    with path.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class ClinicalPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.config = {"analysis": {"age_range": [50,75], "clinical": {"visits": [3, 6, 12], "exclude_conflicts": True,
                       "use_nominal_visit_dates_for_death": False}},
                       "variables": {"clinical": {"fields": {"participant_id": "code_n"},
                       "missing_codes": ["99"], "value_maps": {},
                       "medications": [{"visit_month": 3, "drug": "aspirin", "status": "aspirin"}],
                       "events": [{"type": "ICH", "confirmed": "ich", "date": "ich_date"}]}}}
        self.image = self.root / "image.tsv"
        self.clinical = self.root / "clinical.tsv"
        self.out = self.root / "output"
        write_table(self.image, [{"participant_id": "0007", "subject_key": "C3:0007", "cohort": "C3",
             "age": 60, "sex": 1, "icv_ml": 1400, "wmh_total_ml": 6, "gmv_total_ml": 600,
             "WMH_PC01": 0.3, "GMV_PC01": -0.2, "ACUTE_PC01": 0.1, "acute_volume_ml": 2}])

    def tearDown(self):
        self.temp.cleanup()

    def clinical_row(self, **changes):
        row = {"code_n": "0007", "prestroke_mrs": 0, "nihss": 5, "prior_stroke": 0,
               "hypertension": 1, "diabetes": 0, "af": 0, "ivt": 0, "evt": 0,
               "mrs_m3": "", "mrs_m6": "", "mrs_m12": "", "death_m3": "",
               "death_m6": "", "death_m12": "", "death_date": "", "visit_date_m3": "",
               "visit_date_m6": "", "visit_date_m12": "", "aspirin": "", "ich": "", "ich_date": ""}
        row.update(changes)
        return row

    def run_prepare(self, rows):
        write_table(self.clinical, rows)
        return prepare_clinical(self.image, self.clinical, self.config, self.out, "C3")

    def test_death_propagates_only_forward_and_id_preserved(self):
        self.run_prepare([self.clinical_row(death_m6=1)])
        row = read_table(self.out / "clinical_wide.tsv")[0]
        self.assertEqual(row["participant_id"], "0007")
        self.assertEqual(row["mrs_m3"], "")
        self.assertEqual(row["mrs_m6"], "6")
        self.assertEqual(row["mrs_m12"], "6")
        self.assertEqual(row["eligible_b1_single"], "1")

    def test_known_date_needs_visit_time_and_does_not_fill_earlier(self):
        self.run_prepare([self.clinical_row(death_date="2025-05-01",visit_date_m3="2025-04-01",visit_date_m6="2025-07-01")])
        row = read_table(self.out / "clinical_wide.tsv")[0]
        self.assertEqual(row["mrs_m3"], "")
        self.assertEqual(row["mrs_m6"], "6")
        self.assertEqual(row["mrs_m12"], "6")

    def test_conflicting_alive_score_is_flagged_and_not_silently_eligible(self):
        result = self.run_prepare([self.clinical_row(death_m3=1,mrs_m6=2)])
        self.assertEqual(result["eligible_b1_single"], 0)
        self.assertTrue(read_table(self.out / "conflicts.tsv"))

    def test_missing_event_and_medication_are_not_negative_exposure(self):
        self.run_prepare([self.clinical_row(mrs_m12=2)])
        self.assertEqual(read_table(self.out / "events.tsv")[0]["confirmed_raw"], "")
        medication = read_table(self.out / "medication_snapshots.tsv")[0]
        self.assertEqual(medication["status_raw"], "")
        self.assertEqual(medication["exposure_type"], "snapshot")

    def test_invalid_outcome_and_duplicates_not_accepted(self):
        self.run_prepare([self.clinical_row(mrs_m12=7)])
        audit = {r["variable"]: r for r in read_table(self.out / "variable_audit.tsv")}
        self.assertEqual(audit["mrs_m12"]["invalid"], "1")
        self.run_prepare([self.clinical_row(mrs_m12=2), self.clinical_row(mrs_m12=4)])
        row = read_table(self.out / "clinical_wide.tsv")[0]
        self.assertEqual(row["clinical_joined"], "0")
        self.assertEqual(row["eligible_b1_single"], "0")

    def test_unavailable_clinical_file_still_produces_audit(self):
        result = prepare_clinical(self.image, None, self.config, self.out, "C3")
        self.assertEqual(result["status"], "completed")
        self.assertFalse(result["model_ready"])
        self.assertIn("clinical_table_or_identifier", result["missing_required_columns"])


if __name__ == "__main__":
    unittest.main()
