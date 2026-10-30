import json
import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from wmh_gmv.validation import _merge_scores, _r, _CheckpointStop, run_internal_validation


class InternalValidationTests(unittest.TestCase):
    def test_runtime_change_reuses_scientific_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            config={"project":{"project_root":str(root)},"paths":{"rscript":"Rscript"},
                    "analysis":{"seed":1},"variables":{},
                    "runtime":{"execution":{"blas_threads_per_job":1,"cores":4}}}
            def execute(command,**kwargs):
                request=json.loads(Path(command[-1]).read_text())
                out=Path(request["output_dir"])
                (out/"value.txt").write_text("verified")
                (out/"result.json").write_text(json.dumps({"status":"completed","artifacts":["value.txt"]}))
            with patch("wmh_gmv.validation.subprocess.run",side_effect=execute) as runner:
                _r(config,root/"checkpoint",action="test",inputs={})
                adjusted=copy.deepcopy(config);adjusted["runtime"]["execution"]["cores"]=8
                _r(adjusted,root/"checkpoint",action="test",inputs={})
                self.assertEqual(runner.call_count,1)
                adjusted["analysis"]["seed"]=2
                _r(adjusted,root/"checkpoint",action="test",inputs={})
                self.assertEqual(runner.call_count,2)

    def test_stop_before_subprocess_preserves_checkpoint_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);stop=root/"stop_request.json";stop.write_text("{}")
            with self.assertRaises(_CheckpointStop):
                _r({"_validation_stop_request":str(stop)},root/"new",action="pca_fit")
            (root/"checkpoint.json").write_text("{}")
            with patch("wmh_gmv.validation._run_internal_validation",side_effect=_CheckpointStop("stop")):
                result=run_internal_validation({}, {"run_dir":str(root)},root)
            self.assertEqual(result["status"],"interrupted")
            self.assertIn("checkpoint.json",result["artifacts"])

    def test_duplicate_score_rows_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"scores.tsv"
            pd.DataFrame({"subject_key":["C3:1","C3:1"],"WMH_PC01":[1,2]}).to_csv(path,sep="\t",index=False)
            with self.assertRaises(pd.errors.MergeError):
                _merge_scores(pd.DataFrame({"subject_key":["C3:1"]}),[path])

    def test_missing_inputs_not_fake_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            result=run_internal_validation({"analysis":{"clinical":{}}},{},Path(directory))
            self.assertEqual(result["status"],"waiting_input")


if __name__=="__main__": unittest.main()
