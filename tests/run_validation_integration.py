"""Explicit slow validation smoke: actual R backends, reduced repeat counts."""
from pathlib import Path
import tempfile
import numpy as np
import pandas as pd

from wmh_gmv.config import load_config
from wmh_gmv.validation import run_internal_validation

root=Path(__file__).resolve().parents[1]
cfg=load_config(root/"config/local/smoke/project.yaml")
cfg["analysis"]["pca"].update(max_components=1,cv_repeats=1,bootstrap_repeats=1)
cfg["analysis"]["clinical"]["prediction_folds"]=2
directory=Path(tempfile.mkdtemp(prefix="wmh-validation-integration-"))
metadata=pd.read_csv(root/"data/synthetic/C3/metadata.tsv",sep="\t",dtype={"participant_id":str})
clinical=pd.read_csv(root/"data/synthetic/C3/clinical.tsv",sep="\t",dtype={"participant_id":str})
data=metadata.merge(clinical,on="participant_id",validate="one_to_one")
data["subject_key"]="C3:"+data.participant_id;data["cohort"]="C3"
# Poison columns must disappear before each training fold.
data["WMH_PC99"]=data.mrs_m12*1000
data["GMV_PC99"]=data.mrs_m12*1000
data["JOINT_leaky"]=data.mrs_m12*1000
clinical_path=directory/"clinical.tsv";data.to_csv(clinical_path,sep="\t",index=False)
rng=np.random.default_rng(244)
inputs={"clinical_wide":str(clinical_path)}
for modality in ("wmh","gmv"):
    latent=rng.normal(size=len(data))
    matrix=data[["subject_key","participant_id","cohort"]].copy()
    for index in range(4):matrix[f"{modality}_{index}_dev"]=latent*(index+1)+rng.normal(0,.01,len(data))
    path=directory/f"{modality}.tsv";matrix.to_csv(path,sep="\t",index=False)
    inputs[modality+"_deviations"]=str(path)
result=run_internal_validation(cfg,inputs,directory/"output")
assert result["status"]=="completed",result
pred=pd.read_csv(directory/"output/heldout_predictions.tsv",sep="\t")
assert not pred.duplicated(["subject_key","model"]).any()
assert set(pred.validation)=={"C3_outer_fold_heldout"}
for path in (directory/"output").glob("fold_*/clinical_train.tsv"):
    columns=pd.read_csv(path,sep="\t",nrows=0).columns
    assert "WMH_PC99" not in columns and "GMV_PC99" not in columns and "JOINT_leaky" not in columns
print("Actual-backend patient-fold validation and leakage checks passed:",directory)
