"""Patient-level internal validation with fold-specific imaging refitting."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


class _CheckpointStop(Exception):
    """Explicit cooperative stop after a completed R operation."""


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def _r(config: dict, directory: Path, **request) -> dict:
    stop_request = config.get("_validation_stop_request")
    if stop_request and Path(stop_request).exists():
        raise _CheckpointStop("Stopped at internal validation checkpoint")
    directory.mkdir(parents=True, exist_ok=True)
    root = Path(config["project"]["project_root"])
    request = {**request, "config": config, "output_dir": str(directory)}
    fingerprints = {}
    for name,value in request.get("inputs", {}).items():
        for index,item in enumerate(value if isinstance(value, list) else [value]):
            if isinstance(item, str) and Path(item).is_file(): fingerprints[f"{name}:{index}"] = _sha(Path(item))
    for name in ("common.R", "normative.R", "pca.R", "clinical.R", "run_stage.R"):
        path = root / "R" / name
        if path.is_file(): fingerprints[f"R/{name}"] = _sha(path)
    for name in ("renv.lock","uv.lock"):
        if (root / name).is_file(): fingerprints[name] = _sha(root / name)
    # 资源和监测设置变化不能使已完成的科学重复失效；输入身份用内容哈希。
    scientific_request = {k:v for k,v in request.items() if k not in {"config","inputs","output_dir"}}
    scientific_request["analysis"] = config["analysis"]
    scientific_request["variables"] = config.get("variables",{})
    identity = hashlib.sha256(json.dumps({"request":scientific_request,"inputs":fingerprints},sort_keys=True).encode()).hexdigest()
    receipt = directory / "checkpoint.json"
    if receipt.is_file():
        saved = json.loads(receipt.read_text())
        if saved.get("identity") == identity and all((directory / name).is_file() and _sha(directory / name)==value for name,value in saved.get("artifacts",{}).items()):
            return saved["result"]
    _json(directory / "request.json", request)
    environment = dict(os.environ)
    if config["paths"].get("r_library"):
        environment["R_LIBS"] = config["paths"]["r_library"] + os.pathsep + environment.get("R_LIBS", "")
    # run_stage.R 将项目库放在 .libPaths 首位，同时保留锁文件记录的基础库。
    environment["OPENBLAS_NUM_THREADS"] = str(config["runtime"]["execution"]["blas_threads_per_job"])
    with (directory / "R.log").open("w",encoding="utf-8") as log:
        subprocess.run([config["paths"]["rscript"],str(root / "R" / "run_stage.R"),str(directory / "request.json")],
                       cwd=root,env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
    result = json.loads((directory / "result.json").read_text())
    artifacts = result.get("artifacts",[])
    if isinstance(artifacts,str): artifacts = [artifacts]
    hashes = {name:_sha(directory / name) for name in [*artifacts,"result.json"]}
    _json(receipt,{"identity":identity,"result":result,"artifacts":hashes})
    return result


def _read(path: str) -> pd.DataFrame:
    return pd.read_csv(path,sep="\t",dtype=str,keep_default_na=False)


def _write(data: pd.DataFrame, path: Path) -> str:
    data.to_csv(path,sep="\t",index=False)
    return str(path)


def _merge_scores(data: pd.DataFrame, tables: list[Path]) -> pd.DataFrame:
    for path in tables:
        scores = _read(str(path))
        columns = [c for c in scores if c.startswith(("WMH_PC","GMV_PC"))]
        data = data.merge(scores[["subject_key",*columns]],on="subject_key",how="left",validate="one_to_one")
    return data


def _run_internal_validation(config: dict, inputs: dict, output_dir: Path) -> dict:
    """Refit all C3-dependent representation in each heldout patient fold.

    Frozen LS normative models are intentionally reused. Every R operation has a
    checksum checkpoint; no full-sample PC or JOINT column can enter a fold.
    """
    output_dir = Path(output_dir); output_dir.mkdir(parents=True,exist_ok=True)
    cfg = config["analysis"]["clinical"]
    def finish(status: str, reason: str, artifacts: list[str], **extra) -> dict:
        result = {"status":status,"reason":reason,"artifacts":artifacts,**extra}
        _json(output_dir / "result.json",result);return result
    needed = ["wmh_deviations","gmv_deviations","clinical_wide"]
    absent = [name for name in needed if name not in inputs or not Path(inputs[name]).is_file()]
    if absent: return finish("waiting_input","Missing internal-validation inputs: "+", ".join(absent),[])
    clinical = _read(inputs["clinical_wide"])
    clinical = clinical.drop(columns=[c for c in clinical if c.startswith(("WMH_PC","GMV_PC","ACUTE_PC","JOINT_"))])
    required = [*cfg["covariates"],"icv_ml","wmh_total_ml","gmv_total_ml",f"mrs_m{max(cfg['visits'])}"]
    if any(c not in clinical for c in required):
        return finish("waiting_input","Missing clinical covariates or outcome",[])
    ages = pd.to_numeric(clinical["age"],errors="coerce")
    clinical = clinical[ages.between(*config["analysis"]["age_range"])].copy().reset_index(drop=True)
    if clinical.subject_key.duplicated().any(): raise ValueError("Repeated clinical ID in internal validation")
    n_folds = int(cfg["prediction_folds"])
    if n_folds < 2: raise ValueError("Internal validation requires at least two patient folds")
    if len(clinical) < n_folds: return finish("not_estimable","Too few patients for requested folds",[])
    rng = np.random.default_rng(config["analysis"]["seed"])
    perm = rng.permutation(len(clinical)); fold_id = np.empty(len(clinical),dtype=int)
    fold_id[perm] = np.arange(len(clinical)) % n_folds + 1
    assignment = pd.DataFrame({"subject_key":clinical.subject_key,"fold":fold_id})
    _write(assignment,output_dir / "fold_assignments.tsv")
    matrices = {name:_read(inputs[name+"_deviations"]) for name in ("wmh","gmv")}
    predictions, statuses = [], []
    selected_folds = [int(inputs["fold"])] if "fold" in inputs else list(range(1,n_folds+1))
    for fold in selected_folds:
        fold_dir = output_dir / f"fold_{fold:02d}";fold_dir.mkdir(exist_ok=True)
        _json(output_dir / "progress.json",{"fold":fold,"total_folds":n_folds,"stage":"fold_imaging_refit","completed_folds":len(statuses)})
        train_keys = set(assignment.loc[assignment.fold!=fold,"subject_key"])
        test_keys = set(assignment.loc[assignment.fold==fold,"subject_key"])
        training = clinical[clinical.subject_key.isin(train_keys)].copy()
        heldout = clinical[clinical.subject_key.isin(test_keys)].copy()
        training_scores, heldout_scores = [], []
        blocked = ""
        for modality, matrix in matrices.items():
            features = [c for c in matrix if c.endswith("_dev")]
            if not features:
                blocked = f"No estimable {modality} deviation features";break
            directory = fold_dir / modality;directory.mkdir(exist_ok=True)
            train_path = _write(matrix[matrix.subject_key.isin(train_keys)],directory / "train.tsv")
            test_path = _write(matrix[matrix.subject_key.isin(test_keys)],directory / "heldout.tsv")
            common = {"modality":modality,"features":features}
            reference_dir = directory / "reference"
            reference = _r(config,reference_dir,action="pca_fit",inputs={"train":train_path},k="max",**common)
            if reference["status"]!="completed": blocked = f"{modality} reference: {reference.get('reason','')}";break
            maximum = int(reference["metrics"]["k"])
            cv_paths, bootstrap_paths = [], []
            for repeat in range(1,config["analysis"]["pca"]["cv_repeats"]+1):
                cv_dir = directory / f"cv_{repeat:04d}"
                result = _r(config,cv_dir,action="pca_cv",inputs={"train":train_path},replicate=repeat,**common)
                if result["status"]!="completed": blocked="Incomplete PCA CV";break
                cv_paths.append(str(cv_dir / "criterion.tsv"))
            if blocked: break
            for k in range(1,maximum+1):
                for repeat in range(1,config["analysis"]["pca"]["bootstrap_repeats"]+1):
                    _json(output_dir / "progress.json",{"fold":fold,"modality":modality,"stage":"bootstrap","k":k,"replicate":repeat,
                         "total_replicates":config["analysis"]["pca"]["bootstrap_repeats"],"completed_folds":len(statuses)})
                    boot_dir = directory / f"bootstrap_k{k:02d}_{repeat:04d}"
                    result = _r(config,boot_dir,action="pca_bootstrap",inputs={"train":train_path,"reference":str(reference_dir / "model.rds")},k=k,replicate=repeat,**common)
                    if result["status"]!="completed": blocked="Incomplete PCA bootstrap";break
                    bootstrap_paths.append(str(boot_dir / "bootstrap.tsv"))
                if blocked: break
            if blocked: break
            select_dir = directory / "select"
            result = _r(config,select_dir,action="pca_select",inputs={"cv":cv_paths,"bootstrap":bootstrap_paths},**common)
            if result["status"]!="completed":blocked="PCA selection unavailable";break
            fit_dir = directory / "fit"
            result = _r(config,fit_dir,action="pca_fit",inputs={"train":train_path,"selection":str(select_dir / "selection.json")},**common)
            if result["status"]!="completed": blocked="PCA fit unavailable";break
            project_dir = directory / "heldout"
            result = _r(config,project_dir,action="pca_project",inputs={"data":test_path,"model":str(fit_dir / "model.rds")},**common)
            if result["status"]!="completed":blocked="PCA projection unavailable";break
            training_scores.append(fit_dir / "scores.tsv");heldout_scores.append(project_dir / "scores.tsv")
        if blocked:
            statuses.append({"fold":fold,"status":"not_estimable","reason":blocked,"n_train":len(training),"n_test":len(heldout)})
            continue
        training = _merge_scores(training,training_scores);heldout = _merge_scores(heldout,heldout_scores)
        train_path = _write(training,fold_dir / "clinical_train.tsv")
        test_path = _write(heldout,fold_dir / "clinical_heldout.tsv")
        fit_dir = fold_dir / "clinical_fit"
        result = _r(config,fit_dir,action="clinical_b1_single",inputs={"data":train_path},cohort="C3")
        if result["status"]!="completed":
            statuses.append({"fold":fold,"status":result["status"],"reason":result.get("reason",""),"n_train":len(training),"n_test":len(heldout)})
            continue
        prediction_dir = fold_dir / "clinical_heldout"
        result = _r(config,prediction_dir,action="clinical_predict",inputs={"data":test_path,"reference_models":str(fit_dir / "models.rds")},
                    cohort="C3",validation_label="C3_outer_fold_heldout")
        if result["status"]=="completed":
            frame = pd.read_csv(prediction_dir / "predictions.tsv",sep="\t",dtype={"subject_key":str})
            frame["fold"] = fold;predictions.append(frame)
        statuses.append({"fold":fold,"status":result["status"],"reason":result.get("reason",""),"n_train":len(training),"n_test":len(heldout)})
    _write(pd.DataFrame(statuses),output_dir / "fold_status.tsv")
    artifacts = ["fold_assignments.tsv","fold_status.tsv","progress.json"]
    if not predictions: return finish("not_estimable","No heldout predictions could be estimated",artifacts,n_subjects=len(clinical))
    combined = pd.concat(predictions,ignore_index=True)
    if combined.duplicated(["subject_key","model"]).any(): raise ValueError("Patient appears in more than one heldout fold")
    _write(combined,output_dir / "heldout_predictions.tsv");artifacts.append("heldout_predictions.tsv")
    metrics, calibration = [], []
    for model, group in combined.groupby("model",sort=False):
        probability = group[[f"p_mrs{i}" for i in range(7)]].to_numpy(float)
        y = group.observed_mrs.to_numpy(int)
        if not np.allclose(probability.sum(axis=1),1): raise ValueError("Invalid heldout probability totals")
        cumulative = probability.cumsum(axis=1)[:,:6]
        target = y[:,None] <= np.arange(6)[None,:]
        metrics.extend([{"model":model,"metric":"ranked_probability_score","value":float(np.mean((cumulative-target)**2)),"n":len(group)},
            {"model":model,"metric":"ordinal_log_loss","value":float(-np.log(np.maximum(cfg["prediction_epsilon"],probability[np.arange(len(y)),y])).mean()),"n":len(group)}])
        for threshold in range(6):
            p = cumulative[:,threshold];bins=np.minimum(cfg["calibration_bins"]-1,(p*cfg["calibration_bins"]).astype(int))
            for bin_id in np.unique(bins):
                selected=bins==bin_id
                calibration.append({"model":model,"threshold":threshold,"bin":int(bin_id),"n":int(selected.sum()),
                    "predicted":float(p[selected].mean()),"observed":float((y[selected]<=threshold).mean())})
    _write(pd.DataFrame(metrics),output_dir / "heldout_scores.tsv")
    _write(pd.DataFrame(calibration),output_dir / "heldout_calibration.tsv")
    artifacts.extend(["heldout_scores.tsv","heldout_calibration.tsv"])
    all_done = len(statuses)==n_folds and all(r["status"]=="completed" for r in statuses)
    return finish("completed" if all_done else "not_estimable",
        "All patient-heldout folds completed" if all_done else "Partial heldout results only; one or more folds unavailable",
        artifacts,n_subjects=len(clinical),validation="C3_patient_outer_folds",all_folds_completed=all_done)


def run_internal_validation(config: dict, inputs: dict, output_dir: Path) -> dict:
    local_config = dict(config)
    if inputs.get("run_dir"):
        local_config["_validation_stop_request"] = str(Path(inputs["run_dir"]) / "stop_request.json")
    try:
        return _run_internal_validation(local_config, inputs, output_dir)
    except _CheckpointStop as error:
        checkpoint_files = sorted(str(path.relative_to(output_dir)) for path in Path(output_dir).rglob("*")
                                  if path.is_file() and path != Path(output_dir) / "result.json" and not path.name.endswith(".tmp"))
        result = {"status":"interrupted","reason":str(error),"artifacts":checkpoint_files,
                  "checkpoint_directory":str(output_dir)}
        _json(Path(output_dir) / "result.json",result)
        return result
