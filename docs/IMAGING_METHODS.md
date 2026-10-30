# Imaging engine implementation

The R engine calls GAMLSS/gamlss.dist, FactoMineR, missMDA and clue. It does not substitute a custom PCA or reference-distribution solver. Reference fits save training design data and spline basis for reproducible new-cohort predictions. All identifiers are read as strings.

## Normative actions

`norm_fit` requires `inputs.train`, `feature`, `modality` (`wmh`/`gmv`). Training cohort must be LS. Outputs are `model.rds`, `calibration.tsv`, `metrics.json`, `calibration.png`. WMH requires the configured positive count (also checked within validation folds); all-zero references produce `not_estimable`. ZAGA models the point mass at zero plus conditional Gamma; all-positive references use GA. GMV uses NO on log-volume. Age spline bases are fixed at training. Values outside the configured healthy age interval, missing covariates and invalid volumes are not projected. Prediction-interval level is configured, and output columns are `lower_interval`, `upper_interval`, with `coverage_interval` and `interval_level` in metrics.

`norm_apply` requires `inputs.data` and `inputs.model`. It writes `scores.tsv` with `<feature>_dev`, `<feature>_percentile`, `<feature>_valid`, `<feature>_cdf_clipped`. WMH zeros use the midpoint CDF (nu/2), so the resulting signed score is not exactly standard normal for a zero-inflated distribution. GMV deficit is negative conventional Z. Normative calibration is evaluated and reported; completion does not automatically certify calibration adequacy.

## PCA actions

`pca_fit`: `inputs.train`, explicit ordered `features`, `modality`, and `k` (integer or `max`), or `inputs.selection`. `max` fits the maximum permissible reference rank and does not perform dimension selection. `pca_reference` is an alias. Outputs include model, scores, loadings, eigenvalues, `component_metadata.tsv` and scree/loading PNGs. The component metadata has `component`, `index`, `group_id`, `group_size`, `individually_identifiable` and `score_sd`; downstream interaction caps must not cut a group. No missing-value imputation is allowed. Model scales and axes are learned in C3 only.

`pca_cv`: same training inputs, optionally `replicate`. With a replicate index it executes one missMDA Kfold masked-cell repetition, providing a checkpoint per repetition. Without it the configured number of repeats is executed. Each K uses the same reset seed. `criterion.tsv` contains all candidate K values, convergence flags and diagnostics. The upstream criterion is converted from squared-error sum to masked-cell mean squared error. Reaching the upstream EM iteration limit produces an invalid candidate (`msep=NA`, `converged=FALSE`), not a successful numerical estimate. This is matrix-cell CV, not subject-level clinical validation. The upstream API does not expose its internal EM iteration limit.

`pca_bootstrap`: `inputs.train`, `inputs.reference` (a maximal-rank model), `k`, `replicate`, `modality`. Produces `bootstrap.tsv`. It matches loading columns with Hungarian assignment, aligns signs, computes largest subspace angle and absolute-agreement ICC on out-of-bag subject scores. Near-degenerate components are assessed as subspaces; single-axis criteria do not gate those components. Candidates that cut through a near-degenerate group cannot be selected.

`pca_select`: `inputs.cv` is a list of criterion TSV paths, `inputs.bootstrap` a list of bootstrap TSV paths. It requires all configured repeats for all candidate dimensions. It writes `selection.tsv` and `selection.json`. Missing validation does not silently force K=0. Among stable candidates, choose the lowest MSEP, then the smaller K. K=0 is a valid mean-only candidate.

`pca_project`: `inputs.data`, `inputs.model`. Score columns are `WMH_PC01`, `GMV_PC01` or `ACUTE_PC01`; additional `_C3SD` columns use frozen C3 standard deviations. `*_projection_valid` flags missing features. Positive sign orientation is fixed by the largest absolute loading. GB never refits centers, scales or loadings.

## Joint and acute-location actions

`joint` accepts `inputs.wmh_scores`, `inputs.gmv_scores` and optional `inputs.model`. Without a model, both inputs must be C3. Each retained block is divided by the square root of its total score covariance trace; training scales are saved. `JOINT_WMH_PC01` and `JOINT_GMV_PC01` remain separate coordinates. A zero-dimensional block contributes no coordinates. Both blocks K=0 yields no spatial joint phenotype, not an artificial zero-risk score.

`acute_residualize` accepts `inputs.train`, `features`, optionally `inputs.data`; application accepts `inputs.model` and `inputs.data`. Input requires `acute_volume_ml` and `age`. Both fitting and projection enforce the configured reference age interval; missing age columns return `waiting_input`, and out-of-range/missing age values yield invalid projection rows. The training age interval is stored in the model. Models lacking this provenance require refitting. Regional proportions are regressed on a fixed natural-spline basis of log1p(volume), with C3 coefficients saved. Outputs `residuals.tsv` and `features.json` feed acute PCA. A negative residual means less injury at that location than predicted for overall lesion size; it does not mean negative injury.

Python `extract_acute` records the source mask's declared template, original-grid voxel volume, atlas file/dictionary hashes and uncovered lesion fraction. Both lesion and atlas NIfTI headers must explicitly declare `mm`, matching `analysis.acute.expected_spatial_unit`. Metre, micron and unknown units are rejected, never silently reinterpreted. Its `acute_volume_ml` is the volume in the supplied mask's coordinate space; a normalized mask must not be described as a native-space volume measurement. Grid resampling only aligns sampling grids after matching declared templates. Probability load uses direct weighted sums, not NiftiMapsMasker regression coefficients.

## Verification

`Rscript tests/R/test_imaging.R` uses synthetic data and deliberately reduced validation repeats. It tests reference folds, zero WMH, missing/out-of-range projection, unchanged frozen models, masked CV, bootstrap selection, K=0, block scaling and acute-volume residualization. A separate `Rscript --vanilla` process reloads the acute spline model and verifies frozen predictions, including required S3 namespace registration. Passing it is software verification, not cohort validation.
