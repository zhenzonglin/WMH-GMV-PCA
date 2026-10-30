# Project rules

- Read the approved SAP and current config before changing scientific settings.
- LS is the healthy reference, C3 is discovery/training, GB is frozen representation validation. Formal age range is 50–75.
- Accept existing corrected WMH; do not add GM lesion exclusions, center harmonization/covariates, MRI-onset covariates, SuStaIn, brain age or conditional mismatch.
- Use mature GAMLSS/FactoMineR/missMDA/ordinal implementations. Never replace unavailable libraries with undocumented surrogate algorithms.
- Keep paths and tunable settings in config. Record upstream fixed implementation limits instead of inventing nonfunctional knobs.
- Do not use clinical outcomes to select imaging components. Refit all C3-dependent transformations inside internal prediction folds.
- Keep failed, waiting, non-estimable and nonsignificant findings distinct. Never alter thresholds to obtain clinical significance.
- Do not start real cohort analyses without a user instruction. Synthetic validation must be labeled and use isolated paths.
- Never broadly kill processes. Use the verified stop entrypoint; preserve finished artifacts and failure records.
- Do not commit patient tables, raw images, fitted patient models, logs or private machine configurations.
- Preserve code/config/input/software identities, and maintain implementation and analysis records. Screenshot-based interpretations must state their evidence level.
