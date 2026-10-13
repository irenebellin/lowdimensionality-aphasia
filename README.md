# Code for: "A low dimensional structure of behavioral impairment and anatomical damage in aphasia syndromes"

Bellin I., Bisogno A.L., De Pellegrin S., Facchini S., Vallesi A., Corbetta M.
(manuscript BRAINCOM-2026-026, Brain Communications)

This repository contains the analysis code. Patient data (imaging, clinical scores) are not
included for privacy reasons; see "Data availability" in the manuscript for how to request them.

## Contents

| File | What it does | Manuscript section |
|---|---|---|
| `matlab/ridge_pca_permutation.m` | PCA of the imaging matrix, ridge regression with leave-one-out cross-validation, permutation test, optional multiple-comparison correction | Methods, ridge regression |
| `matlab/load_maps_matrix.m` | Reads one NIfTI per patient into a patients x voxels matrix | Methods |
| `matlab/weights_to_nifti.m` | Back-projects the weights to the brain, scales to [-1, 1], thresholds at 0.05, smooths (`smooth3`, 5 voxels), saves NIfTI | Methods, Supplementary S1 |
| `matlab/run_ridge_lesion.m` | Example runner for lesion maps | Figure 4A |
| `matlab/run_ridge_disconnection.m` | Example runner for disconnection maps | Figure 4B |
| `r/latent_profile_analysis.R` | Latent profile analysis (Gaussian mixture models, `mclust`), Sankey diagrams cohort/syndrome to cluster, radar plots of the cluster profiles | Methods and Results, latent profile analysis; Supplementary S5 |
| `python/outlier_analysis.py` | Mahalanobis-distance outlier analysis of the Pure Aphasia Cohort against the General Stroke Cohort (PC1-PC3 plane and three-dimensional space) | Methods and Results, outlier analysis |
| `python/sensitivity_analysis.py` | Smallest detectable effect (power 80%) for the cohort and syndrome-group comparisons | Methods, sample size and sensitivity analysis |
| `python/cohort_comparison_from_summary.py` | Comparison of the two cohorts from the summary statistics of Tables 1 and 2 | Results, cohort comparison |

## Requirements

- MATLAB with the Statistics and Machine Learning Toolbox and the Image Processing Toolbox
  (`pca`, `zscore`, `niftiread`/`niftiwrite`, `smooth3`).
- R 4.x with the packages `readxl`, `mclust`, `networkD3`, `htmlwidgets`, `jsonlite`, `fmsb` and `scales` (installed automatically by the script).
- Python 3.9+ with the packages in `python/requirements.txt` (`pip install -r python/requirements.txt`).

## Ridge regression (MATLAB)

1. Put one NIfTI per patient in a folder (binary lesion masks, or disconnection maps, same space,
   MNI152 2 mm). Create `scores.csv` with a column `ID` and a column with the behavioral score
   (e.g. `PC1`). Each file name must start with the patient ID.
2. Edit the paths at the top of `run_ridge_lesion.m` or `run_ridge_disconnection.m` and run it.

Settings used in the manuscript: components explaining at least 95% of the variance, z-scored
component scores and behavioral score, ridge penalty chosen from 200 values (`logspace(-2, 4, 200)`)
by leave-one-out cross-validation (maximum cross-validated R2), 10,000 permutations of the
behavioral scores, two-tailed alpha = 0.05.

**Multiple comparisons.** The option `opts.correction` controls how significant components are
selected:

- `'none'`: each component is compared with its own 2.5th/97.5th permutation percentiles. This is
  the original analysis, which applies no correction across components.
- `'maxT'` (default): family-wise error control across components, using the maximum absolute
  coefficient in each permutation.
- `'fdr'`: Benjamini-Hochberg false discovery rate on per-component permutation p-values.

State in the manuscript which option produced the reported maps.

**Reported R2.** `res.r2Loo` is the leave-one-out cross-validated R2. The ridge penalty is selected
with the same cross-validation, so the value is slightly optimistic. `res.pR2` is the permutation
p-value of the cross-validated R2.

## Latent profile analysis (R)

Edit the settings block at the top of `r/latent_profile_analysis.R` (paths, column positions of the 30
behavioral variables, the patient and variable that are excluded, and the number of clusters tried)
and run it. Behavioral scores of both cohorts are inverted where needed (higher = better), z-scored
with the mean and SD of the General Stroke Cohort, and clustered together. Outputs (BIC plot, cluster
means, assignments, Sankey diagrams, radar plots) are written to `lpa_output/`.

Notes: `missing = "zero"` reproduces the original analysis, in which missing values were set to 0;
`missing = "drop"` excludes patients with missing values. The script warns when the best BIC is at the
largest number of clusters tried. The adjusted Rand index between clinical syndromes and clusters is an
optional addition that quantifies what the Sankey diagrams show.

## Additional analyses (Python)

All three scripts run from the command line.

```
python python/sensitivity_analysis.py
python python/cohort_comparison_from_summary.py
python python/outlier_analysis.py general_stroke_scores.txt pure_aphasia_scores.txt
```

`outlier_analysis.py` needs two plain-text files with one row per patient and three columns (PC1, PC2,
PC3), the scores used for Figure 2. Distances are computed with the covariance of the General
Stroke Cohort, and the cutoff is the 95th percentile of that cohort's own distances.

## License

Add a license file before making the repository public (for example MIT).
