# RINENG impurity partitioning

This repository contains the non-confidential data and code used for the study
**“Impurity partitioning during forced-cooling desupersaturation of industrial
wet-process phosphoric acid: a chemometric and mineralogical study.”**

The workflow reproduces the profile-based k-nearest-neighbor (kNN) imputation,
the PCA and Ward.D2 HCA analyses, the sensitivity and leave-one-time-out (LOTO)
assessments, and the data-driven main and supplementary figures.

## Scope and limitations

- Seven sampling times were studied: 0, 12, 24, 36, 48, 60, and 72 h.
- Two independent experimental replicates were used at each time; reported
  uncertainty is the sample standard deviation (n = 2).
- All aliquots originated from one industrial batch. Batch-to-batch variability
  was not evaluated.
- `D_i(t)` denotes the raw apparent solid–liquid partition coefficient.
- `Z_i(t)` denotes the centered and scaled temporal profile used for PCA and HCA.
- Co and Pb were excluded from PCA and HCA because reliable temporal partition
  coefficients could not be established.
- The multivariate results describe statistical associations and profile
  similarities; they do not demonstrate causality or common retention mechanisms.

## Repository structure

```text
data/raw/            Input partition-coefficient workbook
data/external/       Journal workbooks supplied locally and not tracked by Git
data/processed/      Imputed matrices and sensitivity scenarios
notebooks/           Validated profile-based kNN notebook
scripts/             Command-line workflow and figure-generation scripts
results/tables/      PCA, HCA, sensitivity, and LOTO output tables
results/figures/     Reference PDF figures and regenerated graphical outputs
documentation/       Data dictionary and validation notes
```

## Installation

The complete environment can be created with Conda:

```bash
conda env create -f environment.yml
conda activate rineng-impurity-partitioning
```

For Python-only installation:

```bash
python -m pip install -r requirements.txt
```

## Reproduction order

Run all commands from the repository root.

```bash
python scripts/01_profile_knn_imputation.py
Rscript scripts/02_pca_hca_analysis.R
Rscript scripts/03_generate_main_figures.R
Rscript scripts/04_generate_supplementary_figures.R
Rscript scripts/05_capture_session_info.R
```

Before generating Figures 4–7, place local copies of the two workbooks supplied
to the journal at the following paths:

```text
data/external/Supplementary_Data.xlsx
data/external/Granulo_laser_30.xlsx
```

These files are intentionally excluded from GitHub. Their expected SHA-256
checksums are documented in `documentation/external_files.md`.

The notebook `notebooks/01_profile_knn_imputation.ipynb` contains the same
Python analysis in an interactive form. Its saved outputs were cleared so that
local paths and stale execution state are not distributed.

## Validated primary results

- Selected imputation: uniform weighting, `k_NN = 8`.
- Standardized cross-validation RMSE: 0.8008221.
- Standardized cross-validation MAE: 0.6208004.
- PCA: PC1 = 61.1496%, PC2 = 24.0283%, cumulative PC1–PC2 = 85.1780%.
- PC3 = 10.5964%; cumulative PC1–PC3 = 95.7744%.
- HCA: Euclidean distance, Ward.D2 agglomeration, selected solution `G = 3`.

The `G = 3` solution was retained as a parsimonious compromise between the
elbow criterion, average and individual silhouette widths, dendrogram structure,
and group-size balance. It was not selected by maximizing the mean silhouette.

## Data availability

Only non-confidential inputs needed for the core kNN/PCA/HCA analyses are
intended for public release. `Supplementary_Data.xlsx` and the PSD workbook are
supplied to the journal but are not stored on GitHub. Private elemental
mass-balance files and industrially sensitive information are not included.
Raw SEM, XRD, and FTIR instrument files are outside the present repository.

## License

- Python and R source code: MIT License; see `LICENSE`.
- Data and documentation: Creative Commons Attribution 4.0 International
  (CC BY 4.0); see `LICENSE-DATA.md`.

Attribution should cite the associated article and this repository.

## Reproducibility validation

The original Python workflow reports Python 3.13.5. It was also independently
rerun under Python 3.12.13, and all numerical CSV outputs matched the supplied
validated outputs. The complete R workflow was then executed without errors
from an empty global environment under R 4.5.1 on Windows 11. Exact R package
versions are recorded in `sessionInfo_R.txt`; further details are provided in
`documentation/validation_summary.md`.
