# Reproducibility validation summary

## Python workflow

The final notebook was rerun from the supplied raw partition-coefficient
workbook in a clean working directory under Python 3.12.13. All 13 exported CSV
files matched the supplied validated outputs numerically; the maximum absolute
numerical difference was zero at a comparison tolerance of 1 × 10⁻¹⁴. The
original notebook execution reports Python 3.13.5.

The primary scenario reproduced:

- uniform weighting and `k_NN = 8`;
- standardized RMSE = 0.800822060733485;
- standardized MAE = 0.6208003622584358;
- the six imputed values listed in the data dictionary.

## R workflow

The complete cleaned R workflow was executed from the repository root on
Windows 11 using R 4.5.1 after restarting R and clearing the global environment
(`length(ls(all.names = TRUE)) = 0`). The PCA/HCA analysis, main-figure runner,
supplementary-figure runner, and session-information capture all completed
without errors. The reproduced HCA solution used `G = 3`, and the cophenetic
correlation was 0.8671. Exact package versions and locale information are
recorded in the root-level `sessionInfo_R.txt`.

Non-fatal warnings were limited to packages built under later R 4.5 patch
releases and deprecation messages originating within `factoextra`. They did not
alter the numerical calculations. The final analysis script uses a registered
Times New Roman alias on Windows and a portable serif fallback elsewhere.

## External figure inputs

The supplied PSD workbook contains the expected `Granulo-S30%` sheet and the
input ranges used by Figure 6. It is excluded from GitHub at the authors'
request, together with `Supplementary_Data.xlsx`; both files are supplied to
the journal. Their checksums are recorded in `external_files.md`.
