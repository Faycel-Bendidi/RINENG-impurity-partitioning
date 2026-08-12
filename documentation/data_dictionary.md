# Data dictionary

## Core matrices

| File | Content |
|---|---|
| `partition_coefficients_raw.xlsx` | Seven-time matrix of apparent solid–liquid partition coefficients before imputation; Co and Pb are retained in the source sheet but excluded from PCA/HCA. |
| `partition_coefficients_imputed.csv` | Completed 7 × 17 matrix used in the primary analysis. |
| `partition_coefficients_standardized_Z.csv` | Variable-wise centered and scaled profiles used for PCA and HCA. |
| `Supplementary_Data.xlsx` | Journal-supplied external workbook used by Figures 4, 5, and 7; not stored on GitHub. |
| `Granulo_laser_30.xlsx` | Journal-supplied external PSD workbook used by Figure 6; not stored on GitHub. |

## Core notation

| Symbol | Definition |
|---|---|
| `D_i(t)` | Apparent solid–liquid partition coefficient, `C_i,solid(t) / C_i,liquid(t)`. |
| `Z_i(t)` | Standardized temporal profile, `[D_i(t) − mean(D_i)] / s_Di`, using the sample standard deviation across the seven times. |
| `k_NN` | Number of neighbors in the profile-based kNN procedure. |
| `G` | Number of HCA groups. |

## Six imputed values in the primary scenario

| Variable | Time (h) | Imputed `D_i(t)` |
|---|---:|---:|
| K2O | 0 | 0.5430803471253334 |
| Cu | 12 | 0.0268319415227704 |
| Cu | 24 | 0.0248813085668851 |
| Mo | 0 | 0.3851771868618981 |
| Mo | 24 | 0.3408442326745187 |
| Ni | 0 | 0.0762133499118981 |

Imputed values have no experimental error bars because they were estimated by
the validated model rather than obtained from replicated measurements.
