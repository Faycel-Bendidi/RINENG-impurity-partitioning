"""Command-line export of the validated profile-based kNN notebook.

Run this script from the repository root. It uses the same code cells as the
notebook and writes all processed matrices and diagnostics to data/processed.
"""

import matplotlib
matplotlib.use("Agg")

try:
    from IPython.display import display
except ImportError:
    def display(value):
        print(value)

# Step 1 — Import libraries and define reproducibility settings

from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error

# Reproducibility setting for subsequent analyses
RANDOM_SEED = 2026
np.random.seed(RANDOM_SEED)

# Define input and output paths
WORKING_DIRECTORY = Path.cwd()
DATA_PATH = WORKING_DIRECTORY / "data" / "raw" / "partition_coefficients_raw.xlsx"
RESULTS_DIRECTORY = WORKING_DIRECTORY / "data" / "processed"

RESULTS_DIRECTORY.mkdir(parents=True, exist_ok=True)

print("Python version:", sys.version.split()[0])
print("Operating system:", platform.system())
print("Working directory:", WORKING_DIRECTORY)
print("Input file found:", DATA_PATH.exists())
print("Results directory:", RESULTS_DIRECTORY)

# Step 2 — Load and prepare the raw partition-coefficient matrix

SHEET_NAME = "Partition_coefficients"

raw_data = pd.read_excel(
    DATA_PATH,
    sheet_name=SHEET_NAME,
    na_values=["N.D.", "N.D", "ND", ""]
)

# Clean column names
raw_data.columns = [
    str(column).strip().replace("\n", " ")
    for column in raw_data.columns
]

# Variables excluded because complete partition profiles were unavailable
EXCLUDED_VARIABLES = ["DCo", "DPb"]

missing_columns = [
    variable
    for variable in EXCLUDED_VARIABLES
    if variable not in raw_data.columns
]

if missing_columns:
    raise ValueError(
        f"Expected excluded variables were not found: {missing_columns}"
    )

partition_data_raw = (
    raw_data
    .drop(columns=EXCLUDED_VARIABLES)
    .copy()
)

partition_data_raw["Time_h"] = pd.to_numeric(
    partition_data_raw["Time_h"],
    errors="raise"
)

partition_data_raw = (
    partition_data_raw
    .set_index("Time_h")
    .apply(pd.to_numeric, errors="raise")
)

print(
    "Raw analysis matrix dimensions:",
    f"{partition_data_raw.shape[0]} time points × "
    f"{partition_data_raw.shape[1]} variables"
)

display(partition_data_raw)

# Step 3 — Validate time points, variables, and missing entries

EXPECTED_TIME_POINTS = [0, 12, 24, 36, 48, 60, 72]

EXPECTED_MISSING_POSITIONS = {
    (0, "DK2O"),
    (12, "DCu"),
    (24, "DCu"),
    (0, "DMo"),
    (24, "DMo"),
    (0, "DNi"),
}

actual_time_points = partition_data_raw.index.astype(int).tolist()

if actual_time_points != EXPECTED_TIME_POINTS:
    raise ValueError(
        f"Unexpected time points: {actual_time_points}"
    )

if partition_data_raw.shape != (7, 17):
    raise ValueError(
        "The analysis matrix must contain 7 time points and 17 variables."
    )

actual_missing_positions = {
    (int(time_h), variable)
    for time_h in partition_data_raw.index
    for variable in partition_data_raw.columns
    if pd.isna(partition_data_raw.loc[time_h, variable])
}

if actual_missing_positions != EXPECTED_MISSING_POSITIONS:
    raise ValueError(
        "Unexpected missing-value pattern.\n"
        f"Expected: {sorted(EXPECTED_MISSING_POSITIONS)}\n"
        f"Found: {sorted(actual_missing_positions)}"
    )

missing_entries = pd.DataFrame(
    [
        {
            "Time_h": time_h,
            "Variable": variable
        }
        for time_h, variable in sorted(
            actual_missing_positions,
            key=lambda item: (item[0], item[1])
        )
    ]
)

total_entries = partition_data_raw.size
number_of_missing_entries = len(actual_missing_positions)
missing_percentage = (
    100 * number_of_missing_entries / total_entries
)

print("Data validation completed successfully.")
print(f"Total matrix entries: {total_entries}")
print(f"Missing entries: {number_of_missing_entries}")
print(f"Missing-data percentage: {missing_percentage:.2f}%")

display(missing_entries)

# Step 4 — Define the profile-based k-nearest-neighbor imputation function

def profile_knn_predict(
    data,
    target_variable,
    missing_time,
    k,
    weighting="uniform",
    minimum_common_time_points=4
):
    """
    Predict one missing partition-coefficient value using neighboring
    autoscaled temporal profiles.

    Each target profile is standardized using its available observations.
    Each donor profile is standardized using all of its available
    observations, including the value at the sampling time being imputed.
    """

    if target_variable not in data.columns:
        raise KeyError(f"Unknown target variable: {target_variable}")

    if missing_time not in data.index:
        raise KeyError(f"Unknown sampling time: {missing_time}")

    if pd.notna(data.loc[missing_time, target_variable]):
        raise ValueError(
            f"{target_variable} at {missing_time} h is not missing."
        )

    if weighting not in {"uniform", "distance"}:
        raise ValueError(
            "weighting must be either 'uniform' or 'distance'."
        )

    target_profile = data[target_variable].astype(float)
    target_observed_times = target_profile.dropna().index

    target_mean = target_profile.loc[
        target_observed_times
    ].mean()

    target_std = target_profile.loc[
        target_observed_times
    ].std(ddof=1)

    if not np.isfinite(target_std) or target_std == 0:
        raise ValueError(
            f"The observed profile of {target_variable} has "
            "an invalid standard deviation."
        )

    target_z = (
        target_profile - target_mean
    ) / target_std

    donor_records = []

    for donor_variable in data.columns:

        if donor_variable == target_variable:
            continue

        donor_profile = data[donor_variable].astype(float)
        donor_observed_times = donor_profile.dropna().index

        # The donor must have an observed value at the missing time
        if missing_time not in donor_observed_times:
            continue

        common_times = target_observed_times.intersection(
            donor_observed_times
        )

        if len(common_times) < minimum_common_time_points:
            continue

        # Standardize the donor using all of its observed values
        donor_mean = donor_profile.loc[
            donor_observed_times
        ].mean()

        donor_std = donor_profile.loc[
            donor_observed_times
        ].std(ddof=1)

        if not np.isfinite(donor_std) or donor_std == 0:
            continue

        donor_z = (
            donor_profile - donor_mean
        ) / donor_std

        distance = np.sqrt(
            np.sum(
                (
                    target_z.loc[common_times]
                    - donor_z.loc[common_times]
                ) ** 2
            )
        )

        donor_records.append(
            {
                "Donor_variable": donor_variable,
                "Distance": float(distance),
                "Common_time_points": len(common_times),
                "Donor_z_at_missing_time": float(
                    donor_z.loc[missing_time]
                )
            }
        )

    donor_table = (
        pd.DataFrame(donor_records)
        .sort_values(
            ["Distance", "Donor_variable"],
            ascending=[True, True]
        )
        .reset_index(drop=True)
    )

    if donor_table.empty:
        raise ValueError(
            f"No eligible donor profiles were found for "
            f"{target_variable} at {missing_time} h."
        )

    if k > len(donor_table):
        raise ValueError(
            f"k_NN={k} exceeds the number of eligible donors "
            f"({len(donor_table)}) for {target_variable} "
            f"at {missing_time} h."
        )

    selected_donors = donor_table.head(k).copy()

    donor_z_values = selected_donors[
        "Donor_z_at_missing_time"
    ].to_numpy()

    if weighting == "uniform":
        weights = np.ones(k, dtype=float)
    else:
        distances = selected_donors["Distance"].to_numpy()
        weights = 1.0 / np.maximum(distances, 1e-12)

    predicted_z = np.average(
        donor_z_values,
        weights=weights
    )

    predicted_value = (
        target_mean
        + predicted_z * target_std
    )

    selected_donors["Weight"] = (
        weights / weights.sum()
    )

    return {
        "Target_variable": target_variable,
        "Time_h": missing_time,
        "k_NN": k,
        "Weighting": weighting,
        "Predicted_z_score": float(predicted_z),
        "Predicted_value": float(predicted_value),
        "Eligible_donor_count": len(donor_table),
        "Selected_donors": selected_donors,
        "All_eligible_donors": donor_table
    }

# Step 5 — Determine the maximum common candidate value of k

AFFECTED_VARIABLES = ["DK2O", "DCu", "DMo", "DNi"]
MINIMUM_COMMON_TIME_POINTS = 4

donor_count_records = []

# Evaluate the six genuinely missing entries
for missing_time, target_variable in sorted(
    EXPECTED_MISSING_POSITIONS,
    key=lambda item: (item[1], item[0])
):
    result = profile_knn_predict(
        data=partition_data_raw,
        target_variable=target_variable,
        missing_time=missing_time,
        k=1,
        weighting="uniform",
        minimum_common_time_points=MINIMUM_COMMON_TIME_POINTS
    )

    donor_count_records.append(
        {
            "Validation_type": "Actual_missing_entry",
            "Target_variable": target_variable,
            "Time_h": missing_time,
            "Eligible_donor_count": result[
                "Eligible_donor_count"
            ]
        }
    )

# Evaluate each leave-one-observed-value-out validation fold
for target_variable in AFFECTED_VARIABLES:

    observed_times = (
        partition_data_raw[target_variable]
        .dropna()
        .index
    )

    for validation_time in observed_times:

        validation_data = partition_data_raw.copy()

        actual_value = validation_data.loc[
            validation_time,
            target_variable
        ]

        validation_data.loc[
            validation_time,
            target_variable
        ] = np.nan

        result = profile_knn_predict(
            data=validation_data,
            target_variable=target_variable,
            missing_time=validation_time,
            k=1,
            weighting="uniform",
            minimum_common_time_points=MINIMUM_COMMON_TIME_POINTS
        )

        donor_count_records.append(
            {
                "Validation_type": "Artificially_masked_entry",
                "Target_variable": target_variable,
                "Time_h": validation_time,
                "Actual_value": actual_value,
                "Eligible_donor_count": result[
                    "Eligible_donor_count"
                ]
            }
        )

donor_count_table = pd.DataFrame(donor_count_records)

maximum_common_k = int(
    donor_count_table["Eligible_donor_count"].min()
)

K_VALUES = list(range(1, maximum_common_k + 1))

print("Eligible donor-profile counts:")
display(
    donor_count_table[
        [
            "Validation_type",
            "Target_variable",
            "Time_h",
            "Eligible_donor_count"
        ]
    ]
)

print(
    "Minimum number of eligible donors across all cases:",
    maximum_common_k
)

print("Candidate k_NN values:", K_VALUES)

# Step 6 — Summarize eligible donor counts

donor_count_summary = (
    donor_count_table
    .groupby(
        ["Validation_type", "Target_variable"]
    )["Eligible_donor_count"]
    .agg(["min", "max", "mean"])
    .reset_index()
)

print("Summary of eligible donor counts:")
display(donor_count_summary.round(2))

# Step 7 — Perform leave-one-observed-value-out cross-validation

WEIGHTING_METHODS = ["uniform", "distance"]

cross_validation_records = []

for target_variable in AFFECTED_VARIABLES:

    observed_values = (
        partition_data_raw[target_variable]
        .dropna()
        .astype(float)
    )

    # Fixed scale used only to compare prediction errors across variables
    reference_std = observed_values.std(ddof=1)

    if not np.isfinite(reference_std) or reference_std == 0:
        raise ValueError(
            f"Invalid reference standard deviation for {target_variable}."
        )

    for validation_time, actual_value in observed_values.items():

        validation_data = partition_data_raw.copy()

        # Artificially mask one genuinely observed value
        validation_data.loc[
            validation_time,
            target_variable
        ] = np.nan

        for weighting_method in WEIGHTING_METHODS:

            for k in K_VALUES:

                prediction = profile_knn_predict(
                    data=validation_data,
                    target_variable=target_variable,
                    missing_time=validation_time,
                    k=k,
                    weighting=weighting_method,
                    minimum_common_time_points=MINIMUM_COMMON_TIME_POINTS
                )

                predicted_value = prediction["Predicted_value"]
                prediction_error = predicted_value - actual_value
                standardized_error = prediction_error / reference_std

                selected_donors = prediction["Selected_donors"]

                cross_validation_records.append(
                    {
                        "Target_variable": target_variable,
                        "Time_h": validation_time,
                        "Weighting": weighting_method,
                        "k_NN": k,
                        "Actual_value": actual_value,
                        "Predicted_value": predicted_value,
                        "Error": prediction_error,
                        "Absolute_error": abs(prediction_error),
                        "Standardized_error": standardized_error,
                        "Absolute_standardized_error": abs(
                            standardized_error
                        ),
                        "Eligible_donor_count": prediction[
                            "Eligible_donor_count"
                        ],
                        "Selected_donors": ", ".join(
                            selected_donors["Donor_variable"].tolist()
                        )
                    }
                )

cross_validation_results = pd.DataFrame(
    cross_validation_records
)

print(
    "Total number of cross-validation predictions:",
    len(cross_validation_results)
)

display(cross_validation_results.head())

# Step 8 — Summarize global cross-validation performance

def summarize_prediction_errors(group):
    standardized_errors = (
        group["Standardized_error"]
        .dropna()
        .to_numpy()
    )

    return pd.Series(
        {
            "Number_of_predictions": len(standardized_errors),
            "MAE_SD": np.mean(
                np.abs(standardized_errors)
            ),
            "RMSE_SD": np.sqrt(
                np.mean(
                    np.square(standardized_errors)
                )
            ),
            "Median_absolute_error_SD": np.median(
                np.abs(standardized_errors)
            ),
            "Maximum_absolute_error_SD": np.max(
                np.abs(standardized_errors)
            )
        }
    )


summary_records = []

for (weighting_method, k_value), group in (
    cross_validation_results
    .groupby(["Weighting", "k_NN"])
):
    metrics = summarize_prediction_errors(group)

    summary_records.append(
        {
            "Weighting": weighting_method,
            "k_NN": k_value,
            **metrics.to_dict()
        }
    )

cross_validation_summary = (
    pd.DataFrame(summary_records)
    .sort_values(
        ["RMSE_SD", "MAE_SD"]
    )
    .reset_index(drop=True)
)

print("Global cross-validation performance:")
display(cross_validation_summary.round(4))

# Step 9 — Summarize cross-validation performance by variable

variable_summary_records = []

for (
    target_variable,
    weighting_method,
    k_value
), group in cross_validation_results.groupby(
    ["Target_variable", "Weighting", "k_NN"]
):
    metrics = summarize_prediction_errors(group)

    variable_summary_records.append(
        {
            "Target_variable": target_variable,
            "Weighting": weighting_method,
            "k_NN": k_value,
            **metrics.to_dict()
        }
    )

cross_validation_by_variable = pd.DataFrame(
    variable_summary_records
)

best_parameters_by_variable = (
    cross_validation_by_variable
    .sort_values(
        [
            "Target_variable",
            "RMSE_SD",
            "MAE_SD"
        ]
    )
    .groupby(
        "Target_variable",
        as_index=False
    )
    .first()
)

print("Best-performing parameters for each affected variable:")
display(best_parameters_by_variable.round(4))

# Step 10 — Select the global parameter combination

best_global_parameters = (
    cross_validation_summary
    .sort_values(
        ["RMSE_SD", "MAE_SD"],
        ascending=[True, True]
    )
    .iloc[0]
    .copy()
)

SELECTED_WEIGHTING = str(
    best_global_parameters["Weighting"]
)

SELECTED_K = int(
    best_global_parameters["k_NN"]
)

print("Selected global weighting method:", SELECTED_WEIGHTING)
print("Selected global k_NN:", SELECTED_K)
print(
    "Standardized RMSE:",
    round(
        float(best_global_parameters["RMSE_SD"]),
        4
    )
)
print(
    "Standardized MAE:",
    round(
        float(best_global_parameters["MAE_SD"]),
        4
    )
)

# Step 11 — Assess parameter-selection stability by excluding one target variable

parameter_stability_records = []

for excluded_variable in AFFECTED_VARIABLES:

    reduced_results = cross_validation_results[
        cross_validation_results["Target_variable"]
        != excluded_variable
    ].copy()

    reduced_summary_records = []

    for (weighting_method, k_value), group in (
        reduced_results.groupby(["Weighting", "k_NN"])
    ):
        metrics = summarize_prediction_errors(group)

        reduced_summary_records.append(
            {
                "Weighting": weighting_method,
                "k_NN": k_value,
                **metrics.to_dict()
            }
        )

    reduced_summary = (
        pd.DataFrame(reduced_summary_records)
        .sort_values(
            ["RMSE_SD", "MAE_SD"]
        )
        .reset_index(drop=True)
    )

    best_row = reduced_summary.iloc[0]

    parameter_stability_records.append(
        {
            "Excluded_variable": excluded_variable,
            "Selected_weighting": best_row["Weighting"],
            "Selected_k_NN": int(best_row["k_NN"]),
            "RMSE_SD": best_row["RMSE_SD"],
            "MAE_SD": best_row["MAE_SD"]
        }
    )

parameter_selection_stability = pd.DataFrame(
    parameter_stability_records
)

print("Parameter-selection stability after excluding one variable:")
display(parameter_selection_stability.round(4))

# Step 12 — Evaluate zero substitution as a baseline imputation method

# Confirm that the selected k_NN parameters are available
required_objects = [
    "SELECTED_WEIGHTING",
    "SELECTED_K",
    "best_global_parameters"
]

missing_objects = [
    object_name
    for object_name in required_objects
    if object_name not in globals()
]

if missing_objects:
    raise NameError(
        "Run the global parameter-selection step first. "
        f"Missing objects: {missing_objects}"
    )

zero_baseline_records = []

for target_variable in AFFECTED_VARIABLES:

    observed_values = (
        partition_data_raw[target_variable]
        .dropna()
        .astype(float)
    )

    reference_std = observed_values.std(ddof=1)

    if not np.isfinite(reference_std) or reference_std == 0:
        raise ValueError(
            f"Invalid reference standard deviation for "
            f"{target_variable}."
        )

    for time_h, actual_value in observed_values.items():

        # Zero is used only as a baseline prediction for comparison
        predicted_value = 0.0

        prediction_error = predicted_value - actual_value
        standardized_error = prediction_error / reference_std

        zero_baseline_records.append(
            {
                "Target_variable": target_variable,
                "Time_h": time_h,
                "Actual_value": actual_value,
                "Predicted_value": predicted_value,
                "Error": prediction_error,
                "Absolute_error": abs(prediction_error),
                "Standardized_error": standardized_error,
                "Absolute_standardized_error": abs(
                    standardized_error
                )
            }
        )

zero_baseline_results = pd.DataFrame(
    zero_baseline_records
)

zero_standardized_errors = (
    zero_baseline_results["Standardized_error"]
    .to_numpy()
)

zero_baseline_summary = pd.DataFrame(
    [
        {
            "Method": "Zero substitution",
            "Number_of_predictions": len(
                zero_standardized_errors
            ),
            "MAE_SD": np.mean(
                np.abs(zero_standardized_errors)
            ),
            "RMSE_SD": np.sqrt(
                np.mean(
                    np.square(zero_standardized_errors)
                )
            ),
            "Median_absolute_error_SD": np.median(
                np.abs(zero_standardized_errors)
            ),
            "Maximum_absolute_error_SD": np.max(
                np.abs(zero_standardized_errors)
            )
        }
    ]
)

selected_knn_summary = pd.DataFrame(
    [
        {
            "Method": (
                f"Profile k-NN "
                f"({SELECTED_WEIGHTING}, k_NN={SELECTED_K})"
            ),
            "Number_of_predictions": int(
                best_global_parameters[
                    "Number_of_predictions"
                ]
            ),
            "MAE_SD": float(
                best_global_parameters["MAE_SD"]
            ),
            "RMSE_SD": float(
                best_global_parameters["RMSE_SD"]
            ),
            "Median_absolute_error_SD": float(
                best_global_parameters[
                    "Median_absolute_error_SD"
                ]
            ),
            "Maximum_absolute_error_SD": float(
                best_global_parameters[
                    "Maximum_absolute_error_SD"
                ]
            )
        }
    ]
)

baseline_comparison = pd.concat(
    [
        zero_baseline_summary,
        selected_knn_summary
    ],
    ignore_index=True
)

baseline_comparison["RMSE_reduction_vs_zero_percent"] = np.nan
baseline_comparison["MAE_reduction_vs_zero_percent"] = np.nan

zero_rmse = float(
    zero_baseline_summary.loc[0, "RMSE_SD"]
)

zero_mae = float(
    zero_baseline_summary.loc[0, "MAE_SD"]
)

knn_rmse = float(
    selected_knn_summary.loc[0, "RMSE_SD"]
)

knn_mae = float(
    selected_knn_summary.loc[0, "MAE_SD"]
)

baseline_comparison.loc[
    1,
    "RMSE_reduction_vs_zero_percent"
] = 100 * (zero_rmse - knn_rmse) / zero_rmse

baseline_comparison.loc[
    1,
    "MAE_reduction_vs_zero_percent"
] = 100 * (zero_mae - knn_mae) / zero_mae

print("Comparison with zero substitution:")
display(
    baseline_comparison.round(4)
)

# Step 13 — Summarize zero-substitution errors by variable

zero_summary_by_variable_records = []

for target_variable, group in (
    zero_baseline_results
    .groupby("Target_variable")
):

    standardized_errors = (
        group["Standardized_error"]
        .to_numpy()
    )

    zero_summary_by_variable_records.append(
        {
            "Target_variable": target_variable,
            "Number_of_predictions": len(
                standardized_errors
            ),
            "MAE_SD": np.mean(
                np.abs(standardized_errors)
            ),
            "RMSE_SD": np.sqrt(
                np.mean(
                    np.square(standardized_errors)
                )
            ),
            "Median_absolute_error_SD": np.median(
                np.abs(standardized_errors)
            ),
            "Maximum_absolute_error_SD": np.max(
                np.abs(standardized_errors)
            )
        }
    )

zero_summary_by_variable = pd.DataFrame(
    zero_summary_by_variable_records
)

print("Zero-substitution performance by variable:")
display(
    zero_summary_by_variable.round(4)
)

# Step 14 — Impute the six missing entries

# Preserve the raw matrix and create a separate completed matrix
partition_data_imputed = partition_data_raw.copy()

imputation_records = []
selected_donor_details = {}

for missing_time, target_variable in sorted(
    EXPECTED_MISSING_POSITIONS,
    key=lambda item: (item[0], item[1])
):

    prediction = profile_knn_predict(
        data=partition_data_raw,
        target_variable=target_variable,
        missing_time=missing_time,
        k=SELECTED_K,
        weighting=SELECTED_WEIGHTING,
        minimum_common_time_points=MINIMUM_COMMON_TIME_POINTS
    )

    predicted_value = prediction["Predicted_value"]
    selected_donors = prediction["Selected_donors"].copy()

    # Insert the predicted value into the completed matrix
    partition_data_imputed.loc[
        missing_time,
        target_variable
    ] = predicted_value

    # Store detailed donor information for the next analysis step
    selected_donor_details[
        (int(missing_time), target_variable)
    ] = selected_donors

    imputation_records.append(
        {
            "Target_variable": target_variable,
            "Time_h": int(missing_time),
            "Selected_weighting": SELECTED_WEIGHTING,
            "Selected_k_NN": SELECTED_K,
            "Eligible_donor_count": prediction[
                "Eligible_donor_count"
            ],
            "Selected_donor_count": len(selected_donors),
            "Predicted_z_score": prediction[
                "Predicted_z_score"
            ],
            "Imputed_value": predicted_value,
            "Nearest_donor": selected_donors.iloc[0][
                "Donor_variable"
            ],
            "Nearest_donor_distance": selected_donors.iloc[0][
                "Distance"
            ],
            "Nearest_donor_weight": selected_donors.iloc[0][
                "Weight"
            ],
            "Selected_donors": ", ".join(
                selected_donors["Donor_variable"].tolist()
            )
        }
    )

imputed_values_summary = (
    pd.DataFrame(imputation_records)
    .sort_values(
        ["Time_h", "Target_variable"]
    )
    .reset_index(drop=True)
)

# Confirm that all six missing entries were imputed
remaining_missing_entries = int(
    partition_data_imputed.isna().sum().sum()
)

if remaining_missing_entries != 0:
    raise ValueError(
        f"The completed matrix still contains "
        f"{remaining_missing_entries} missing entries."
    )

# Confirm that originally observed values were not modified
observed_mask = partition_data_raw.notna().to_numpy()

raw_observed_values = (
    partition_data_raw
    .to_numpy()[observed_mask]
)

completed_observed_values = (
    partition_data_imputed
    .to_numpy()[observed_mask]
)

if not np.allclose(
    raw_observed_values,
    completed_observed_values
):
    raise ValueError(
        "At least one originally observed value was modified."
    )

print("Selected imputation method:")
print("Weighting:", SELECTED_WEIGHTING)
print("Number of neighbors (k_NN):", SELECTED_K)

print("\nImputed values:")
display(
    imputed_values_summary.round(6)
)

print("\nCompleted partition-coefficient matrix:")
display(
    partition_data_imputed.round(6)
)

# Step 15 — Inspect selected donors, distances, and weights

donor_detail_records = []

for (
    time_h,
    target_variable
), donor_table in selected_donor_details.items():

    for donor_rank, (_, donor_row) in enumerate(
        donor_table.iterrows(),
        start=1
    ):

        donor_detail_records.append(
            {
                "Target_variable": target_variable,
                "Time_h": time_h,
                "Donor_rank": donor_rank,
                "Donor_variable": donor_row[
                    "Donor_variable"
                ],
                "Distance": donor_row["Distance"],
                "Weight": donor_row["Weight"],
                "Common_time_points": donor_row[
                    "Common_time_points"
                ],
                "Donor_z_at_missing_time": donor_row[
                    "Donor_z_at_missing_time"
                ]
            }
        )

selected_donor_table = (
    pd.DataFrame(donor_detail_records)
    .sort_values(
        [
            "Time_h",
            "Target_variable",
            "Donor_rank"
        ]
    )
    .reset_index(drop=True)
)

print("Selected donor profiles for each imputed entry:")
display(
    selected_donor_table.round(6)
)

# Step 16 — Verify donor weights and effective donor counts

donor_weight_diagnostics = (
    selected_donor_table
    .groupby(
        ["Target_variable", "Time_h"],
        as_index=False
    )
    .agg(
        Weight_sum=("Weight", "sum"),
        Maximum_weight=("Weight", "max"),
        Minimum_weight=("Weight", "min"),
        Number_of_selected_donors=(
            "Donor_variable",
            "count"
        )
    )
)

effective_donor_counts = (
    selected_donor_table
    .assign(
        Squared_weight=lambda table: table["Weight"] ** 2
    )
    .groupby(
        ["Target_variable", "Time_h"],
        as_index=False
    )["Squared_weight"]
    .sum()
)

effective_donor_counts[
    "Effective_number_of_donors"
] = (
    1.0
    / effective_donor_counts["Squared_weight"]
)

donor_weight_diagnostics = (
    donor_weight_diagnostics
    .merge(
        effective_donor_counts[
            [
                "Target_variable",
                "Time_h",
                "Effective_number_of_donors"
            ]
        ],
        on=["Target_variable", "Time_h"],
        how="left"
    )
)

donor_weight_diagnostics[
    "Valid_weight_sum"
] = np.isclose(
    donor_weight_diagnostics["Weight_sum"],
    1.0
)

donor_weight_diagnostics[
    "High_weight_concentration"
] = (
    donor_weight_diagnostics["Maximum_weight"]
    > 0.40
)

print("Donor-weight diagnostics:")
display(
    donor_weight_diagnostics.round(6)
)

if not donor_weight_diagnostics[
    "Valid_weight_sum"
].all():
    raise ValueError(
        "At least one set of donor weights does not sum to one."
    )

# Step 17 — Assess imputation sensitivity around the selected k_NN

SENSITIVITY_K_VALUES = sorted(
    {
        max(1, SELECTED_K - 1),
        SELECTED_K,
        min(maximum_common_k, SELECTED_K + 1)
    }
)

local_sensitivity_records = []

for missing_time, target_variable in sorted(
    EXPECTED_MISSING_POSITIONS,
    key=lambda item: (item[0], item[1])
):

    for k_value in SENSITIVITY_K_VALUES:

        prediction = profile_knn_predict(
            data=partition_data_raw,
            target_variable=target_variable,
            missing_time=missing_time,
            k=k_value,
            weighting=SELECTED_WEIGHTING,
            minimum_common_time_points=MINIMUM_COMMON_TIME_POINTS
        )

        selected_donors = prediction["Selected_donors"]

        local_sensitivity_records.append(
            {
                "Target_variable": target_variable,
                "Time_h": int(missing_time),
                "Weighting": SELECTED_WEIGHTING,
                "k_NN": k_value,
                "Imputed_value": prediction["Predicted_value"],
                "Maximum_weight": selected_donors["Weight"].max(),
                "Effective_number_of_donors": (
                    1.0
                    / np.sum(
                        np.square(
                            selected_donors["Weight"].to_numpy()
                        )
                    )
                ),
                "Selected_donors": ", ".join(
                    selected_donors["Donor_variable"].tolist()
                )
            }
        )

local_sensitivity_results = (
    pd.DataFrame(local_sensitivity_records)
    .sort_values(
        ["Time_h", "Target_variable", "k_NN"]
    )
    .reset_index(drop=True)
)

print("Imputed values around the selected k_NN:")
display(
    local_sensitivity_results.round(6)
)

# Summarize the variation in imputed values across k_NN = 7, 8, and 9

selected_k_values = (
    local_sensitivity_results[
        local_sensitivity_results["k_NN"] == SELECTED_K
    ]
    [
        [
            "Target_variable",
            "Time_h",
            "Imputed_value"
        ]
    ]
    .rename(
        columns={
            "Imputed_value": "Selected_k_NN_value"
        }
    )
)

local_sensitivity_summary = (
    local_sensitivity_results
    .groupby(
        ["Target_variable", "Time_h"],
        as_index=False
    )
    .agg(
        Minimum_value=("Imputed_value", "min"),
        Maximum_value=("Imputed_value", "max"),
        Mean_value=("Imputed_value", "mean"),
        Standard_deviation=("Imputed_value", "std")
    )
    .merge(
        selected_k_values,
        on=["Target_variable", "Time_h"],
        how="left"
    )
)

local_sensitivity_summary[
    "Absolute_range"
] = (
    local_sensitivity_summary["Maximum_value"]
    - local_sensitivity_summary["Minimum_value"]
)

local_sensitivity_summary[
    "Relative_range_percent"
] = (
    100
    * local_sensitivity_summary["Absolute_range"]
    / local_sensitivity_summary["Selected_k_NN_value"].abs()
)

print("Local k_NN-sensitivity summary:")
display(
    local_sensitivity_summary.round(6)
)

# Step 18 — Assess leave-one-donor-out influence on imputed values

leave_one_donor_out_records = []

for missing_time, target_variable in sorted(
    EXPECTED_MISSING_POSITIONS,
    key=lambda item: (item[0], item[1])
):

    baseline_prediction = profile_knn_predict(
        data=partition_data_raw,
        target_variable=target_variable,
        missing_time=missing_time,
        k=SELECTED_K,
        weighting=SELECTED_WEIGHTING,
        minimum_common_time_points=MINIMUM_COMMON_TIME_POINTS
    )

    baseline_value = baseline_prediction["Predicted_value"]

    baseline_donors = (
        baseline_prediction["Selected_donors"][
            "Donor_variable"
        ]
        .tolist()
    )

    for excluded_donor in baseline_donors:

        reduced_data = partition_data_raw.drop(
            columns=[excluded_donor]
        )

        reduced_prediction = profile_knn_predict(
            data=reduced_data,
            target_variable=target_variable,
            missing_time=missing_time,
            k=SELECTED_K,
            weighting=SELECTED_WEIGHTING,
            minimum_common_time_points=MINIMUM_COMMON_TIME_POINTS
        )

        reduced_value = reduced_prediction["Predicted_value"]

        absolute_change = abs(
            reduced_value - baseline_value
        )

        relative_change_percent = (
            100
            * absolute_change
            / abs(baseline_value)
            if baseline_value != 0
            else np.nan
        )

        replacement_donors = (
            reduced_prediction["Selected_donors"][
                "Donor_variable"
            ]
            .tolist()
        )

        newly_included_donors = sorted(
            set(replacement_donors)
            - set(baseline_donors)
        )

        leave_one_donor_out_records.append(
            {
                "Target_variable": target_variable,
                "Time_h": int(missing_time),
                "Excluded_donor": excluded_donor,
                "Baseline_value": baseline_value,
                "Recalculated_value": reduced_value,
                "Absolute_change": absolute_change,
                "Relative_change_percent": (
                    relative_change_percent
                ),
                "Newly_included_donor": ", ".join(
                    newly_included_donors
                )
            }
        )

leave_one_donor_out_results = (
    pd.DataFrame(leave_one_donor_out_records)
    .sort_values(
        [
            "Target_variable",
            "Time_h",
            "Relative_change_percent"
        ],
        ascending=[True, True, False]
    )
    .reset_index(drop=True)
)

print("Leave-one-donor-out influence analysis:")
display(
    leave_one_donor_out_results.round(6)
)

# Summarize the maximum donor influence for each imputed entry

donor_influence_summary = (
    leave_one_donor_out_results
    .groupby(
        ["Target_variable", "Time_h"],
        as_index=False
    )
    .agg(
        Most_influential_donor=(
            "Excluded_donor",
            "first"
        ),
        Maximum_relative_change_percent=(
            "Relative_change_percent",
            "first"
        ),
        Maximum_absolute_change=(
            "Absolute_change",
            "first"
        )
    )
)

print("Maximum donor influence for each imputed entry:")
display(
    donor_influence_summary.round(6)
)

# Step 19 — Export the imputed matrix and diagnostic tables for R analysis

export_tables = {
    "partition_coefficients_raw.csv": (
        partition_data_raw.reset_index()
    ),
    "partition_coefficients_imputed.csv": (
        partition_data_imputed.reset_index()
    ),
    "imputed_values_summary.csv": (
        imputed_values_summary
    ),
    "knn_cross_validation_summary.csv": (
        cross_validation_summary
    ),
    "donor_weight_diagnostics.csv": (
        donor_weight_diagnostics
    ),
    "local_k_sensitivity_summary.csv": (
        local_sensitivity_summary
    ),
    "donor_influence_summary.csv": (
        donor_influence_summary
    )
}

for file_name, table in export_tables.items():

    output_path = RESULTS_DIRECTORY / file_name

    table.to_csv(
        output_path,
        index=False
    )

    print(f"Exported: {output_path}")

## Step 19B — Summarize cross-validation performance for all candidate settings

from pathlib import Path
import pandas as pd

summary_records = []

for (weighting_method, k_value), group in cross_validation_results.groupby(
    ["Weighting", "k_NN"]
):
    metrics = summarize_prediction_errors(group)

    summary_records.append(
        {
            "Weighting": weighting_method,
            "k_NN": int(k_value),
            "RMSE_SD": metrics["RMSE_SD"],
            "MAE_SD": metrics["MAE_SD"],
        }
    )

knn_cross_validation_summary = pd.DataFrame(summary_records)

# Use uniform k_NN = 8 as the primary reference
reference_rmse = knn_cross_validation_summary.loc[
    (
        knn_cross_validation_summary["Weighting"].eq("uniform")
        & knn_cross_validation_summary["k_NN"].eq(8)
    ),
    "RMSE_SD",
].iloc[0]

knn_cross_validation_summary[
    "Delta_RMSE_percent_vs_uniform_k8"
] = (
    100
    * (
        knn_cross_validation_summary["RMSE_SD"]
        - reference_rmse
    )
    / reference_rmse
)

# Display uniform weighting before inverse-distance weighting
weighting_order = pd.CategoricalDtype(
    categories=["uniform", "distance"],
    ordered=True,
)

knn_cross_validation_summary["Weighting"] = (
    knn_cross_validation_summary["Weighting"].astype(
        weighting_order
    )
)

knn_cross_validation_summary = (
    knn_cross_validation_summary
    .sort_values(["Weighting", "k_NN"])
    .reset_index(drop=True)
)

print(
    "Cross-validation performance for all candidate settings:"
)

display(
    knn_cross_validation_summary.round(
        {
            "RMSE_SD": 4,
            "MAE_SD": 4,
            "Delta_RMSE_percent_vs_uniform_k8": 2,
        }
    )
)

# Export the complete table
results_directory = RESULTS_DIRECTORY
results_directory.mkdir(
    parents=True,
    exist_ok=True,
)

output_path = (
    results_directory
    / "knn_cross_validation_rmse_mae_delta.csv"
)

knn_cross_validation_summary.to_csv(
    output_path,
    index=False,
)

print("\nTable exported to:")
print(output_path)

# Step 20 — Export imputation-sensitivity scenarios for R analysis

SCENARIO_K_VALUES = [
    SELECTED_K - 1,
    SELECTED_K,
    SELECTED_K + 1
]

scenario_export_records = []

for k_value in SCENARIO_K_VALUES:

    scenario_matrix = partition_data_raw.copy()

    for missing_time, target_variable in sorted(
        EXPECTED_MISSING_POSITIONS,
        key=lambda item: (item[0], item[1])
    ):

        prediction = profile_knn_predict(
            data=partition_data_raw,
            target_variable=target_variable,
            missing_time=missing_time,
            k=k_value,
            weighting=SELECTED_WEIGHTING,
            minimum_common_time_points=MINIMUM_COMMON_TIME_POINTS
        )

        scenario_matrix.loc[
            missing_time,
            target_variable
        ] = prediction["Predicted_value"]

    if scenario_matrix.isna().any().any():
        raise ValueError(
            f"The scenario matrix for k_NN={k_value} "
            "still contains missing values."
        )

    scenario_name = (
        f"{SELECTED_WEIGHTING}_k{k_value}"
    )

    output_file_name = (
        f"partition_coefficients_{scenario_name}.csv"
    )

    output_path = (
        RESULTS_DIRECTORY / output_file_name
    )

    scenario_matrix.reset_index().to_csv(
        output_path,
        index=False
    )

    scenario_export_records.append(
        {
            "Scenario": scenario_name,
            "Weighting": SELECTED_WEIGHTING,
            "k_NN": k_value,
            "Number_of_time_points": (
                scenario_matrix.shape[0]
            ),
            "Number_of_variables": (
                scenario_matrix.shape[1]
            ),
            "Output_file": output_file_name
        }
    )

    print(f"Exported: {output_path}")


# Export a complete-case variable scenario by excluding
# all four variables affected by missing values

excluded_variables_matrix = (
    partition_data_raw
    .drop(columns=AFFECTED_VARIABLES)
    .copy()
)

if excluded_variables_matrix.isna().any().any():
    raise ValueError(
        "The excluded-variable scenario still contains missing values."
    )

excluded_output_file_name = (
    "partition_coefficients_affected_variables_excluded.csv"
)

excluded_output_path = (
    RESULTS_DIRECTORY
    / excluded_output_file_name
)

excluded_variables_matrix.reset_index().to_csv(
    excluded_output_path,
    index=False
)

scenario_export_records.append(
    {
        "Scenario": "affected_variables_excluded",
        "Weighting": "Not applicable",
        "k_NN": np.nan,
        "Number_of_time_points": (
            excluded_variables_matrix.shape[0]
        ),
        "Number_of_variables": (
            excluded_variables_matrix.shape[1]
        ),
        "Output_file": excluded_output_file_name
    }
)

scenario_manifest = pd.DataFrame(
    scenario_export_records
)

scenario_manifest.to_csv(
    RESULTS_DIRECTORY
    / "imputation_sensitivity_scenarios.csv",
    index=False
)

print("\nExported sensitivity-analysis scenarios:")
display(scenario_manifest)
