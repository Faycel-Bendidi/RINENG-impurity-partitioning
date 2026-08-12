source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("FactoMineR", "cluster", "factoextra", "ggplot2"))

# Step 1 — Define paths, packages, and analysis settings
input_path <- file.path(
  "data",
  "processed",
  "partition_coefficients_imputed.csv"
)

tables_directory <- file.path(
  "results",
  "tables"
)

figures_directory <- file.path(
  "results",
  "figures"
)


processed_data_directory <- file.path(
  "data",
  "processed"
)

dir.create(
  processed_data_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Working directory:", getwd(), "\n")
cat("Input file:", input_path, "\n")
cat("Input file found:", file.exists(input_path), "\n")

# Primary imputation setting selected by cross-validation

primary_weighting_method <- "uniform"
primary_k_nn <- 8L
primary_scenario_name <- "uniform_k8"

# Number of HCA groups, denoted G, will be selected later

selected_number_of_groups <- NA_integer_

cat(
  "Primary imputation setting:",
  primary_weighting_method,
  "with k_NN =",
  primary_k_nn,
  "\n"
)

cat(
  "Primary sensitivity scenario:",
  primary_scenario_name,
  "\n"
)

# Step 2 — Import and validate the completed partition-coefficient matrix

data_with_time <- read.csv(
  input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expected_time_points <- c(
  0, 12, 24, 36, 48, 60, 72
)

expected_variables <- c(
  "DP2O5",
  "DAl2O3",
  "DCaO",
  "DFe2O3",
  "DK2O",
  "DMgO",
  "DNa2O",
  "DSO3",
  "DSiO2",
  "DCd",
  "DCr",
  "DCu",
  "DMn",
  "DMo",
  "DNi",
  "DZn",
  "DF"
)

expected_columns <- c(
  "Time_h",
  expected_variables
)

if (!identical(
  names(data_with_time),
  expected_columns
)) {
  stop(
    paste(
      "Unexpected column structure.\nFound:",
      paste(
        names(data_with_time),
        collapse = ", "
      )
    )
  )
}

if (!identical(
  as.numeric(data_with_time$Time_h),
  expected_time_points
)) {
  stop("Unexpected sampling-time sequence.")
}

partition_data <- data_with_time[
  expected_variables
]

partition_data[] <- lapply(
  partition_data,
  as.numeric
)

rownames(partition_data) <- paste0(
  "t",
  data_with_time$Time_h
)

if (!identical(
  dim(partition_data),
  c(7L, 17L)
)) {
  stop(
    paste(
      "Unexpected matrix dimensions:",
      paste(
        dim(partition_data),
        collapse = " x "
      )
    )
  )
}

if (anyNA(partition_data)) {
  stop(
    "The completed analysis matrix still contains missing values."
  )
}

cat(
  "Matrix dimensions:",
  nrow(partition_data),
  "time points x",
  ncol(partition_data),
  "variables\n"
)

cat(
  "Number of missing values:",
  sum(is.na(partition_data)),
  "\n"
)

print(
  round(partition_data, 6)
)

# Step 3 — Center and scale each variable across sampling times

scaled_data <- scale(
  partition_data,
  center = TRUE,
  scale = TRUE
)

# Rename standardized variables from D_i to Z_i

colnames(scaled_data) <- sub(
  "^D",
  "Z",
  colnames(scaled_data)
)

scaled_means <- colMeans(
  scaled_data
)

scaled_standard_deviations <- apply(
  scaled_data,
  2,
  sd
)

scaling_diagnostics <- data.frame(
  D_variable = colnames(partition_data),
  Z_variable = colnames(scaled_data),
  Original_mean = colMeans(partition_data),
  Original_standard_deviation = apply(
    partition_data,
    2,
    sd
  ),
  Scaled_mean = scaled_means,
  Scaled_standard_deviation = scaled_standard_deviations,
  row.names = NULL
)

# Verify that every scaled variable has mean 0 and standard deviation 1
if (!all(
  abs(scaled_means) < 1e-10
)) {
  stop(
    "At least one scaled variable does not have a mean of zero."
  )
}

if (!all(
  abs(
    scaled_standard_deviations - 1
  ) < 1e-10
)) {
  stop(
    paste(
      "At least one scaled variable does not have",
      "a standard deviation of one."
    )
  )
}

# Export scaling diagnostics
write.csv(
  scaling_diagnostics,
  file.path(
    tables_directory,
    "scaling_diagnostics.csv"
  ),
  row.names = FALSE
)

# Export the autoscaled matrix used by PCA and HCA
scaled_data_export <- data.frame(
  Time = rownames(scaled_data),
  as.data.frame(scaled_data),
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  scaled_data_export,
  file.path(
    processed_data_directory,
    "partition_coefficients_standardized_Z.csv"
  ),
  row.names = FALSE
)

cat(
  "Centering and scaling completed successfully.\n"
)

cat(
  "Maximum absolute scaled mean:",
  max(abs(scaled_means)),
  "\n"
)

cat(
  "Maximum absolute deviation from unit standard deviation:",
  max(
    abs(
      scaled_standard_deviations - 1
    )
  ),
  "\n"
)
# Create a display copy and round only numeric columns

scaling_diagnostics_display <- scaling_diagnostics

numeric_columns <- vapply(
  scaling_diagnostics_display,
  is.numeric,
  logical(1)
)

scaling_diagnostics_display[
  numeric_columns
] <- lapply(
  scaling_diagnostics_display[
    numeric_columns
  ],
  round,
  digits = 6
)

print(
  scaling_diagnostics_display,
  row.names = FALSE
)

# Step 4 — Perform PCA on the centered and scaled matrix

# The data have already been centered and scaled in Step 3.
# Therefore, scale.unit is set to FALSE to avoid a second scaling.

pca_result <- FactoMineR::PCA(
  as.data.frame(scaled_data),
  scale.unit = FALSE,
  ncp = 6,
  graph = FALSE
)

pca_eigenvalues <- data.frame(
  Dimension = rownames(pca_result$eig),
  Eigenvalue = pca_result$eig[, 1],
  Variance_percent = pca_result$eig[, 2],
  Cumulative_variance_percent = pca_result$eig[, 3],
  row.names = NULL,
  check.names = FALSE
)

# Export the PCA eigenvalues and explained variance
write.csv(
  pca_eigenvalues,
  file.path(
    tables_directory,
    "pca_eigenvalues.csv"
  ),
  row.names = FALSE
)

# Prepare a rounded copy for display
pca_eigenvalues_display <- pca_eigenvalues

numeric_columns <- vapply(
  pca_eigenvalues_display,
  is.numeric,
  logical(1)
)

pca_eigenvalues_display[
  numeric_columns
] <- lapply(
  pca_eigenvalues_display[
    numeric_columns
  ],
  round,
  digits = 4
)

cat(
  "PCA completed successfully.\n"
)

cat(
  "Variance explained by PC1 and PC2:",
  round(
    sum(
      pca_result$eig[1:2, 2]
    ),
    4
  ),
  "%\n"
)

print(
  pca_eigenvalues_display,
  row.names = FALSE
)

# Step 5 — Save PCA coordinates, contributions, and cos2 values

# Individual coordinates: sampling times
pca_individual_coordinates <- data.frame(
  Time = rownames(
    pca_result$ind$coord
  ),
  pca_result$ind$coord,
  row.names = NULL,
  check.names = FALSE
)

# Variable coordinates on the principal components
pca_variable_coordinates <- data.frame(
  Variable = rownames(
    pca_result$var$coord
  ),
  pca_result$var$coord,
  row.names = NULL,
  check.names = FALSE
)

# Variable contributions to each principal component
pca_variable_contributions <- data.frame(
  Variable = rownames(
    pca_result$var$contrib
  ),
  pca_result$var$contrib,
  row.names = NULL,
  check.names = FALSE
)

# Quality of representation of each variable
pca_variable_cos2 <- data.frame(
  Variable = rownames(
    pca_result$var$cos2
  ),
  pca_result$var$cos2,
  row.names = NULL,
  check.names = FALSE
)

# Export the result tables
write.csv(
  pca_individual_coordinates,
  file.path(
    tables_directory,
    "pca_individual_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_variable_coordinates,
  file.path(
    tables_directory,
    "pca_variable_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_variable_contributions,
  file.path(
    tables_directory,
    "pca_variable_contributions.csv"
  ),
  row.names = FALSE
)

write.csv(
  pca_variable_cos2,
  file.path(
    tables_directory,
    "pca_variable_cos2.csv"
  ),
  row.names = FALSE
)

# Prepare a compact table for interpretation of PC1, PC2, and PC3
pca_variable_summary <- data.frame(
  Variable = rownames(
    pca_result$var$coord
  ),
  PC1_coordinate = pca_result$var$coord[, 1],
  PC2_coordinate = pca_result$var$coord[, 2],
  PC3_coordinate = pca_result$var$coord[, 3],
  PC1_contribution = pca_result$var$contrib[, 1],
  PC2_contribution = pca_result$var$contrib[, 2],
  PC3_contribution = pca_result$var$contrib[, 3],
  PC1_PC2_cos2 = (
    pca_result$var$cos2[, 1]
    + pca_result$var$cos2[, 2]
  ),
  PC1_PC2_PC3_cos2 = (
    pca_result$var$cos2[, 1]
    + pca_result$var$cos2[, 2]
    + pca_result$var$cos2[, 3]
  ),
  row.names = NULL
)

write.csv(
  pca_variable_summary,
  file.path(
    tables_directory,
    "pca_variable_summary_pc1_pc2_pc3.csv"
  ),
  row.names = FALSE
)

# Round only numeric columns for display
pca_variable_summary_display <- pca_variable_summary

numeric_columns <- vapply(
  pca_variable_summary_display,
  is.numeric,
  logical(1)
)

pca_variable_summary_display[
  numeric_columns
] <- lapply(
  pca_variable_summary_display[
    numeric_columns
  ],
  round,
  digits = 4
)

cat(
  "PCA coordinates and diagnostic tables exported successfully.\n"
)

print(
  pca_variable_summary_display,
  row.names = FALSE
)
# Step 6 — Import and validate imputation-sensitivity scenarios

scenario_manifest_path <- file.path(
  processed_data_directory,
  "imputation_sensitivity_scenarios.csv"
)

if (!file.exists(scenario_manifest_path)) {
  stop(
    paste(
      "Scenario manifest not found:",
      scenario_manifest_path
    )
  )
}

scenario_manifest <- read.csv(
  scenario_manifest_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_manifest_columns <- c(
  "Scenario",
  "Weighting",
  "k_NN",
  "Number_of_time_points",
  "Number_of_variables",
  "Output_file"
)

if (!identical(
  names(scenario_manifest),
  required_manifest_columns
)) {
  stop(
    paste(
      "Unexpected scenario-manifest structure.\nFound:",
      paste(names(scenario_manifest), collapse = ", ")
    )
  )
}

affected_variables <- c(
  "DK2O",
  "DCu",
  "DMo",
  "DNi"
)

expected_excluded_variables <- setdiff(
  expected_variables,
  affected_variables
)

imputation_scenarios <- list()
scenario_validation_records <- list()

for (row_index in seq_len(nrow(scenario_manifest))) {
  
  scenario_name <- scenario_manifest$Scenario[row_index]
  
  scenario_path <- file.path(
    processed_data_directory,
    scenario_manifest$Output_file[row_index]
  )
  
  if (!file.exists(scenario_path)) {
    stop(
      paste(
        "Scenario file not found:",
        scenario_path
      )
    )
  }
  
  scenario_data_with_time <- read.csv(
    scenario_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!isTRUE(
    all.equal(
      as.numeric(scenario_data_with_time$Time_h),
      expected_time_points
    )
  )) {
    stop(
      paste(
        "Unexpected time-point sequence in scenario:",
        scenario_name
      )
    )
  }
  
  if (scenario_name == "affected_variables_excluded") {
    
    scenario_expected_variables <- (
      expected_excluded_variables
    )
    
  } else {
    
    scenario_expected_variables <- expected_variables
  }
  
  expected_scenario_columns <- c(
    "Time_h",
    scenario_expected_variables
  )
  
  if (!identical(
    names(scenario_data_with_time),
    expected_scenario_columns
  )) {
    stop(
      paste(
        "Unexpected variables in scenario:",
        scenario_name,
        "\nFound:",
        paste(
          names(scenario_data_with_time),
          collapse = ", "
        )
      )
    )
  }
  
  scenario_data <- scenario_data_with_time[
    scenario_expected_variables
  ]
  
  scenario_data[] <- lapply(
    scenario_data,
    as.numeric
  )
  
  rownames(scenario_data) <- paste0(
    "t",
    scenario_data_with_time$Time_h
  )
  
  if (anyNA(scenario_data)) {
    stop(
      paste(
        "Missing values remain in scenario:",
        scenario_name
      )
    )
  }
  
  expected_number_of_variables <- as.integer(
    scenario_manifest$Number_of_variables[row_index]
  )
  
  if (
    nrow(scenario_data) != 7L ||
    ncol(scenario_data) != expected_number_of_variables
  ) {
    stop(
      paste(
        "Unexpected dimensions for scenario:",
        scenario_name
      )
    )
  }
  
  imputation_scenarios[[scenario_name]] <- scenario_data
  
  scenario_validation_records[[row_index]] <- data.frame(
    Scenario = scenario_name,
    Number_of_time_points = nrow(scenario_data),
    Number_of_variables = ncol(scenario_data),
    Number_of_missing_values = sum(
      is.na(scenario_data)
    ),
    File_found = file.exists(scenario_path),
    stringsAsFactors = FALSE
  )
}

scenario_validation_summary <- do.call(
  rbind,
  scenario_validation_records
)

# Confirm that the main k_NN = 8 scenario matches the matrix
# previously analyzed in Steps 2–5.

main_scenario_comparison <- all.equal(
  imputation_scenarios[["uniform_k8"]],
  partition_data,
  tolerance = 1e-12,
  check.attributes = TRUE
)

if (!isTRUE(main_scenario_comparison)) {
  stop(
    paste(
      "The uniform_k8 scenario does not match",
      "the reference matrix from Steps 2–5:",
      main_scenario_comparison
    )
  )
}

write.csv(
  scenario_validation_summary,
  file.path(
    tables_directory,
    "imputation_scenario_validation.csv"
  ),
  row.names = FALSE
)

cat(
  "All imputation-sensitivity scenarios were imported",
  "and validated successfully.\n"
)

cat(
  "The uniform_k8 scenario matches the reference matrix:",
  isTRUE(main_scenario_comparison),
  "\n"
)

print(
  scenario_validation_summary,
  row.names = FALSE
)
# Step 6B — Import the Python cross-validation summary

cross_validation_summary_path <- file.path(
  processed_data_directory,
  "knn_cross_validation_rmse_mae_delta.csv"
)

if (!file.exists(cross_validation_summary_path)) {
  stop(
    "The k-NN cross-validation summary file was not found."
  )
}

knn_cross_validation_summary <- read.csv(
  cross_validation_summary_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "The Python k-NN cross-validation summary was imported successfully.\n"
)

print(
  knn_cross_validation_summary,
  row.names = FALSE
)

# Step 7 — Run PCA for all imputation-sensitivity scenarios

if (!exists("imputation_scenarios")) {
  stop(
    "The object 'imputation_scenarios' was not found. "
  )
}

if (!primary_scenario_name %in% names(imputation_scenarios)) {
  stop(
    paste(
      "The primary scenario was not found:",
      primary_scenario_name
    )
  )
}

pca_sensitivity_results <- list()
pca_variance_records <- list()

scenario_index <- 1L

for (scenario_name in names(imputation_scenarios)) {
  
  scenario_data <- imputation_scenarios[[scenario_name]]
  
  # Convert all variables to numeric values
  numeric_scenario_data <- as.data.frame(
    scenario_data,
    check.names = FALSE
  )
  
  numeric_scenario_data[] <- lapply(
    numeric_scenario_data,
    as.numeric
  )
  
  rownames(numeric_scenario_data) <- rownames(
    scenario_data
  )
  
  if (anyNA(numeric_scenario_data)) {
    stop(
      paste(
        "Missing values were detected in scenario:",
        scenario_name
      )
    )
  }
  
  variable_standard_deviations <- apply(
    numeric_scenario_data,
    2,
    sd
  )
  
  if (
    any(!is.finite(variable_standard_deviations)) ||
    any(variable_standard_deviations == 0)
  ) {
    stop(
      paste(
        "A constant or invalid variable was detected in scenario:",
        scenario_name
      )
    )
  }
  
  # Center and scale each variable across sampling times
  scaled_scenario_data <- scale(
    numeric_scenario_data,
    center = TRUE,
    scale = TRUE
  )
  
  colnames(scaled_scenario_data) <- sub(
    "^D",
    "Z",
    colnames(scaled_scenario_data)
  )
  
  maximum_number_of_components <- min(
    nrow(scaled_scenario_data) - 1L,
    ncol(scaled_scenario_data)
  )
  
  # Run PCA on the standardized matrix
  scenario_pca <- FactoMineR::PCA(
    as.data.frame(scaled_scenario_data),
    scale.unit = FALSE,
    ncp = maximum_number_of_components,
    graph = FALSE
  )
  
  eigenvalue_table <- as.data.frame(
    scenario_pca$eig,
    check.names = FALSE
  )
  
  eigenvalue_table$Principal_component <- rownames(
    eigenvalue_table
  )
  
  rownames(eigenvalue_table) <- NULL
  
  eigenvalue_table$Scenario <- scenario_name
  
  eigenvalue_table <- eigenvalue_table[
    ,
    c(
      "Scenario",
      "Principal_component",
      setdiff(
        names(eigenvalue_table),
        c(
          "Scenario",
          "Principal_component"
        )
      )
    )
  ]
  
  individual_coordinates <- as.data.frame(
    scenario_pca$ind$coord,
    check.names = FALSE
  )
  
  individual_coordinates$Time <- rownames(
    individual_coordinates
  )
  
  individual_coordinates$Scenario <- scenario_name
  
  rownames(individual_coordinates) <- NULL
  
  variable_coordinates <- as.data.frame(
    scenario_pca$var$coord,
    check.names = FALSE
  )
  
  variable_coordinates$Variable <- rownames(
    variable_coordinates
  )
  
  variable_coordinates$Scenario <- scenario_name
  
  rownames(variable_coordinates) <- NULL
  
  variable_contributions <- as.data.frame(
    scenario_pca$var$contrib,
    check.names = FALSE
  )
  
  variable_contributions$Variable <- rownames(
    variable_contributions
  )
  
  variable_contributions$Scenario <- scenario_name
  
  rownames(variable_contributions) <- NULL
  
  variable_cos2 <- as.data.frame(
    scenario_pca$var$cos2,
    check.names = FALSE
  )
  
  variable_cos2$Variable <- rownames(
    variable_cos2
  )
  
  variable_cos2$Scenario <- scenario_name
  
  rownames(variable_cos2) <- NULL
  
  pca_sensitivity_results[[scenario_name]] <- list(
    original_data = numeric_scenario_data,
    scaled_data = scaled_scenario_data,
    pca_result = scenario_pca,
    eigenvalues = eigenvalue_table,
    individual_coordinates = individual_coordinates,
    variable_coordinates = variable_coordinates,
    variable_contributions = variable_contributions,
    variable_cos2 = variable_cos2
  )
  
  pca_variance_records[[scenario_index]] <- data.frame(
    Scenario = scenario_name,
    Number_of_variables = ncol(
      numeric_scenario_data
    ),
    PC1_percent = scenario_pca$eig[1, 2],
    PC2_percent = scenario_pca$eig[2, 2],
    PC3_percent = scenario_pca$eig[3, 2],
    PC1_PC2_cumulative_percent = scenario_pca$eig[2, 3],
    PC1_PC2_PC3_cumulative_percent = scenario_pca$eig[3, 3],
    stringsAsFactors = FALSE
  )
  
  write.csv(
    eigenvalue_table,
    file.path(
      tables_directory,
      paste0(
        "pca_eigenvalues_",
        scenario_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    individual_coordinates,
    file.path(
      tables_directory,
      paste0(
        "pca_individual_coordinates_",
        scenario_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    variable_coordinates,
    file.path(
      tables_directory,
      paste0(
        "pca_variable_coordinates_",
        scenario_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  scenario_index <- scenario_index + 1L
}

pca_variance_summary <- do.call(
  rbind,
  pca_variance_records
)

rownames(pca_variance_summary) <- NULL

write.csv(
  pca_variance_summary,
  file.path(
    tables_directory,
    "pca_variance_summary_all_scenarios.csv"
  ),
  row.names = FALSE
)

cat(
  "PCA was completed successfully for all sensitivity scenarios.\n"
)

print(
  pca_variance_summary,
  row.names = FALSE
)
# Verify consistency with the primary PCA from Steps 4-5

if (exists("pca_result")) {
  
  primary_pca_consistency <- all.equal(
    pca_sensitivity_results[[
      primary_scenario_name
    ]]$pca_result$eig,
    pca_result$eig,
    tolerance = 1e-10,
    check.attributes = TRUE
  )
  
  cat(
    "The primary-scenario PCA matches the PCA from Steps 4-5:",
    isTRUE(primary_pca_consistency),
    "\n"
  )
  
  if (!isTRUE(primary_pca_consistency)) {
    print(primary_pca_consistency)
  }
}

# Step 8 — Quantify PCA robustness across scenarios

if (!exists("pca_sensitivity_results")) {
  stop(
    "The object 'pca_sensitivity_results' was not found. Run Step 7 first."
  )
}

if (!exists("pca_variance_summary")) {
  stop(
    "The object 'pca_variance_summary' was not found. Run Step 7 first."
  )
}

if (!exists("primary_scenario_name")) {
  stop(
    "The object 'primary_scenario_name' was not found. Run Step 1 first."
  )
}

reference_scenario_name <- primary_scenario_name

if (!reference_scenario_name %in% names(pca_sensitivity_results)) {
  stop(
    paste(
      "The reference PCA scenario was not found:",
      reference_scenario_name
    )
  )
}


# Extract PCA coordinates and assign the labels as row names

extract_coordinate_matrix <- function(
    coordinate_table,
    label_column,
    number_of_dimensions
) {
  
  dimension_columns <- paste0(
    "Dim.",
    seq_len(number_of_dimensions)
  )
  
  required_columns <- c(
    label_column,
    dimension_columns
  )
  
  if (!all(required_columns %in% names(coordinate_table))) {
    stop(
      paste(
        "Required PCA columns were not found:",
        paste(required_columns, collapse = ", ")
      )
    )
  }
  
  coordinate_matrix <- as.matrix(
    coordinate_table[
      ,
      dimension_columns,
      drop = FALSE
    ]
  )
  
  storage.mode(coordinate_matrix) <- "double"
  
  rownames(coordinate_matrix) <- coordinate_table[
    ,
    label_column,
    drop = TRUE
  ]
  
  return(coordinate_matrix)
}


# Compare two PCA configurations after matching common labels

compare_pca_configurations <- function(
    reference_matrix,
    comparison_matrix
) {
  
  common_labels <- intersect(
    rownames(reference_matrix),
    rownames(comparison_matrix)
  )
  
  if (length(common_labels) < 3L) {
    stop(
      "At least three common labels are required."
    )
  }
  
  reference_common <- reference_matrix[
    common_labels,
    ,
    drop = FALSE
  ]
  
  comparison_common <- comparison_matrix[
    common_labels,
    ,
    drop = FALSE
  ]
  
  # Compare the pairwise-distance structure
  
  reference_distances <- as.vector(
    dist(reference_common)
  )
  
  comparison_distances <- as.vector(
    dist(comparison_common)
  )
  
  distance_correlation <- cor(
    reference_distances,
    comparison_distances,
    method = "pearson"
  )
  
  # Center both PCA configurations
  
  reference_centered <- scale(
    reference_common,
    center = TRUE,
    scale = FALSE
  )
  
  comparison_centered <- scale(
    comparison_common,
    center = TRUE,
    scale = FALSE
  )
  
  reference_norm <- sqrt(
    sum(reference_centered^2)
  )
  
  comparison_norm <- sqrt(
    sum(comparison_centered^2)
  )
  
  if (
    !is.finite(reference_norm) ||
    !is.finite(comparison_norm) ||
    reference_norm == 0 ||
    comparison_norm == 0
  ) {
    stop(
      "A PCA configuration has zero or invalid total variation."
    )
  }
  
  reference_normalized <- (
    reference_centered /
      reference_norm
  )
  
  comparison_normalized <- (
    comparison_centered /
      comparison_norm
  )
  
  # Orthogonal Procrustes alignment
  
  svd_result <- svd(
    crossprod(
      comparison_normalized,
      reference_normalized
    )
  )
  
  rotation_matrix <- (
    svd_result$u %*%
      t(svd_result$v)
  )
  
  comparison_aligned <- (
    comparison_normalized %*%
      rotation_matrix
  )
  
  procrustes_similarity <- sum(
    svd_result$d
  )
  
  procrustes_similarity <- min(
    1,
    max(
      0,
      procrustes_similarity
    )
  )
  
  normalized_residual_rmse <- sqrt(
    mean(
      (
        reference_normalized -
          comparison_aligned
      )^2
    )
  )
  
  return(
    list(
      number_of_common_labels = length(common_labels),
      distance_correlation = distance_correlation,
      procrustes_similarity = procrustes_similarity,
      normalized_residual_rmse = normalized_residual_rmse
    )
  )
}


# Calculate explained-variance differences relative to uniform_k8

variance_columns <- c(
  "PC1_percent",
  "PC2_percent",
  "PC3_percent",
  "PC1_PC2_cumulative_percent",
  "PC1_PC2_PC3_cumulative_percent"
)

reference_variance_row <- pca_variance_summary[
  pca_variance_summary$Scenario == reference_scenario_name,
  ,
  drop = FALSE
]

if (nrow(reference_variance_row) != 1L) {
  stop(
    "The reference variance row could not be identified uniquely."
  )
}

pca_variance_differences <- pca_variance_summary

for (variance_column in variance_columns) {
  
  difference_column_name <- paste0(
    "Delta_",
    variance_column,
    "_vs_",
    reference_scenario_name
  )
  
  pca_variance_differences[
    ,
    difference_column_name
  ] <- (
    pca_variance_differences[
      ,
      variance_column
    ] -
      reference_variance_row[
        1,
        variance_column
      ]
  )
}

write.csv(
  pca_variance_differences,
  file.path(
    tables_directory,
    "pca_variance_differences_vs_uniform_k8.csv"
  ),
  row.names = FALSE
)

cat(
  "Explained-variance differences were calculated successfully.\n"
)

variance_display <- pca_variance_differences

numeric_variance_columns <- vapply(
  variance_display,
  is.numeric,
  logical(1)
)

variance_display[
  numeric_variance_columns
] <- lapply(
  variance_display[
    numeric_variance_columns
  ],
  round,
  digits = 4
)

print(
  variance_display,
  row.names = FALSE
)


# Extract the reference PCA configurations

reference_results <- pca_sensitivity_results[[reference_scenario_name]]

reference_time_coordinates_2d <- extract_coordinate_matrix(
  reference_results$individual_coordinates,
  "Time",
  2L
)

reference_time_coordinates_3d <- extract_coordinate_matrix(
  reference_results$individual_coordinates,
  "Time",
  3L
)

reference_variable_coordinates_2d <- extract_coordinate_matrix(
  reference_results$variable_coordinates,
  "Variable",
  2L
)

reference_variable_coordinates_3d <- extract_coordinate_matrix(
  reference_results$variable_coordinates,
  "Variable",
  3L
)


# Compare all PCA configurations with uniform_k8

configuration_comparison_records <- list()

scenario_index <- 1L

for (scenario_name in names(pca_sensitivity_results)) {
  
  scenario_results <- pca_sensitivity_results[[scenario_name]]
  
  scenario_time_coordinates_2d <- extract_coordinate_matrix(
    scenario_results$individual_coordinates,
    "Time",
    2L
  )
  
  scenario_time_coordinates_3d <- extract_coordinate_matrix(
    scenario_results$individual_coordinates,
    "Time",
    3L
  )
  
  scenario_variable_coordinates_2d <- extract_coordinate_matrix(
    scenario_results$variable_coordinates,
    "Variable",
    2L
  )
  
  scenario_variable_coordinates_3d <- extract_coordinate_matrix(
    scenario_results$variable_coordinates,
    "Variable",
    3L
  )
  
  time_comparison_2d <- compare_pca_configurations(
    reference_time_coordinates_2d,
    scenario_time_coordinates_2d
  )
  
  time_comparison_3d <- compare_pca_configurations(
    reference_time_coordinates_3d,
    scenario_time_coordinates_3d
  )
  
  variable_comparison_2d <- compare_pca_configurations(
    reference_variable_coordinates_2d,
    scenario_variable_coordinates_2d
  )
  
  variable_comparison_3d <- compare_pca_configurations(
    reference_variable_coordinates_3d,
    scenario_variable_coordinates_3d
  )
  
  configuration_comparison_records[[scenario_index]] <- data.frame(
    Scenario = scenario_name,
    
    Number_of_common_times = (
      time_comparison_3d$number_of_common_labels
    ),
    
    Time_distance_correlation_2D = (
      time_comparison_2d$distance_correlation
    ),
    
    Time_Procrustes_similarity_2D = (
      time_comparison_2d$procrustes_similarity
    ),
    
    Time_residual_RMSE_2D = (
      time_comparison_2d$normalized_residual_rmse
    ),
    
    Time_distance_correlation_3D = (
      time_comparison_3d$distance_correlation
    ),
    
    Time_Procrustes_similarity_3D = (
      time_comparison_3d$procrustes_similarity
    ),
    
    Time_residual_RMSE_3D = (
      time_comparison_3d$normalized_residual_rmse
    ),
    
    Number_of_common_variables = (
      variable_comparison_3d$number_of_common_labels
    ),
    
    Variable_distance_correlation_2D = (
      variable_comparison_2d$distance_correlation
    ),
    
    Variable_Procrustes_similarity_2D = (
      variable_comparison_2d$procrustes_similarity
    ),
    
    Variable_residual_RMSE_2D = (
      variable_comparison_2d$normalized_residual_rmse
    ),
    
    Variable_distance_correlation_3D = (
      variable_comparison_3d$distance_correlation
    ),
    
    Variable_Procrustes_similarity_3D = (
      variable_comparison_3d$procrustes_similarity
    ),
    
    Variable_residual_RMSE_3D = (
      variable_comparison_3d$normalized_residual_rmse
    ),
    
    stringsAsFactors = FALSE
  )
  
  scenario_index <- scenario_index + 1L
}

pca_configuration_comparison <- do.call(
  rbind,
  configuration_comparison_records
)

rownames(pca_configuration_comparison) <- NULL

write.csv(
  pca_configuration_comparison,
  file.path(
    tables_directory,
    "pca_configuration_similarity_vs_uniform_k8.csv"
  ),
  row.names = FALSE
)

cat(
  "PCA configuration comparisons were completed successfully.\n"
)

configuration_display <- pca_configuration_comparison

numeric_configuration_columns <- vapply(
  configuration_display,
  is.numeric,
  logical(1)
)

configuration_display[numeric_configuration_columns] <- lapply(
  configuration_display[numeric_configuration_columns],
  round,
  digits = 4
)

print(
  configuration_display,
  row.names = FALSE
)

# Step 9 — Evaluate the number of HCA groups G
# using the primary uniform_k8 scenario

if (!exists("pca_sensitivity_results")) {
  stop(
    "The object 'pca_sensitivity_results' was not found. Run Step 7 first."
  )
}

if (!exists("primary_scenario_name")) {
  stop(
    "The object 'primary_scenario_name' was not found. Run Step 1 first."
  )
}

if (!requireNamespace("cluster", quietly = TRUE)) {
  stop(
    "The 'cluster' package is required."
  )
}

if (!primary_scenario_name %in% names(pca_sensitivity_results)) {
  stop(
    paste(
      "The primary scenario was not found:",
      primary_scenario_name
    )
  )
}

# Extract the standardized matrix of the primary scenario

primary_scaled_matrix <- pca_sensitivity_results[[primary_scenario_name]]$scaled_data

# Variables are clustered according to their standardized temporal profiles

primary_variable_profiles <- t(
  primary_scaled_matrix
)

number_of_variables <- nrow(
  primary_variable_profiles
)

cat(
  "Primary HCA scenario:",
  primary_scenario_name,
  "\n"
)

cat(
  "Number of variables:",
  number_of_variables,
  "\n"
)

# Calculate Euclidean distances between variable profiles

primary_variable_distance <- dist(
  primary_variable_profiles,
  method = "euclidean"
)

# Perform Ward hierarchical clustering

primary_hca_result <- hclust(
  primary_variable_distance,
  method = "ward.D2"
)

# Calculate the cophenetic correlation

primary_cophenetic_distance <- cophenetic(
  primary_hca_result
)

primary_cophenetic_correlation <- cor(
  as.vector(primary_variable_distance),
  as.vector(primary_cophenetic_distance),
  method = "pearson"
)

cat(
  "Cophenetic correlation:",
  round(primary_cophenetic_correlation, 4),
  "\n"
)

# Define candidate numbers of groups G

elbow_group_numbers <- 1:10

silhouette_group_numbers <- 2:10

# Calculate the k-means elbow diagnostic

elbow_records <- list()

for (number_of_groups in elbow_group_numbers) {
  
  set.seed(2026)
  
  kmeans_result <- kmeans(
    primary_variable_profiles,
    centers = number_of_groups,
    nstart = 100,
    iter.max = 100
  )
  
  elbow_records[[as.character(number_of_groups)]] <- data.frame(
    Scenario = primary_scenario_name,
    Number_of_groups_G = number_of_groups,
    Total_within_cluster_SS = kmeans_result$tot.withinss,
    stringsAsFactors = FALSE
  )
}

elbow_summary <- do.call(
  rbind,
  elbow_records
)

rownames(elbow_summary) <- NULL

# Calculate Ward HCA silhouette diagnostics

silhouette_records <- list()

for (number_of_groups in silhouette_group_numbers) {
  
  hca_membership <- cutree(
    primary_hca_result,
    k = number_of_groups
  )
  
  silhouette_result <- cluster::silhouette(
    hca_membership,
    primary_variable_distance
  )
  
  silhouette_widths <- silhouette_result[, "sil_width"]
  
  silhouette_records[[as.character(number_of_groups)]] <- data.frame(
    Scenario = primary_scenario_name,
    Number_of_groups_G = number_of_groups,
    Mean_silhouette_width = mean(silhouette_widths),
    Minimum_silhouette_width = min(silhouette_widths),
    Number_of_negative_silhouettes = sum(silhouette_widths < 0),
    stringsAsFactors = FALSE
  )
}

silhouette_summary <- do.call(
  rbind,
  silhouette_records
)

rownames(silhouette_summary) <- NULL

# Export the diagnostic tables

write.csv(
  elbow_summary,
  file.path(
    tables_directory,
    "hca_elbow_diagnostics_uniform_k8.csv"
  ),
  row.names = FALSE
)

write.csv(
  silhouette_summary,
  file.path(
    tables_directory,
    "hca_silhouette_diagnostics_uniform_k8.csv"
  ),
  row.names = FALSE
)

cophenetic_summary <- data.frame(
  Scenario = primary_scenario_name,
  Number_of_variables = number_of_variables,
  Cophenetic_correlation = primary_cophenetic_correlation,
  stringsAsFactors = FALSE
)

write.csv(
  cophenetic_summary,
  file.path(
    tables_directory,
    "hca_cophenetic_correlation_uniform_k8.csv"
  ),
  row.names = FALSE
)

cat(
  "HCA group-number diagnostics were completed successfully.\n"
)

cat(
  "\nElbow diagnostic:\n"
)

elbow_display <- elbow_summary

elbow_display$Total_within_cluster_SS <- round(
  elbow_display$Total_within_cluster_SS,
  4
)

print(
  elbow_display,
  row.names = FALSE
)

cat(
  "\nHCA silhouette diagnostic:\n"
)

silhouette_display <- silhouette_summary

silhouette_display$Mean_silhouette_width <- round(
  silhouette_display$Mean_silhouette_width,
  4
)

silhouette_display$Minimum_silhouette_width <- round(
  silhouette_display$Minimum_silhouette_width,
  4
)

print(
  silhouette_display,
  row.names = FALSE
)
# Step 9B — Display elbow and silhouette diagnostic plots

if (!exists("elbow_summary")) {
  stop(
    "The object 'elbow_summary' was not found. Run Step 9 first."
  )
}

if (!exists("silhouette_summary")) {
  stop(
    "The object 'silhouette_summary' was not found. Run Step 9 first."
  )
}

elbow_plot <- ggplot2::ggplot(
  elbow_summary,
  ggplot2::aes(
    x = Number_of_groups_G,
    y = Total_within_cluster_SS
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::scale_x_continuous(
    breaks = elbow_group_numbers
  ) +
  ggplot2::labs(
    title = "Elbow diagnostic for variable profiles",
    x = "Number of groups (G)",
    y = "Total within-cluster sum of squares"
  ) +
  ggplot2::theme_bw()

print(
  elbow_plot
)

silhouette_plot <- ggplot2::ggplot(
  silhouette_summary,
  ggplot2::aes(
    x = Number_of_groups_G,
    y = Mean_silhouette_width
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::scale_x_continuous(
    breaks = silhouette_group_numbers
  ) +
  ggplot2::labs(
    title = "Ward HCA silhouette diagnostic",
    x = "Number of groups (G)",
    y = "Mean silhouette width"
  ) +
  ggplot2::theme_bw()

print(
  silhouette_plot
)

cat(
  "Elbow and silhouette plots were generated successfully.\n"
)

# Step 10 — Examine the Ward dendrogram and individual silhouettes

if (!exists("primary_hca_result")) {
  stop(
    "The object 'primary_hca_result' was not found. Run Step 9 first."
  )
}

if (!exists("primary_variable_distance")) {
  stop(
    "The object 'primary_variable_distance' was not found. Run Step 9 first."
  )
}

if (!requireNamespace("cluster", quietly = TRUE)) {
  stop(
    "The 'cluster' package is required."
  )
}

if (!requireNamespace("factoextra", quietly = TRUE)) {
  stop(
    "The 'factoextra' package is required."
  )
}

if (!exists("figure_font_family")) {
  if (.Platform$OS.type == "windows") {
    windowsFonts(Times = windowsFont("Times New Roman"))
    figure_font_family <- "Times"
  } else {
    figure_font_family <- "serif"
  }
}

candidate_G_values <- 2:4


# Create a table comparing memberships for G = 2, 3, and 4

hca_membership_comparison <- data.frame(
  Variable = primary_hca_result$labels,
  stringsAsFactors = FALSE
)

hca_group_size_records <- list()
individual_silhouette_records <- list()
candidate_dendrogram_plots <- list()
candidate_silhouette_plots <- list()

record_index <- 1L


# Helper function to recover variable names from silhouette results

extract_silhouette_variable_names <- function(
    silhouette_result,
    membership
) {
  
  silhouette_ids <- rownames(
    silhouette_result
  )
  
  if (is.null(silhouette_ids)) {
    return(
      names(membership)
    )
  }
  
  numeric_ids <- suppressWarnings(
    as.integer(silhouette_ids)
  )
  
  numeric_ids_are_valid <- (
    all(!is.na(numeric_ids)) &&
      all(numeric_ids >= 1L) &&
      all(numeric_ids <= length(membership))
  )
  
  if (numeric_ids_are_valid) {
    return(
      names(membership)[numeric_ids]
    )
  }
  
  return(
    silhouette_ids
  )
}


for (G in candidate_G_values) {
  
  hca_membership <- cutree(
    primary_hca_result,
    k = G
  )
  
  membership_column_name <- paste0(
    "Group_G",
    G
  )
  
  hca_membership_comparison[[membership_column_name]] <- as.integer(
    hca_membership[hca_membership_comparison$Variable]
  )
  
  
  # Calculate group sizes
  
  group_sizes <- table(
    hca_membership
  )
  
  hca_group_size_records[[paste0("G_", G)]] <- data.frame(
    Number_of_groups_G = G,
    Group = as.integer(names(group_sizes)),
    Number_of_variables = as.integer(group_sizes),
    stringsAsFactors = FALSE
  )
  
  
  # Calculate individual silhouette widths
  
  silhouette_result <- cluster::silhouette(
    hca_membership,
    primary_variable_distance
  )
  
  silhouette_variable_names <- extract_silhouette_variable_names(
    silhouette_result,
    hca_membership
  )
  
  individual_silhouette_records[[paste0("G_", G)]] <- data.frame(
    Variable = silhouette_variable_names,
    Number_of_groups_G = G,
    Group = as.integer(silhouette_result[, "cluster"]),
    Neighbor_group = as.integer(silhouette_result[, "neighbor"]),
    Silhouette_width = as.numeric(silhouette_result[, "sil_width"]),
    stringsAsFactors = FALSE
  )
  
  
  # Generate a diagnostic dendrogram
  
  candidate_dendrogram_plots[[paste0("G_", G)]] <- factoextra::fviz_dend(
    primary_hca_result,
    k = G,
    rect = TRUE,
    rect_fill = FALSE,
    show_labels = TRUE,
    cex = 0.8,
    main = paste0(
      "Ward.D2 HCA dendrogram — G = ",
      G
    ),
    xlab = "",
    ylab = "Ward.D2 height"
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = figure_font_family
      )
    )
  
  
  # Generate an individual-silhouette plot
  
  candidate_silhouette_plots[[paste0("G_", G)]] <- factoextra::fviz_silhouette(
    silhouette_result
  ) +
    ggplot2::labs(
      title = paste0(
        "Individual silhouette widths — G = ",
        G
      ),
      x = "Variables",
      y = "Silhouette width"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = figure_font_family
      )
    )
  
  record_index <- record_index + 1L
}


# Combine and export group-size results

hca_group_size_summary <- do.call(
  rbind,
  hca_group_size_records
)

rownames(hca_group_size_summary) <- NULL


# Combine and export individual-silhouette results

individual_silhouette_summary <- do.call(
  rbind,
  individual_silhouette_records
)

rownames(individual_silhouette_summary) <- NULL


# Sort membership table by the G = 3 partition

hca_membership_comparison <- hca_membership_comparison[
  order(
    hca_membership_comparison$Group_G3,
    hca_membership_comparison$Variable
  ),
  ,
  drop = FALSE
]

rownames(hca_membership_comparison) <- NULL


# Export diagnostic tables

write.csv(
  hca_membership_comparison,
  file.path(
    tables_directory,
    "hca_membership_comparison_G2_G3_G4.csv"
  ),
  row.names = FALSE
)

write.csv(
  hca_group_size_summary,
  file.path(
    tables_directory,
    "hca_group_sizes_G2_G3_G4.csv"
  ),
  row.names = FALSE
)

write.csv(
  individual_silhouette_summary,
  file.path(
    tables_directory,
    "hca_individual_silhouettes_G2_G3_G4.csv"
  ),
  row.names = FALSE
)


cat(
  "Ward dendrogram and individual-silhouette diagnostics were completed successfully.\n"
)

cat(
  "\nCluster memberships for G = 2, 3, and 4:\n"
)

print(
  hca_membership_comparison,
  row.names = FALSE
)

cat(
  "\nGroup sizes:\n"
)

print(
  hca_group_size_summary,
  row.names = FALSE
)

cat(
  "\nVariables with the smallest silhouette widths for each value of G:\n"
)

smallest_silhouettes <- do.call(
  rbind,
  lapply(
    candidate_G_values,
    function(G) {
      
      scenario_table <- individual_silhouette_summary[
        individual_silhouette_summary$Number_of_groups_G == G,
        ,
        drop = FALSE
      ]
      
      scenario_table <- scenario_table[
        order(scenario_table$Silhouette_width),
        ,
        drop = FALSE
      ]
      
      head(
        scenario_table,
        5L
      )
    }
  )
)

smallest_silhouettes$Silhouette_width <- round(
  smallest_silhouettes$Silhouette_width,
  4
)

print(
  smallest_silhouettes,
  row.names = FALSE
)

# Step 10B — Display the uncut Ward dendrogram

if (!exists("primary_hca_result")) {
  stop(
    "The object 'primary_hca_result' was not found. Run Step 9 first."
  )
}

uncut_dendrogram_plot <- factoextra::fviz_dend(
  primary_hca_result,
  k = 1,
  color_labels_by_k = FALSE,
  rect = FALSE,
  show_labels = TRUE,
  cex = 0.8,
  main = "Ward.D2 hierarchical clustering of variable profiles",
  xlab = "",
  ylab = "Ward.D2 height"
)

print(
  uncut_dendrogram_plot
)

cat(
  "The uncut Ward dendrogram was generated successfully.\n"
)

# Step 10C — Display candidate dendrograms

print(candidate_dendrogram_plots[["G_2"]])
print(candidate_dendrogram_plots[["G_3"]])
print(candidate_dendrogram_plots[["G_4"]])

# Step 10D — Display individual-silhouette plots

print(candidate_silhouette_plots[["G_2"]])
print(candidate_silhouette_plots[["G_3"]])
print(candidate_silhouette_plots[["G_4"]])

# Selected HCA solution: parsimonious compromise among all diagnostics

selected_number_of_groups <- 3L

cat(
  "Selected number of HCA groups G:",
  selected_number_of_groups,
  "\n"
)

# Step 11 — Compare HCA stability across imputation scenarios

if (!exists("pca_sensitivity_results")) {
  stop(
    "The object 'pca_sensitivity_results' was not found. Run Step 7 first."
  )
}

if (!exists("selected_number_of_groups")) {
  stop(
    "The selected number of HCA groups was not found."
  )
}

if (is.na(selected_number_of_groups)) {
  stop(
    "The selected number of HCA groups G has not been defined."
  )
}

if (!requireNamespace("cluster", quietly = TRUE)) {
  stop(
    "The 'cluster' package is required."
  )
}


# Calculate the Adjusted Rand Index without an additional package

adjusted_rand_index <- function(
    reference_membership,
    comparison_membership
) {
  
  contingency_table <- table(
    reference_membership,
    comparison_membership
  )
  
  choose_two <- function(x) {
    x * (x - 1) / 2
  }
  
  sum_cells <- sum(
    choose_two(contingency_table)
  )
  
  row_totals <- rowSums(
    contingency_table
  )
  
  column_totals <- colSums(
    contingency_table
  )
  
  sum_rows <- sum(
    choose_two(row_totals)
  )
  
  sum_columns <- sum(
    choose_two(column_totals)
  )
  
  number_of_items <- sum(
    contingency_table
  )
  
  total_pairs <- choose_two(
    number_of_items
  )
  
  expected_index <- (
    sum_rows * sum_columns / total_pairs
  )
  
  maximum_index <- (
    0.5 * (sum_rows + sum_columns)
  )
  
  denominator <- (
    maximum_index - expected_index
  )
  
  if (denominator == 0) {
    return(NA_real_)
  }
  
  adjusted_index <- (
    sum_cells - expected_index
  ) / denominator
  
  return(adjusted_index)
}


# Generate all permutations of group labels

generate_permutations <- function(values) {
  
  if (length(values) == 1L) {
    return(
      matrix(
        values,
        nrow = 1L
      )
    )
  }
  
  permutation_list <- lapply(
    seq_along(values),
    function(index) {
      
      remaining_values <- values[-index]
      
      remaining_permutations <- generate_permutations(
        remaining_values
      )
      
      cbind(
        values[index],
        remaining_permutations
      )
    }
  )
  
  return(
    do.call(
      rbind,
      permutation_list
    )
  )
}


# Align arbitrary group labels with the reference partition

align_group_labels <- function(
    reference_membership,
    comparison_membership
) {
  
  common_variables <- names(reference_membership)[
    names(reference_membership) %in% names(comparison_membership)
  ]
  
  reference_common <- reference_membership[
    common_variables
  ]
  
  comparison_common <- comparison_membership[
    common_variables
  ]
  
  reference_groups <- sort(
    unique(reference_common)
  )
  
  comparison_groups <- sort(
    unique(comparison_common)
  )
  
  if (length(reference_groups) != length(comparison_groups)) {
    stop(
      "The two partitions do not contain the same number of groups."
    )
  }
  
  candidate_mappings <- generate_permutations(
    reference_groups
  )
  
  best_agreement <- -1L
  best_mapped_membership <- NULL
  best_mapping <- NULL
  
  for (row_index in seq_len(nrow(candidate_mappings))) {
    
    candidate_labels <- candidate_mappings[
      row_index,
      ,
      drop = TRUE
    ]
    
    mapped_membership <- candidate_labels[
      match(
        comparison_common,
        comparison_groups
      )
    ]
    
    names(mapped_membership) <- common_variables
    
    agreement_count <- sum(
      mapped_membership == reference_common
    )
    
    if (agreement_count > best_agreement) {
      
      best_agreement <- agreement_count
      best_mapped_membership <- mapped_membership
      
      best_mapping <- data.frame(
        Original_group = comparison_groups,
        Aligned_group = candidate_labels,
        stringsAsFactors = FALSE
      )
    }
  }
  
  changed_variables <- common_variables[
    best_mapped_membership != reference_common
  ]
  
  return(
    list(
      common_variables = common_variables,
      reference_membership = reference_common,
      aligned_membership = best_mapped_membership,
      agreement_percent = (
        100 * best_agreement / length(common_variables)
      ),
      changed_variables = changed_variables,
      mapping = best_mapping
    )
  )
}


# Perform HCA for every sensitivity scenario

hca_sensitivity_results <- list()
hca_diagnostic_records <- list()
hca_membership_records <- list()
hca_group_size_records <- list()

scenario_index <- 1L

for (scenario_name in names(pca_sensitivity_results)) {
  
  scaled_scenario_data <- pca_sensitivity_results[[scenario_name]]$scaled_data
  
  variable_profiles <- t(
    scaled_scenario_data
  )
  
  variable_distance <- dist(
    variable_profiles,
    method = "euclidean"
  )
  
  hca_result <- hclust(
    variable_distance,
    method = "ward.D2"
  )
  
  group_membership <- cutree(
    hca_result,
    k = selected_number_of_groups
  )
  
  silhouette_result <- cluster::silhouette(
    group_membership,
    variable_distance
  )
  
  silhouette_widths <- silhouette_result[
    ,
    "sil_width"
  ]
  
  cophenetic_distance <- cophenetic(
    hca_result
  )
  
  cophenetic_correlation <- cor(
    as.vector(variable_distance),
    as.vector(cophenetic_distance),
    method = "pearson"
  )
  
  hca_sensitivity_results[[scenario_name]] <- list(
    variable_profiles = variable_profiles,
    distance_matrix = variable_distance,
    hca_result = hca_result,
    membership = group_membership,
    silhouette_result = silhouette_result,
    cophenetic_correlation = cophenetic_correlation
  )
  
  hca_diagnostic_records[[scenario_index]] <- data.frame(
    Scenario = scenario_name,
    Number_of_variables = length(group_membership),
    Number_of_groups_G = selected_number_of_groups,
    Mean_silhouette_width = mean(silhouette_widths),
    Minimum_silhouette_width = min(silhouette_widths),
    Number_of_negative_silhouettes = sum(silhouette_widths < 0),
    Cophenetic_correlation = cophenetic_correlation,
    stringsAsFactors = FALSE
  )
  
  hca_membership_records[[scenario_index]] <- data.frame(
    Scenario = scenario_name,
    Variable = names(group_membership),
    Group = as.integer(group_membership),
    stringsAsFactors = FALSE
  )
  
  group_sizes <- table(
    group_membership
  )
  
  hca_group_size_records[[scenario_index]] <- data.frame(
    Scenario = scenario_name,
    Group = as.integer(names(group_sizes)),
    Number_of_variables = as.integer(group_sizes),
    stringsAsFactors = FALSE
  )
  
  scenario_index <- scenario_index + 1L
}


# Combine HCA results

hca_diagnostic_summary <- do.call(
  rbind,
  hca_diagnostic_records
)

rownames(hca_diagnostic_summary) <- NULL

hca_membership_summary <- do.call(
  rbind,
  hca_membership_records
)

rownames(hca_membership_summary) <- NULL

hca_group_size_summary_all_scenarios <- do.call(
  rbind,
  hca_group_size_records
)

rownames(hca_group_size_summary_all_scenarios) <- NULL


# Compare every partition with the primary uniform_k8 partition

reference_membership <- hca_sensitivity_results[[primary_scenario_name]]$membership

hca_comparison_records <- list()

scenario_index <- 1L

for (scenario_name in names(hca_sensitivity_results)) {
  
  comparison_membership <- hca_sensitivity_results[[scenario_name]]$membership
  
  common_variables <- names(reference_membership)[
    names(reference_membership) %in% names(comparison_membership)
  ]
  
  reference_common <- reference_membership[
    common_variables
  ]
  
  comparison_common <- comparison_membership[
    common_variables
  ]
  
  ari_value <- adjusted_rand_index(
    reference_common,
    comparison_common
  )
  
  label_alignment <- align_group_labels(
    reference_membership,
    comparison_membership
  )
  
  changed_variable_text <- if (
    length(label_alignment$changed_variables) == 0L
  ) {
    "None"
  } else {
    paste(
      label_alignment$changed_variables,
      collapse = "; "
    )
  }
  
  hca_comparison_records[[scenario_index]] <- data.frame(
    Scenario = scenario_name,
    Reference_scenario = primary_scenario_name,
    Number_of_common_variables = length(common_variables),
    Adjusted_Rand_Index = ari_value,
    Membership_agreement_percent = label_alignment$agreement_percent,
    Number_of_changed_variables = length(label_alignment$changed_variables),
    Changed_variables = changed_variable_text,
    stringsAsFactors = FALSE
  )
  
  scenario_index <- scenario_index + 1L
}

hca_stability_comparison <- do.call(
  rbind,
  hca_comparison_records
)

rownames(hca_stability_comparison) <- NULL


# Export all HCA sensitivity tables

write.csv(
  hca_diagnostic_summary,
  file.path(
    tables_directory,
    "hca_diagnostics_all_scenarios_G3.csv"
  ),
  row.names = FALSE
)

write.csv(
  hca_membership_summary,
  file.path(
    tables_directory,
    "hca_memberships_all_scenarios_G3.csv"
  ),
  row.names = FALSE
)

write.csv(
  hca_group_size_summary_all_scenarios,
  file.path(
    tables_directory,
    "hca_group_sizes_all_scenarios_G3.csv"
  ),
  row.names = FALSE
)

write.csv(
  hca_stability_comparison,
  file.path(
    tables_directory,
    "hca_stability_vs_uniform_k8_G3.csv"
  ),
  row.names = FALSE
)


# Display rounded diagnostic results

hca_diagnostic_display <- hca_diagnostic_summary

hca_diagnostic_display$Mean_silhouette_width <- round(
  hca_diagnostic_display$Mean_silhouette_width,
  4
)

hca_diagnostic_display$Minimum_silhouette_width <- round(
  hca_diagnostic_display$Minimum_silhouette_width,
  4
)

hca_diagnostic_display$Cophenetic_correlation <- round(
  hca_diagnostic_display$Cophenetic_correlation,
  4
)

hca_stability_display <- hca_stability_comparison

hca_stability_display$Adjusted_Rand_Index <- round(
  hca_stability_display$Adjusted_Rand_Index,
  4
)

hca_stability_display$Membership_agreement_percent <- round(
  hca_stability_display$Membership_agreement_percent,
  2
)

cat(
  "HCA sensitivity analysis was completed successfully for G =",
  selected_number_of_groups,
  "\n"
)

cat(
  "\nHCA diagnostic results:\n"
)

print(
  hca_diagnostic_display,
  row.names = FALSE
)

cat(
  "\nHCA stability relative to uniform_k8:\n"
)

print(
  hca_stability_display,
  row.names = FALSE
)

cat(
  "\nGroup sizes for all sensitivity scenarios:\n"
)

print(
  hca_group_size_summary_all_scenarios,
  row.names = FALSE
)

# Step 12 — Leave-one-time-point-out sensitivity analysis

if (!exists("pca_sensitivity_results")) {
  stop(
    "The object 'pca_sensitivity_results' was not found. Run Step 7 first."
  )
}

if (!exists("primary_scenario_name")) {
  stop(
    "The object 'primary_scenario_name' was not found."
  )
}

if (!exists("selected_number_of_groups")) {
  stop(
    "The selected number of HCA groups was not found."
  )
}

if (is.na(selected_number_of_groups)) {
  stop(
    "The selected number of HCA groups G has not been defined."
  )
}

if (!exists("compare_pca_configurations")) {
  stop(
    "The function 'compare_pca_configurations' was not found. Run Step 8 first."
  )
}

if (!exists("adjusted_rand_index")) {
  stop(
    "The function 'adjusted_rand_index' was not found. Run Step 11 first."
  )
}

if (!exists("align_group_labels")) {
  stop(
    "The function 'align_group_labels' was not found. Run Step 11 first."
  )
}

if (!requireNamespace("cluster", quietly = TRUE)) {
  stop(
    "The 'cluster' package is required."
  )
}


# Extract the complete primary matrix

complete_primary_data <- pca_sensitivity_results[[primary_scenario_name]]$original_data

if (nrow(complete_primary_data) != 7L) {
  stop(
    "The primary matrix does not contain the expected seven time points."
  )
}

time_points <- rownames(
  complete_primary_data
)

if (is.null(time_points)) {
  stop(
    "Time-point labels were not found in the primary matrix."
  )
}


# Extract the complete-reference PCA coordinates

complete_reference_pca <- pca_sensitivity_results[[primary_scenario_name]]$pca_result

reference_time_coordinates_2d <- as.matrix(
  complete_reference_pca$ind$coord[, 1:2, drop = FALSE]
)

reference_time_coordinates_3d <- as.matrix(
  complete_reference_pca$ind$coord[, 1:3, drop = FALSE]
)

reference_variable_coordinates_2d <- as.matrix(
  complete_reference_pca$var$coord[, 1:2, drop = FALSE]
)

reference_variable_coordinates_3d <- as.matrix(
  complete_reference_pca$var$coord[, 1:3, drop = FALSE]
)


# Extract the complete-reference HCA partition

reference_hca_membership <- cutree(
  primary_hca_result,
  k = selected_number_of_groups
)


# Prepare result containers

leave_one_time_out_results <- list()
loto_pca_records <- list()
loto_hca_records <- list()
loto_membership_records <- list()

record_index <- 1L


# Repeat the analysis after removing each time point

for (omitted_time in time_points) {
  
  reduced_data <- complete_primary_data[
    rownames(complete_primary_data) != omitted_time,
    ,
    drop = FALSE
  ]
  
  if (nrow(reduced_data) != 6L) {
    stop(
      paste(
        "Unexpected number of retained time points after excluding:",
        omitted_time
      )
    )
  }
  
  
  # Check variable variability after time-point exclusion
  
  reduced_standard_deviations <- apply(
    reduced_data,
    2,
    sd
  )
  
  if (
    any(!is.finite(reduced_standard_deviations)) ||
    any(reduced_standard_deviations == 0)
  ) {
    stop(
      paste(
        "A constant or invalid variable was detected after excluding:",
        omitted_time
      )
    )
  }
  
  
  # Re-standardize the reduced matrix
  
  reduced_scaled_data <- scale(
    reduced_data,
    center = TRUE,
    scale = TRUE
  )
  
  colnames(reduced_scaled_data) <- sub(
    "^D",
    "Z",
    colnames(reduced_scaled_data)
  )
  
  # Recalculate PCA
  
  maximum_number_of_components <- min(
    nrow(reduced_scaled_data) - 1L,
    ncol(reduced_scaled_data)
  )
  
  reduced_pca <- FactoMineR::PCA(
    as.data.frame(reduced_scaled_data),
    scale.unit = FALSE,
    ncp = maximum_number_of_components,
    graph = FALSE
  )
  
  
  # Extract reduced PCA configurations
  
  reduced_time_coordinates_2d <- as.matrix(
    reduced_pca$ind$coord[, 1:2, drop = FALSE]
  )
  
  reduced_time_coordinates_3d <- as.matrix(
    reduced_pca$ind$coord[, 1:3, drop = FALSE]
  )
  
  reduced_variable_coordinates_2d <- as.matrix(
    reduced_pca$var$coord[, 1:2, drop = FALSE]
  )
  
  reduced_variable_coordinates_3d <- as.matrix(
    reduced_pca$var$coord[, 1:3, drop = FALSE]
  )
  
  
  # Compare PCA configurations with the complete analysis
  
  time_comparison_2d <- compare_pca_configurations(
    reference_time_coordinates_2d,
    reduced_time_coordinates_2d
  )
  
  time_comparison_3d <- compare_pca_configurations(
    reference_time_coordinates_3d,
    reduced_time_coordinates_3d
  )
  
  variable_comparison_2d <- compare_pca_configurations(
    reference_variable_coordinates_2d,
    reduced_variable_coordinates_2d
  )
  
  variable_comparison_3d <- compare_pca_configurations(
    reference_variable_coordinates_3d,
    reduced_variable_coordinates_3d
  )
  
  
  # Recalculate Ward HCA
  
  reduced_variable_profiles <- t(
    reduced_scaled_data
  )
  
  reduced_variable_distance <- dist(
    reduced_variable_profiles,
    method = "euclidean"
  )
  
  reduced_hca <- hclust(
    reduced_variable_distance,
    method = "ward.D2"
  )
  
  reduced_membership <- cutree(
    reduced_hca,
    k = selected_number_of_groups
  )
  
  
  # Calculate HCA diagnostics
  
  reduced_silhouette <- cluster::silhouette(
    reduced_membership,
    reduced_variable_distance
  )
  
  reduced_silhouette_widths <- reduced_silhouette[, "sil_width"]
  
  reduced_cophenetic_distance <- cophenetic(
    reduced_hca
  )
  
  reduced_cophenetic_correlation <- cor(
    as.vector(reduced_variable_distance),
    as.vector(reduced_cophenetic_distance),
    method = "pearson"
  )
  
  
  # Compare the reduced HCA partition with the complete partition
  
  ari_value <- adjusted_rand_index(
    reference_hca_membership,
    reduced_membership
  )
  
  label_alignment <- align_group_labels(
    reference_hca_membership,
    reduced_membership
  )
  
  changed_variable_text <- if (
    length(label_alignment$changed_variables) == 0L
  ) {
    "None"
  } else {
    paste(
      label_alignment$changed_variables,
      collapse = "; "
    )
  }
  
  
  # Store detailed results
  
  leave_one_time_out_results[[omitted_time]] <- list(
    omitted_time = omitted_time,
    reduced_data = reduced_data,
    scaled_data = reduced_scaled_data,
    pca_result = reduced_pca,
    hca_result = reduced_hca,
    membership = reduced_membership,
    silhouette_result = reduced_silhouette,
    cophenetic_correlation = reduced_cophenetic_correlation
  )
  
  
  # Store PCA summary
  
  loto_pca_records[[record_index]] <- data.frame(
    Omitted_time = omitted_time,
    Number_of_retained_times = nrow(reduced_data),
    
    PC1_percent = reduced_pca$eig[1, 2],
    PC2_percent = reduced_pca$eig[2, 2],
    PC3_percent = reduced_pca$eig[3, 2],
    
    PC1_PC2_cumulative_percent = reduced_pca$eig[2, 3],
    PC1_PC2_PC3_cumulative_percent = reduced_pca$eig[3, 3],
    
    Time_distance_correlation_2D = time_comparison_2d$distance_correlation,
    Time_Procrustes_similarity_2D = time_comparison_2d$procrustes_similarity,
    Time_residual_RMSE_2D = time_comparison_2d$normalized_residual_rmse,
    
    Time_distance_correlation_3D = time_comparison_3d$distance_correlation,
    Time_Procrustes_similarity_3D = time_comparison_3d$procrustes_similarity,
    Time_residual_RMSE_3D = time_comparison_3d$normalized_residual_rmse,
    
    Variable_distance_correlation_2D = variable_comparison_2d$distance_correlation,
    Variable_Procrustes_similarity_2D = variable_comparison_2d$procrustes_similarity,
    Variable_residual_RMSE_2D = variable_comparison_2d$normalized_residual_rmse,
    
    Variable_distance_correlation_3D = variable_comparison_3d$distance_correlation,
    Variable_Procrustes_similarity_3D = variable_comparison_3d$procrustes_similarity,
    Variable_residual_RMSE_3D = variable_comparison_3d$normalized_residual_rmse,
    
    stringsAsFactors = FALSE
  )
  
  
  # Store HCA summary
  
  loto_hca_records[[record_index]] <- data.frame(
    Omitted_time = omitted_time,
    Number_of_variables = length(reduced_membership),
    Number_of_groups_G = selected_number_of_groups,
    
    Mean_silhouette_width = mean(reduced_silhouette_widths),
    Minimum_silhouette_width = min(reduced_silhouette_widths),
    Number_of_negative_silhouettes = sum(reduced_silhouette_widths < 0),
    
    Cophenetic_correlation = reduced_cophenetic_correlation,
    Adjusted_Rand_Index = ari_value,
    Membership_agreement_percent = label_alignment$agreement_percent,
    
    Number_of_changed_variables = length(
      label_alignment$changed_variables
    ),
    
    Changed_variables = changed_variable_text,
    
    stringsAsFactors = FALSE
  )
  
  
  # Store aligned memberships
  
  loto_membership_records[[record_index]] <- data.frame(
    Omitted_time = omitted_time,
    Variable = names(label_alignment$reference_membership),
    Reference_group = as.integer(
      label_alignment$reference_membership
    ),
    Leave_one_out_group = as.integer(
      label_alignment$aligned_membership
    ),
    Group_changed = (
      label_alignment$reference_membership !=
        label_alignment$aligned_membership
    ),
    stringsAsFactors = FALSE
  )
  
  record_index <- record_index + 1L
}


# Combine all results

loto_pca_summary <- do.call(
  rbind,
  loto_pca_records
)

rownames(loto_pca_summary) <- NULL

loto_hca_summary <- do.call(
  rbind,
  loto_hca_records
)

rownames(loto_hca_summary) <- NULL

loto_membership_summary <- do.call(
  rbind,
  loto_membership_records
)

rownames(loto_membership_summary) <- NULL


# Export the leave-one-time-out tables

write.csv(
  loto_pca_summary,
  file.path(
    tables_directory,
    "leave_one_time_out_pca_sensitivity.csv"
  ),
  row.names = FALSE
)

write.csv(
  loto_hca_summary,
  file.path(
    tables_directory,
    "leave_one_time_out_hca_sensitivity_G3.csv"
  ),
  row.names = FALSE
)

write.csv(
  loto_membership_summary,
  file.path(
    tables_directory,
    "leave_one_time_out_hca_memberships_G3.csv"
  ),
  row.names = FALSE
)


# Prepare rounded display tables

loto_pca_display <- loto_pca_summary

numeric_pca_columns <- vapply(
  loto_pca_display,
  is.numeric,
  logical(1)
)

loto_pca_display[numeric_pca_columns] <- lapply(
  loto_pca_display[numeric_pca_columns],
  round,
  digits = 4
)

loto_hca_display <- loto_hca_summary

numeric_hca_columns <- vapply(
  loto_hca_display,
  is.numeric,
  logical(1)
)

loto_hca_display[numeric_hca_columns] <- lapply(
  loto_hca_display[numeric_hca_columns],
  round,
  digits = 4
)


cat(
  "Leave-one-time-point-out sensitivity analysis was completed successfully.\n"
)

cat(
  "\nPCA leave-one-time-point-out results:\n"
)

print(
  loto_pca_display,
  row.names = FALSE
)

cat(
  "\nHCA leave-one-time-point-out results:\n"
)

print(
  loto_hca_display,
  row.names = FALSE
)

# Step 13 — Summarize variable-level HCA stability

if (!exists("loto_membership_summary")) {
  stop(
    "The object 'loto_membership_summary' was not found. Run Step 12 first."
  )
}

if (!exists("tables_directory")) {
  stop(
    "The object 'tables_directory' was not found. Run Step 1 first."
  )
}

dir.create(
  tables_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

variable_names <- sort(
  unique(loto_membership_summary$Variable)
)

variable_stability_records <- list()

record_index <- 1L

for (variable_name in variable_names) {
  
  variable_results <- loto_membership_summary[
    loto_membership_summary$Variable == variable_name,
    ,
    drop = FALSE
  ]
  
  changed_results <- variable_results[
    variable_results$Group_changed,
    ,
    drop = FALSE
  ]
  
  number_of_changes <- nrow(
    changed_results
  )
  
  omitted_times_causing_change <- if (number_of_changes == 0L) {
    "None"
  } else {
    paste(
      changed_results$Omitted_time,
      collapse = "; "
    )
  }
  
  group_transitions <- if (number_of_changes == 0L) {
    "None"
  } else {
    paste(
      paste0(
        changed_results$Omitted_time,
        ": ",
        changed_results$Reference_group,
        " -> ",
        changed_results$Leave_one_out_group
      ),
      collapse = "; "
    )
  }
  
  sensitivity_class <- if (number_of_changes == 0L) {
    "Very stable"
  } else if (number_of_changes == 1L) {
    "Moderately sensitive"
  } else {
    "Sensitive"
  }
  
  variable_stability_records[[record_index]] <- data.frame(
    Variable = variable_name,
    Reference_group = unique(variable_results$Reference_group),
    Number_of_omissions_tested = nrow(variable_results),
    Number_of_group_changes = number_of_changes,
    Stability_percent = 100 * (
      1 - number_of_changes / nrow(variable_results)
    ),
    Sensitivity_class = sensitivity_class,
    Omitted_times_causing_change = omitted_times_causing_change,
    Group_transitions = group_transitions,
    stringsAsFactors = FALSE
  )
  
  record_index <- record_index + 1L
}

variable_hca_stability_summary <- do.call(
  rbind,
  variable_stability_records
)

rownames(variable_hca_stability_summary) <- NULL

variable_hca_stability_summary <- variable_hca_stability_summary[
  order(
    -variable_hca_stability_summary$Number_of_group_changes,
    variable_hca_stability_summary$Variable
  ),
  ,
  drop = FALSE
]

# Define a simple and valid export path

variable_stability_output_path <- file.path(
  tables_directory,
  "leave_one_time_out_variable_hca_stability_G3.csv"
)

# Export the summary table

write.csv(
  variable_hca_stability_summary,
  file = variable_stability_output_path,
  row.names = FALSE
)

# Prepare the displayed table

variable_hca_stability_display <- variable_hca_stability_summary

variable_hca_stability_display$Stability_percent <- round(
  variable_hca_stability_display$Stability_percent,
  2
)

cat(
  "Variable-level HCA stability summary was completed successfully.\n"
)

cat(
  "Output file:",
  variable_stability_output_path,
  "\n"
)

cat(
  "Output file found:",
  file.exists(variable_stability_output_path),
  "\n"
)

print(
  variable_hca_stability_display,
  row.names = FALSE
)

# Step 13B — Compare candidate G values after each time omission

if (!exists("leave_one_time_out_results")) {
  stop(
    "The object 'leave_one_time_out_results' was not found. Run Step 12 first."
  )
}

if (!exists("tables_directory")) {
  stop(
    "The object 'tables_directory' was not found. Run Step 1 first."
  )
}

if (!requireNamespace("cluster", quietly = TRUE)) {
  stop(
    "The 'cluster' package is required."
  )
}

dir.create(
  tables_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# Candidate numbers of HCA groups

candidate_G_values_loto <- 2:4

loto_G_diagnostic_records <- list()

record_index <- 1L

# Repeat the comparison for every omitted time point

for (omitted_time in names(leave_one_time_out_results)) {
  
  reduced_scaled_data <- leave_one_time_out_results[[omitted_time]]$scaled_data
  
  if (anyNA(reduced_scaled_data)) {
    stop(
      paste(
        "Missing values were detected after excluding:",
        omitted_time
      )
    )
  }
  
  # Variables are clustered according to their reduced temporal profiles
  
  reduced_variable_profiles <- t(
    reduced_scaled_data
  )
  
  reduced_variable_distance <- dist(
    reduced_variable_profiles,
    method = "euclidean"
  )
  
  reduced_hca_result <- hclust(
    reduced_variable_distance,
    method = "ward.D2"
  )
  
  for (G in candidate_G_values_loto) {
    
    reduced_membership <- cutree(
      reduced_hca_result,
      k = G
    )
    
    reduced_silhouette <- cluster::silhouette(
      reduced_membership,
      reduced_variable_distance
    )
    
    silhouette_widths <- reduced_silhouette[, "sil_width"]
    
    group_sizes <- table(
      reduced_membership
    )
    
    loto_G_diagnostic_records[[record_index]] <- data.frame(
      Omitted_time = omitted_time,
      Number_of_groups_G = G,
      Mean_silhouette_width = mean(silhouette_widths),
      Minimum_silhouette_width = min(silhouette_widths),
      Number_of_negative_silhouettes = sum(silhouette_widths < 0),
      Smallest_group_size = min(group_sizes),
      Largest_group_size = max(group_sizes),
      stringsAsFactors = FALSE
    )
    
    record_index <- record_index + 1L
  }
}

# Combine all candidate-G diagnostics

loto_G_diagnostic_summary <- do.call(
  rbind,
  loto_G_diagnostic_records
)

rownames(loto_G_diagnostic_summary) <- NULL

# Identify the highest-silhouette G for each omitted time

omitted_time_values <- unique(
  loto_G_diagnostic_summary$Omitted_time
)

best_G_by_omission_records <- list()

record_index <- 1L

for (omitted_time in omitted_time_values) {
  
  diagnostic_table <- loto_G_diagnostic_summary[
    loto_G_diagnostic_summary$Omitted_time == omitted_time,
    ,
    drop = FALSE
  ]
  
  best_row_index <- which.max(
    diagnostic_table$Mean_silhouette_width
  )
  
  best_G_by_omission_records[[record_index]] <- diagnostic_table[
    best_row_index,
    ,
    drop = FALSE
  ]
  
  record_index <- record_index + 1L
}

best_G_by_omission <- do.call(
  rbind,
  best_G_by_omission_records
)

rownames(best_G_by_omission) <- NULL

# Extract the G = 3 results for direct verification

selected_G_loto_summary <- loto_G_diagnostic_summary[
  loto_G_diagnostic_summary$Number_of_groups_G == selected_number_of_groups,
  ,
  drop = FALSE
]

rownames(selected_G_loto_summary) <- NULL

# Define output paths

loto_G_diagnostics_output_path <- file.path(
  tables_directory,
  "leave_one_time_out_candidate_G_diagnostics.csv"
)

best_G_output_path <- file.path(
  tables_directory,
  "leave_one_time_out_best_silhouette_G.csv"
)

selected_G_output_path <- file.path(
  tables_directory,
  "leave_one_time_out_selected_G3_diagnostics.csv"
)

# Export the tables

write.csv(
  loto_G_diagnostic_summary,
  file = loto_G_diagnostics_output_path,
  row.names = FALSE
)

write.csv(
  best_G_by_omission,
  file = best_G_output_path,
  row.names = FALSE
)

write.csv(
  selected_G_loto_summary,
  file = selected_G_output_path,
  row.names = FALSE
)

# Prepare rounded tables for display

loto_G_diagnostic_display <- loto_G_diagnostic_summary

loto_G_diagnostic_display$Mean_silhouette_width <- round(
  loto_G_diagnostic_display$Mean_silhouette_width,
  4
)

loto_G_diagnostic_display$Minimum_silhouette_width <- round(
  loto_G_diagnostic_display$Minimum_silhouette_width,
  4
)

best_G_by_omission_display <- best_G_by_omission

best_G_by_omission_display$Mean_silhouette_width <- round(
  best_G_by_omission_display$Mean_silhouette_width,
  4
)

best_G_by_omission_display$Minimum_silhouette_width <- round(
  best_G_by_omission_display$Minimum_silhouette_width,
  4
)

selected_G_loto_display <- selected_G_loto_summary

selected_G_loto_display$Mean_silhouette_width <- round(
  selected_G_loto_display$Mean_silhouette_width,
  4
)

selected_G_loto_display$Minimum_silhouette_width <- round(
  selected_G_loto_display$Minimum_silhouette_width,
  4
)

# Display the results

cat(
  "Candidate-G diagnostics were completed for all time omissions.\n"
)

cat(
  "\nDiagnostics for G = 2, G = 3, and G = 4:\n"
)

print(
  loto_G_diagnostic_display,
  row.names = FALSE
)

cat(
  "\nHighest-silhouette G for each omitted time:\n"
)

print(
  best_G_by_omission_display,
  row.names = FALSE
)

cat(
  "\nDiagnostics for the selected G =",
  selected_number_of_groups,
  "after each time omission:\n"
)

print(
  selected_G_loto_display,
  row.names = FALSE
)

# Verify the exported files

cat(
  "\nCandidate-G output file:",
  loto_G_diagnostics_output_path,
  "\n"
)

cat(
  "Candidate-G file found:",
  file.exists(loto_G_diagnostics_output_path),
  "\n"
)

cat(
  "Best-G output file:",
  best_G_output_path,
  "\n"
)

cat(
  "Best-G file found:",
  file.exists(best_G_output_path),
  "\n"
)

cat(
  "Selected-G output file:",
  selected_G_output_path,
  "\n"
)

cat(
  "Selected-G file found:",
  file.exists(selected_G_output_path),
  "\n"
)
