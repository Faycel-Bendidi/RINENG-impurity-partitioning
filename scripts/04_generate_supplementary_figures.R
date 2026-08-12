source(file.path("scripts", "_helpers.R"))
assert_repository_root()

required_inputs <- c(
  file.path("results", "tables", "pca_individual_coordinates_uniform_k8.csv"),
  file.path("results", "tables", "hca_elbow_diagnostics_uniform_k8.csv"),
  file.path("results", "tables", "hca_silhouette_diagnostics_uniform_k8.csv"),
  file.path("data", "processed", "knn_cross_validation_rmse_mae_delta.csv"),
  file.path("data", "processed", "local_k_sensitivity_summary.csv"),
  file.path("results", "tables", "leave_one_time_out_pca_sensitivity.csv"),
  file.path("results", "tables", "leave_one_time_out_hca_sensitivity_G3.csv"),
  file.path("results", "tables", "leave_one_time_out_variable_hca_stability_G3.csv")
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    paste("Missing supplementary-figure input files:", paste(missing_inputs, collapse = ", ")),
    call. = FALSE
  )
}

for (script_name in c("Figure_S1.R", "Figure_S2.R", "Figure_S3.R", "Figure_S4.R")) {
  message("Running ", script_name)
  source(
    file.path("scripts", "figures", script_name),
    local = new.env(parent = globalenv())
  )
}
