source(file.path("scripts", "_helpers.R"))
assert_repository_root()

required_inputs <- c(
  file.path("data", "external", "Supplementary_Data.xlsx"),
  file.path("data", "external", "Granulo_laser_30.xlsx"),
  file.path("results", "tables", "pca_variable_coordinates_uniform_k8.csv"),
  file.path("data", "processed", "partition_coefficients_standardized_Z.csv")
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    paste("Missing figure input files:", paste(missing_inputs, collapse = ", ")),
    call. = FALSE
  )
}

figure_scripts <- c(
  "Figure_3.R",
  "Figure_4.R",
  "Figure_5.R",
  "Figure_6.R",
  "Figure_7.R",
  "Figure_8.R",
  "Figure_11.R"
)

for (script_name in figure_scripts) {
  message("Running ", script_name)
  source(
    file.path("scripts", "figures", script_name),
    local = new.env(parent = globalenv())
  )
}
