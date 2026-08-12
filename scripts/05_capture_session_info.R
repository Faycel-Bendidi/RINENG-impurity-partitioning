source(file.path("scripts", "_helpers.R"))
assert_repository_root()

packages <- c(
  "FactoMineR", "cluster", "factoextra", "ggplot2", "readxl", "tidyr",
  "dplyr", "readr", "patchwork", "ggrepel", "dendextend", "scales"
)
load_required_packages(packages)

output_path <- "sessionInfo_R.txt"
sink(output_path)
print(sessionInfo())
cat("
Explicit package versions:
")
for (package_name in packages) {
  cat(package_name, as.character(packageVersion(package_name)), "
")
}
sink()
message("R session information written to ", output_path)
