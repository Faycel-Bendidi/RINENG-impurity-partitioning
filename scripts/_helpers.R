require_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0L) {
    stop(
      paste0(
        "Missing R packages: ",
        paste(missing_packages, collapse = ", "),
        ". Install them before running this script."
      ),
      call. = FALSE
    )
  }

  invisible(packages)
}

load_required_packages <- function(packages) {
  require_packages(packages)
  invisible(
    lapply(packages, library, character.only = TRUE)
  )
}

assert_repository_root <- function() {
  required_paths <- c("README.md", "data", "scripts")
  missing_paths <- required_paths[!file.exists(required_paths)]

  if (length(missing_paths) > 0L) {
    stop(
      "Run this script from the repository root.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
