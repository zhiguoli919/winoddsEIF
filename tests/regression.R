# Numerical regression tests for winoddsEIF.
# R CMD check automatically runs each .R file in tests/.
# The original simulation drivers are NOT run: the bundled validation
# scripts load only function definitions and use two small fixed datasets.

library(winoddsEIF)

run_regression <- function(filename) {
  script <- system.file(
    "validation", filename,
    package = "winoddsEIF",
    mustWork = TRUE
  )
  cat("\nRunning numerical regression test: ", filename, "\n", sep = "")
  # An error in the validation script propagates to R CMD check.
  sys.source(script, envir = new.env(parent = globalenv()))
  invisible(TRUE)
}

run_regression("validate_winoddsEIF_scalar.R")
run_regression("validate_winoddsEIF_semicomp.R")

cat("\nBoth numerical regression tests passed.\n")
