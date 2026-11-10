# Validate winoddsEIF 0.1.0 against the original semi-competing-risk code.
# Corrected loader (original script functions are evaluated explicitly).
# This comparison uses identical data, cross-fitting seeds, ridge penalty,
# and conditional-projection grids. It does NOT run the simulation study.
# Run from RStudio with: source(file.choose())

library(winoddsEIF)

original_path <- system.file(
  "reproducibility", "semicomp_ridge_all_scenarios.R",
  package = "winoddsEIF"
)
if (!nzchar(original_path)) {
  stop("Original semi-competing-risk script not found in installed package")
}

# Important: The original script finishes by automatically running thousands
# of simulation replications. Load only definitions, not its final run commands.
lines <- readLines(original_path, warn = FALSE)
run_start <- which(grepl("^semi_sanity_check\\(\\)$", trimws(lines)))
if (length(run_start) != 1L || run_start < 2L) {
  stop("Could not safely locate the original script's run commands")
}
# Load just the function definitions in an isolated environment.
# Do not execute the simulation commands at the end of the source file.
original <- new.env(parent = globalenv())
definitions <- parse(
  text = paste(lines[seq_len(run_start - 1L)], collapse = "\n"),
  keep.source = FALSE
)
for (expression in definitions) {
  eval(expression, envir = original)
}
for (name in c("semi_parameters", "semi_calibrate_censoring",
               "semi_generate", "semi_one_step")) {
  if (!exists(name, envir = original, inherits = FALSE) ||
      !is.function(get(name, envir = original, inherits = FALSE))) {
    stop("Failed to load original function: ", name)
  }
}
cat("Original simulation functions loaded successfully.\n")

# For code-equivalence checks, 16 grid points are enough and fast.
# Change to 256 when you want to check the paper's integration setting.
M_Q <- 16L
N <- 240L
FOLDS <- 5L
FOLD_SEED <- 2027L
SCENARIOS_TO_CHECK <- c("moderate", "strongest")
TOL <- 1e-8

validate_semicomp_package <- function() {
  for (scenario in SCENARIOS_TO_CHECK) {
  cat("\n=== Scenario:", scenario, "===\n")
  par <- original$semi_parameters(scenario)
  c0 <- original$semi_calibrate_censoring(
    par, M = 10000L, seed = 20261103L
  )
  set.seed(2031L)
  d_original <- original$semi_generate(N, par, c0)
  attr(d_original, "tau") <- par$tau

  # Original estimator, with all nuisance fits based on training folds.
  fit_original <- original$semi_one_step(
    d_original,
    V = FOLDS, M_q = M_Q, seed = FOLD_SEED,
    piA = par$piA, mu_method = "model", oracle = FALSE
  )

  # Same observed data, reshaped only to the public package interface.
  d_package <- data.frame(
    A = d_original$A,
    X = d_original$X,
    time1 = d_original$T1tilde,
    status1 = d_original$Delta1,
    time2 = d_original$T2tilde,
    status2 = d_original$Delta2
  )
  fit_package <- winprob_semicomp(
    d_package,
    tau = par$tau,
    treatment_prob = par$piA,
    folds = FOLDS,
    seed = FOLD_SEED,
    m_q = M_Q,
    mu_method = "model"
  )

  reference <- c(
    psi = fit_original$full,
    se_psi = fit_original$se_full,
    plugin = fit_original$plugin,
    ipcw = fit_original$ipcw,
    se_ipcw = fit_original$se_ipcw,
    min_fitted_G = fit_original$minG
  )
  package <- c(
    psi = fit_package$psi,
    se_psi = fit_package$se_psi,
    plugin = fit_package$plugin,
    ipcw = fit_package$ipcw,
    se_ipcw = fit_package$se_ipcw,
    min_fitted_G = fit_package$diagnostics$min_fitted_G
  )
  comparison <- data.frame(
    quantity = names(reference),
    original = as.numeric(reference),
    package = as.numeric(package),
    absolute_difference = abs(as.numeric(reference - package)),
    row.names = NULL
  )
  print(comparison, row.names = FALSE, digits = 12)
  if (!isTRUE(all.equal(reference, package, tolerance = TOL))) {
    stop("Original and package disagree in scenario: ", scenario)
  }
  cat("PASS: all six quantities agree within tolerance ", TOL, "\n", sep = "")
}

  cat("\nAll semi-competing-risk comparisons passed.\n")
  invisible(TRUE)
}

# Run the entire validation; the success message above is printed ONLY
# if every comparison completes without error.
validate_semicomp_package()
