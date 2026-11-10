# winoddsEIF

**Cross-fitted estimation of marginal win probability and win odds with censored time-to-event outcomes**

`winoddsEIF` is a research R package accompanying the manuscript *Semiparametric Inference for Marginal Win Odds with Censored Time-to-Event Outcomes*. It implements cross-fitted one-step estimators for (i) scalar right-censored survival outcomes and (ii) prioritized semi-competing-risk outcomes. It also provides the original simulation programs, numerical comparison tests, and the separate penalized-spline analysis of a breast cancer trial.

## Installation

From a terminal in the parent directory of the package source folder:

```sh
R CMD INSTALL winoddsEIF
```

Alternatively, install the source archive from R (using the archive's actual path):

```r
install.packages("winoddsEIF_0.1.0.tar.gz", repos = NULL, type = "source")
library(winoddsEIF)
```

## Scalar right-censored survival outcomes

The input data frame has one row per participant and columns `A` (0 = control, 1 = treatment), `X` (one numeric baseline covariate), `time` (observed event/censoring time, administratively restricted at `tau`), and `status` (1 if the event occurred and was observed by `tau`, 0 otherwise).

```r
library(winoddsEIF)

set.seed(3)
n <- 240
tau <- 5
A <- rbinom(n, 1, 0.5)
X <- runif(n, -1, 1)
T <- rexp(n, rate = 0.10 * exp(-0.4 * A + 0.7 * X))
C <- rexp(n, rate = 0.09 * exp(0.5 * X))

d <- data.frame(
  A = A,
  X = X,
  time = pmin(T, C, tau),
  status = as.integer(T <= C & T <= tau)
)

fit <- winprob_scalar(d, tau = tau, treatment_prob = 0.5,
                      folds = 5, seed = 123)
print(fit)
fit$diagnostics
```

`status = 0` includes both censoring before `tau` and event-free follow-up through `tau`. The latter yields a fully observed restricted outcome at the horizon.

## Prioritized semi-competing risks

Here the **terminal event has higher priority** than the nonterminal event, with a common restriction horizon `tau`.

| Column | Meaning |
|---|---|
| `A` | 0 = control, 1 = treatment |
| `X` | One numeric baseline covariate |
| `time1` | Minimum of terminal-event time, censoring time, and follow-up end |
| `status1` | 1 if a terminal event was observed, 0 otherwise |
| `time2` | Minimum of nonterminal-event time, terminal-event time, censoring time, and follow-up end |
| `status2` | 1 if a nonterminal event was observed before the terminal event, 0 otherwise |

A complete synthetic-data example:

```r
library(winoddsEIF)

set.seed(2027)
n <- 240L
tau <- 5
A <- rbinom(n, 1, 0.5)
X <- runif(n, -1, 1)

T1 <- rexp(n, rate = 0.10 * exp(-0.35 * A + 0.7 * X))
T2 <- rexp(n, rate = 0.20 * exp(-0.30 * A + 0.4 * X))
C  <- rexp(n, rate = 0.08 * exp(0.2 * X))

d <- data.frame(
  A = A, X = X,
  time1 = pmin(T1, C, tau),
  status1 = as.integer(T1 <= C & T1 <= tau),
  time2 = pmin(T2, T1, C, tau),
  status2 = as.integer(T2 < T1 & T2 <= C & T2 <= tau)
)

fit <- winprob_semicomp(d, tau = tau, treatment_prob = 0.5,
                        folds = 5, seed = 2027, m_q = 16)
print(fit)
fit$psi
fit$win_odds
fit$ci_win_odds
fit$diagnostics
```

The example uses `m_q = 16` to reduce computation time. The original semi-competing-risk simulation code uses `m_q = 256`. The function restricts outcomes to `tau` even if input follow-up extends beyond it.

## Reproduce the tamoxifen real-data analysis

The package includes two standalone programs under `inst/analysis/`:

- **Main analysis:** `tamoxifen_penalized_spline_fullq.R`, with penalized cubic B-spline baseline log-hazards, five-fold cross-fitting, and a full estimated censoring augmentation.
- **Supporting analysis:** `tamoxifen_exponential_fullq.R`, with exponential working hazards.

These scripts use the publicly available `tamoxifen` dataset in the `collett` package. The main spline analysis uses age, log tumor size, and hormone-receptor status as baseline regressors and analyzes 4- and 5-year horizons. It is a **separate, reproducible analysis implementation**, not a call to the one-covariate exported `winprob_semicomp()` function.

To run the main analysis:

```r
install.packages("collett")  # Only if it is not already installed
library(winoddsEIF)

script <- system.file(
  "analysis", "tamoxifen_penalized_spline_fullq.R",
  package = "winoddsEIF", mustWork = TRUE
)
e <- new.env(parent = globalenv())
sys.source(script, envir = e)
print(e$spline_results$table)

# Record the software environment for reproducibility:
sessionInfo()
packageVersion("collett")
```

The 5-year output previously checked against the manuscript analysis was:

| Estimator | Win probability | SE | Win odds | 95% CI for win odds |
|---|---:|---:|---:|---|
| IPCW (spline censoring) | 0.536470 | 0.0153624 | 1.15736 | (1.02537, 1.30633) |
| One-step (spline full q) | 0.529994 | 0.0132644 | 1.12763 | (1.01588, 1.25168) |

**Note:** Sourcing either real-data script executes the analysis and writes output files to the current working directory. Neither analysis is run automatically by installing or loading the package. See `inst/analysis/README.md` for further details.

## Simulation code and numerical validation

The original simulation programs are provided under `inst/reproducibility/`:

- `scalar_ridge_all_scenarios.R`
- `semicomp_ridge_all_scenarios.R`

**Caution:** These standalone simulation programs contain top-level commands that can run large simulation studies and write CSV/RDS results when sourced. Review their simulation settings before running them. They are **not** executed by package installation, package examples, or the regression tests.

The numerical comparison scripts are stored under `inst/validation/`:

- `validate_winoddsEIF_scalar.R`
- `validate_winoddsEIF_semicomp.R`

The automated test file `tests/regression.R` runs both comparisons during `R CMD check`. In the moderate and strongest scenarios, all nine scalar quantities and all six semi-competing-risk quantities agreed with the original implementations within the specified tolerance of `1e-8`. The package also includes `tests/smoke.R`.

From the parent directory of the package source:

```sh
R CMD build winoddsEIF
R CMD check --no-manual winoddsEIF_0.1.0.tar.gz
```

## Scope and limitations

- The two exported estimators support **one numeric baseline covariate**, binary treatment, and a **known constant randomization probability**. This is a narrower interface than the manuscript's general semiparametric framework.
- The supplied nuisance models use ridge-regularized hazard fitting with fixed penalty `kappa = 2` to reproduce the simulation implementations; user-selectable tuning is not yet available.
- The predictable residual projections are model-based. Statistical efficiency and valid one-step inference depend on the assumptions and nuisance-rate conditions discussed in the manuscript.
- The IPCW standard error (`se_ipcw`) is a first-order model-based approximation; failures to compute it are reported in diagnostics.
- The exported estimator functions do not implement a bootstrap. The `collett` data are not redistributed with this package.
- Archived simulation scripts are retained for reproducibility; diagnostic variants in those scripts are not additional publicly exported estimators.

Reference:

Li, Z. Semiparametric Inference for Marginal Win Odds with Censored Time-to-Event Outcomes. Manuscript under review.
