# winoddsEIF (research prototype)

An initial R package adapted from the user's **scalar** and **semi-competing-risk** simulation programs for the manuscript *Semiparametric Inference for Marginal Win Odds with Censored Time-to-Event Outcomes*.

**Status:** Development version `0.1.0`, unchanged until finalization.
RStudio installation and functional examples succeeded. The scalar
and semi-competing-risk estimators passed numerical equivalence
comparisons with their original simulations in two scenarios each.
The *separate tamoxifen spline analysis* and formal `R CMD check` have
not yet been validated; do not release until those checks are complete.

## Install locally

From the parent directory containing `winoddsEIF/`:

```r
install.packages("winoddsEIF", repos = NULL, type = "source")
library(winoddsEIF)
```

Or use the accompanying `.tar.gz` source archive (downloaded separately):

```r
install.packages("winoddsEIF_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Scalar survival outcomes

```r
set.seed(3)
n <- 240
A <- rbinom(n, 1, .5)
X <- runif(n, -1, 1)
T <- rexp(n, .1*exp(-.4*A + .7*X))
C <- rexp(n, .09*exp(.5*X))
d <- data.frame(A, X, time=pmin(T,C,5),
                status=as.integer(T <= C & T <= 5))
fit <- winprob_scalar(d, tau=5, treatment_prob=.5, folds=5, seed=123)
print(fit)
fit$diagnostics
```

`status=1` indicates that the event was observed at/before the analysis horizon; `status=0` indicates that the event was not observed. A subject with `status=0` but observed follow-up through `tau` contributes a fully revealed outcome at `tau`.

## Prioritized semi-competing risks

Provide one row per participant with six columns:

| Column | Meaning |
|---|---|
| `A` | 0 = control, 1 = treatment |
| `X` | Numeric baseline covariate (one covariate supported for now) |
| `time1` | min(terminal event, censoring, follow-up end) |
| `status1` | 1 when terminal event observed, 0 otherwise |
| `time2` | min(nonterminal event, terminal event, censoring, follow-up end) |
| `status2` | 1 when nonterminal event observed **before** terminal event, 0 otherwise |

**Reproducible example** (synthetic trial data; `m_q=16` for a short demo,
use `m_q=256` to match the numerical integration setting in the simulations):

```r
set.seed(2027)
n <- 240L
tau <- 5
A <- rbinom(n, 1, 0.5)
X <- runif(n, -1, 1)
# Terminal event has higher priority than nonterminal event.
T1 <- rexp(n, rate = 0.10 * exp(-0.35*A + 0.7*X))
T2 <- rexp(n, rate = 0.20 * exp(-0.30*A + 0.4*X))
C  <- rexp(n, rate = 0.08 * exp(0.2*X))
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
fit$diagnostics
fit$psi
fit$win_odds
fit$ci_win_odds
```

In this version, the primary (terminal) event has priority over the secondary (nonterminal) event, and the same fixed horizon `tau` is used for both. The observed follow-up columns can extend beyond `tau` and are administratively truncated by the function.

## Research scope and limitations

- Exactly **one numeric baseline covariate**, binary treatment, **constant known treatment assignment probability**. This is a narrower API than the paper's general semiparametric model.
- Fixed ridge penalty `kappa=2` to match the supplied simulations. No user-selectable tuning yet.
- Model-based predictable residual projections, as in the supplied scripts. Claims of nonparametric efficiency require the assumptions and nuisance-rate conditions in the manuscript; they do not follow merely from use of this package.
- `se_ipcw` is a **first-order model-based approximation**; consult the paper's numerical checks. A failure to compute it is reported in diagnostics.
- No bootstrap. No separate dataset is bundled.
- No q=0 method is exposed publicly; the original unmodified simulation scripts are included solely as archival reproduction material.

## Reproduce the manuscript's simulations

The exact supplied scripts are preserved under `inst/reproducibility/`:

- `scalar_ridge_all_scenarios.R`
- `semicomp_ridge_all_scenarios.R`

**Warning:** The scripts contain top-level commands that **automatically execute many replications** and write CSV/RDS results when sourced. Review their `SIM_N`, `SIM_B`, `RUN_SCENARIOS`, and (for scalar) `RUN_MAIN` settings before running. They are not invoked by package installation or function calls.

## Reproducing the tamoxifen real-data analysis

The package bundles **both** author-supplied scripts in `inst/analysis/`:

- **Main**: `tamoxifen_penalized_spline_fullq.R`, with penalized cubic
  B-spline baseline log-hazards, five-fold cross-fitting, and full q.
- **Supporting**: `tamoxifen_exponential_fullq.R`, with exponential
  working hazards.

The main spline script is **not** the same estimator implementation as
the simple one-covariate public `winprob_semicomp()` simulation interface:
it includes age, log tumor size, and hormone receptor status and spline
baseline hazards. It is packaged as reproducible standalone analysis
code, not presented as a new general-purpose public estimator.

To reproduce the 4- and 5-year spline results:

```r
install.packages("collett")  # one time, if not installed
library(winoddsEIF)
f <- system.file("analysis", "tamoxifen_penalized_spline_fullq.R",
                 package = "winoddsEIF", mustWork = TRUE)
e <- new.env(parent = globalenv())
sys.source(f, envir = e)
print(e$spline_results$table)
```

**Warning:** sourcing either standalone script *executes* its analysis
and creates CSV/RDS files in `getwd()`. Neither runs automatically on
package installation or when calling `library(winoddsEIF)`. See
`inst/analysis/README.md` for usage and methodological distinctions.

Both successful numerical equivalence scripts are retained under
`inst/validation/`.

## Validation checklist before GitHub release or journal submission

1. Run `R CMD build winoddsEIF` and `R CMD check --as-cran winoddsEIF_0.1.0.tar.gz` on a machine with R installed.
2. Run `tests/smoke.R` and compare the one-step estimates against the original programs using identical generated data and folds.
3. Run the scalar and semi-competing replication studies with their published configurations, checking all reported tables, numerical failure counts, and random seeds.
4. Run the included **tamoxifen spline analysis** and compare its numeric results with the manuscript; record the `collett` version and `sessionInfo()`.
5. Add license review, a release tag, session information, and a DOI via Zenodo when ready.

## Acknowledgment of provenance

The internal numerical functions are adapted from the two user-supplied simulation scripts without claiming new independent validation of their numerical implementation.
