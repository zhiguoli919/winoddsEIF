\name{winprob_scalar}
\alias{winprob_scalar}
\title{Cross-fitted marginal win probability for scalar survival outcomes}
\usage{winprob_scalar(data, tau, treatment_prob = 0.5, folds = 5L, seed = 2027L)}
\arguments{
  \item{data}{Data frame with \code{A} (0/1 treatment), numeric \code{X},
  \code{time} (minimum of event, censoring, and administrative follow-up),
  and \code{status} (1 if event occurred, 0 otherwise).}
  \item{tau}{Positive restriction horizon.}
  \item{treatment_prob}{Known, constant probability of assignment to treatment.}
  \item{folds}{Number of cross-fitting folds.}
  \item{seed}{Integer seed for fold assignment.}
}
\description{Estimates marginal win probability using coherent IPCW
kernel projections and a cross-fitted efficient-influence-function correction.
The current implementation uses one numeric covariate and exponential
hazards, with arm-specific event models and a pooled censoring model.
It uses a fixed ridge penalty of 2 matching the study simulation code.}
\details{The observed \code{time} may exceed \code{tau}; values past \code{tau}
are administratively truncated in the analysis. A status-0 subject observed
through \code{tau} has a fully revealed restricted outcome. This method
assumes randomized assignment with known constant allocation probability.
The standard error uses empirical cross-fitted influence-function values.
The optional IPCW comparator has a separate approximate variance calculation.}
\value{An object of class \code{winoddsEIF} with \code{psi}, \code{se_psi},
\code{ci_psi}, \code{win_odds}, \code{ci_win_odds}, \code{pvalue_win_odds},
\code{plugin}, \code{ipcw}, \code{se_ipcw}, and \code{diagnostics}.}
\examples{\dontrun{
set.seed(12)
n <- 200
A <- rbinom(n, 1, .5); X <- runif(n, -1, 1)
T <- rexp(n, .1 * exp(-.3*A + .7*X))
C <- rexp(n, .09 * exp(.4*X))
d <- data.frame(A=A, X=X, time=pmin(T,C,5), status=as.integer(T<=C & T<=5))
fit <- winprob_scalar(d, tau=5)
print(fit)
}}
