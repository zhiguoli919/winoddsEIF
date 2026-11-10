\name{winprob_semicomp}
\alias{winprob_semicomp}
\title{Cross-fitted marginal win probability for semi-competing risks}
\usage{winprob_semicomp(data, tau, treatment_prob = 0.5, folds = 5L,
    seed = 2027L, m_q = 256L, mu_method = c("model", "quadratic"),
    g_floor = 1e-8)}
\arguments{
  \item{data}{Data frame containing \code{A}, \code{X}, \code{time1},
    \code{status1}, \code{time2}, \code{status2}. The primary terminal-event
    observation \code{time1} is the minimum of terminal event, censoring,
    and administrative follow-up; \code{status1} indicates terminal event.
    The secondary observation \code{time2} is the minimum of nonterminal
    event, terminal event, censoring, and administrative follow-up;
    \code{status2} indicates a nonterminal event observed strictly before
    terminal event or censoring.}
  \item{tau}{Positive restriction horizon.}
  \item{treatment_prob}{Known constant probability of assignment to treatment.}
  \item{folds}{Number of cross-fitting folds.}
  \item{seed}{Integer seed for fold assignment.}
  \item{m_q}{Number of quadrature grid points for the predictable projection.}
  \item{mu_method}{Baseline regression method: \code{model} or \code{quadratic}.}
  \item{g_floor}{Fitted censoring survival values below this numerical floor cause a fitting error; weights are not truncated.}
}
\description{Estimates marginal win probability and win odds for a prioritized
terminal/nonterminal outcome with event-history-dependent censoring using
cross-fitted one-step inference. The terminal event takes priority.}
\details{The first version supports one continuous numeric baseline covariate,
known constant treatment assignment probability, and piecewise-exponential
history-dependent event and censoring models. The ridge penalty is fixed at 2
for consistency with the supplied simulations.}
\value{An object of class \code{winoddsEIF}.}
\examples{\dontrun{
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
}}
