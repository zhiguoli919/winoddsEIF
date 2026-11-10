library(winoddsEIF)

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
