# Lightweight smoke tests. Requires an installed winoddsEIF package.
library(winoddsEIF)
set.seed(991)
n <- 220L
A <- rbinom(n, 1, .5); X <- runif(n, -1, 1)
T <- rexp(n, .12 * exp(-.35*A + .7*X))
C <- rexp(n, .09 * exp(.3*X))
scalar <- data.frame(A=A, X=X, time=pmin(T,C,5),
                     status=as.integer(T<=C & T<=5))
fit <- winprob_scalar(scalar, 5, seed=7)
stopifnot(inherits(fit,"winoddsEIF"), is.finite(fit$psi),
          is.finite(fit$se_psi), is.finite(fit$ipcw),
          length(fit$ci_win_odds)==2L)
cat("Scalar smoke test passed.\n")

# Minimal simulated semi-competing histories (terminal and secondary hazards
# here are independent, which the method allows as a special case).
set.seed(992)
A <- rbinom(n, 1, .5); X <- runif(n, -1, 1)
T1 <- rexp(n, .1*exp(-.35*A + .7*X))
T2 <- rexp(n, .2*exp(-.3*A + .4*X))
C  <- rexp(n, .08*exp(.2*X))
semi <- data.frame(A=A, X=X,
                   time1=pmin(T1,C,5), status1=as.integer(T1<=C & T1<=5),
                   time2=pmin(T2,T1,C,5),
                   status2=as.integer(T2<T1 & T2<=C & T2<=5))
fit2 <- winprob_semicomp(semi, 5, seed=8, m_q=16L)
stopifnot(inherits(fit2,"winoddsEIF"), is.finite(fit2$psi),
          is.finite(fit2$se_psi), is.finite(fit2$ipcw))
cat("Semi-competing smoke test passed.\n")
