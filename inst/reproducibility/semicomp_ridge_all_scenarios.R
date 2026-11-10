# Prioritized semi-competing risks: simulation and estimator comparison
# --------------------------------------------------------------------
# Primary outcome: terminal event T1; secondary outcome: nonterminal T2. The terminal and censoring hazards may change after the nonterminal event.
# Estimators: cross-fitted one-step (FULL_EIF), normalized IPCW, plug-in, and an oracle-hazard benchmark. Ridge-regularized piecewise-exponential
# nuisance models use fixed kappa=2; the intercept is unpenalized. Reference win probability is estimated by an independent Monte Carlo run;
# the reported reference Monte Carlo SE should be considered in bias checks.
#
# Default n=800, B=5000 per scenario, M_q=256. 

# Outputs are saved to getwd(). 
# --------------------------------------------------------------------

SIM_N <- 200L
SIM_B <- 5000L            # Monte Carlo replications PER scenario
SIM_V <- 5L
SIM_SEED <- 2031L
RIDGE_KAPPA <- 2.0
SIM_M_Q <- 256L            # Conditional-projection integration grid (16 x 16)
SIM_M_REF <- 500000L       # Independent reference-truth sample size
SIM_ORACLE <- TRUE
RUN_SCENARIOS <- c("moderate","stronger","strongest")

rX <- function(n) qnorm(runif(n, pnorm(-2), pnorm(2)))

semi_parameters <- function(scenario = "moderate") {
  if (scenario == "moderate") {
    ## Moderate heterogeneity and censoring
    return(list(tau=5, lambda10=.10, beta1A=-.35, beta1X=1.25,
      beta1N2=1.25, lambda20=.20, beta2A=-.30, beta2X=1.25,
      piA=.5, gammaX=.80, gammaN2=1.25, target_censoring=.30))
  }
  if (scenario == "stronger") {
    ## Censoring DGP: gammaX=0.50, gammaN2=0.50
    ## Event hazards unchanged. Calibrate c0 to 40% censoring.
    ## The conservative true G lower bound is about 0.024.
    return(list(tau=5, lambda10=.10, beta1A=-.35, beta1X=1.25,
      beta1N2=1.60, lambda20=.20, beta2A=-.30, beta2X=1.25,
      piA=.5, gammaX=.50, gammaN2=0.50, target_censoring=.40))
  }
  if (scenario == "strongest") {
    ## Censoring DGP: gammaX=0.40, gammaN2=0.40
    ## Event hazards unchanged. Calibrate c0 to 50% censoring.
    ## The conservative true G lower bound is about 0.015.
    return(list(tau=5, lambda10=.10, beta1A=-.35, beta1X=1.25,
      beta1N2=2.00, lambda20=.20, beta2A=-.30, beta2X=1.25,
      piA=.5, gammaX=.40, gammaN2=0.40, target_censoring=.50))
  }
  stop("Unknown scenario. Use 'moderate', 'stronger', or 'strongest'.")
}

semi_h <- function(p1,s1,p0,s0,tol=1e-10) {
  out <- numeric(length(p1))
  gt <- p1 > p0+tol
  tiep <- abs(p1-p0) <= tol
  out[gt] <- 1
  out[tiep & (s1 > s0+tol)] <- 1
  out[tiep & (abs(s1-s0)<=tol)] <- .5
  out
}

semi_rate_T2 <- function(a,x,par)
  par$lambda20 * exp(par$beta2A*a + par$beta2X*x)
semi_rate_T1_pre <- function(a,x,par)
  par$lambda10 * exp(par$beta1A*a + par$beta1X*x)

semi_T1_given_T2 <- function(a,x,T2,par) {
  lam0 <- semi_rate_T1_pre(a,x,par)
  lam1 <- lam0 * exp(par$beta1N2)
  E <- rexp(length(x))
  before <- E <= lam0*T2
  ifelse(before, E/lam0, T2+(E-lam0*T2)/lam1)
}

semi_generate <- function(n,par,c0) {
  X <- rX(n)
  A <- rbinom(n,1,par$piA)
  T2 <- rexp(n,rate=semi_rate_T2(A,X,par))
  T1 <- semi_T1_given_T2(A,X,T2,par)
  baseC <- c0*exp(par$gammaX*X)
  E <- rexp(n)
  before <- E <= baseC*T2
  C <- ifelse(before,E/baseC,
              T2+(E-baseC*T2)/(baseC*exp(par$gammaN2)))
  K <- pmin(T1,par$tau)
  R <- as.integer(C>=K)
  tstop <- pmin(K,C)
  d1 <- as.integer(T1 <= C & T1 <= par$tau)
  d2 <- as.integer(T2 <= T1 & T2 <= C & T2 <= par$tau)
  data.frame(X=X,A=A,obs_time=tstop,R=R,cens_event=1L-R,
    T1tilde=pmin(T1,C,par$tau),Delta1=d1,
    T2tilde=pmin(T2,T1,C,par$tau),Delta2=d2,
    Yp=ifelse(R==1L,K,NA_real_),
    Ys=ifelse(R==1L,pmin(T2,T1,par$tau),NA_real_))
  ## Restricted outcomes for censored subjects are masked as NA;
  ## the estimator cannot accidentally access latent post-censoring data.
}

semi_reference_truth <- function(par, M=500000, seed=20261102) {
  set.seed(seed)
  x1 <- rX(M); x0 <- rX(M)
  t21 <- rexp(M,rate=semi_rate_T2(1,x1,par))
  t20 <- rexp(M,rate=semi_rate_T2(0,x0,par))
  t11 <- semi_T1_given_T2(1,x1,t21,par)
  t10 <- semi_T1_given_T2(0,x0,t20,par)
  pair_score <- semi_h(pmin(t11,par$tau),pmin(t21,t11,par$tau),
                       pmin(t10,par$tau),pmin(t20,t10,par$tau))
  ans <- mean(pair_score)
  attr(ans,"mcse") <- sd(pair_score)/sqrt(length(pair_score))
  ans
}

semi_calibrate_censoring <- function(par,M=100000,seed=20261103) {
  set.seed(seed)
  X <- rX(M); A <- rbinom(M,1,par$piA)
  T2 <- rexp(M,rate=semi_rate_T2(A,X,par))
  T1 <- semi_T1_given_T2(A,X,T2,par)
  K <- pmin(T1,par$tau)
  exposure <- exp(par$gammaX*X) *
    (pmin(K,T2)+exp(par$gammaN2)*pmax(K-T2,0))
  uniroot(function(c0) mean(-expm1(-c0*exposure)) -
                       par$target_censoring,
          interval=c(1e-8,20))$root
}

# Conservative lower bound of true G(t | X, N2 history) over
# X in [-2, 2] and t in [0, tau], assuming gammaX,gammaN2 >= 0.
# It allows an arbitrarily early nonterminal event (worst-case history).
semi_true_min_G <- function(par,c0) {
  stopifnot(par$gammaX >= 0, par$gammaN2 >= 0)
  exp(-c0 * par$tau * exp(2*par$gammaX + par$gammaN2))
}

## Piecewise exponential likelihood via Poisson likelihood.
## All event-history covariates are built from OBSERVED data.
semi_ped <- function(d, type) {
  ans <- vector("list",2*nrow(d)); k<-0L
  for (i in seq_len(nrow(d))) {
    end <- if (type=="C") d$obs_time[i] else d$T1tilde[i]
    event <- if (type=="C") d$cens_event[i] else d$Delta1[i]
    t2_obs <- if (d$Delta2[i]==1) d$T2tilde[i] else Inf
    t0 <- min(end,t2_obs)
    if (t0>0) {
      k<-k+1L
      ans[[k]]<-data.frame(exposure=t0,event=as.integer(event==1 && end<=t2_obs),
                            X=d$X[i],N2=0)
    }
    t1 <- end-t2_obs
    if (is.finite(t2_obs) && t1>0) {
      k<-k+1L
      ans[[k]]<-data.frame(exposure=t1,event=as.integer(event==1 && end>t2_obs),
                            X=d$X[i],N2=1)
    }
  }
  if (k==0L) return(data.frame(exposure=numeric(),event=numeric(),
                                 X=numeric(),N2=numeric()))
  do.call(rbind,ans[seq_len(k)])
}

## Penalized piecewise-exponential likelihood (Poisson PED),
## log hazard = alpha + bX*X + bN2*N2.
## Maximize the TOTAL log likelihood - kappa*(bX^2+bN2^2)/2.
## Profile out the intercept to reduce numerical optimization.
semi_fit_rate <- function(d,useN2=TRUE,kappa=RIDGE_KAPPA) {
  # Profiled ridge likelihood is strictly convex in its slopes when
  # kappa > 0. Use analytic-Hessian Newton steps with backtracking,
  # rather than BFGS plus an overly strict post-hoc gradient check.
  if (nrow(d)<6L || any(!is.finite(d$exposure)) ||
      any(d$exposure<=0) || !is.finite(kappa) || kappa<=0 ||
      any(!is.finite(d$event)) || any(!is.finite(d$X)) ||
      (useN2 && any(!is.finite(d$N2))))
    stop("Invalid at-risk data / ridge penalty")
  E <- sum(d$event)
  if (E < 1L) {
    # No observed events: report the existing transparent fallback.
    rate <- max(E/sum(d$exposure),1e-12)
    return(list(alpha=log(rate),bX=0,bN2=0,fallback=TRUE))
  }
  Z <- if (useN2) cbind(d$X,d$N2) else matrix(d$X,ncol=1L)
  storage.mode(Z) <- "double"
  p <- ncol(Z)
  eventZ <- as.numeric(crossprod(Z,d$event))
  log_exposure <- log(d$exposure)
  score_tol <- 1e-8 * max(1,E)

  evaluate <- function(beta, need_hessian=TRUE) {
    eta <- as.numeric(Z %*% beta)+log_exposure
    shift <- max(eta)
    ew <- exp(eta-shift)
    prob <- ew/sum(ew)
    logtotal <- shift+log(sum(ew))
    zmean <- as.numeric(crossprod(Z,prob))
    objective <- E*logtotal - sum(eventZ*beta) +
      .5*kappa*sum(beta^2)
    gradient <- E*zmean-eventZ+kappa*beta
    if (!need_hessian)
      return(list(value=objective,gradient=gradient,logtotal=logtotal))
    Zc <- sweep(Z,2L,zmean,FUN="-")
    H <- E*crossprod(Zc, Zc*as.numeric(prob)) +
      diag(rep(kappa,p),nrow=p)
    list(value=objective,gradient=gradient,hessian=H,
         logtotal=logtotal)
  }

  beta <- rep(0,p)
  converged <- FALSE
  for (iter in seq_len(100L)) {
    cur <- evaluate(beta)
    if (max(abs(cur$gradient)) <= score_tol) {
      converged <- TRUE
      break
    }
    # Positive-definite Hessian because kappa > 0.
    step <- as.numeric(solve(cur$hessian,cur$gradient))
    directional <- sum(cur$gradient*step)
    if (!is.finite(directional) || directional<=0)
      stop("Invalid ridge Newton direction")
    scale <- 1
    accepted <- FALSE
    for (backtrack in seq_len(50L)) {
      candidate <- beta-scale*step
      next_fit <- evaluate(candidate,need_hessian=FALSE)
      if (is.finite(next_fit$value) &&
          next_fit$value <= cur$value-1e-4*scale*directional) {
        accepted <- TRUE
        break
      }
      scale <- scale/2
    }
    if (!accepted) {
      # The objective can be indistinguishable at machine precision
      # even when the absolute gradient is a little above the target.
      if (max(abs(cur$gradient)) < 1e-6*max(1,E)) {
        converged <- TRUE
        break
      }
      stop("Ridge Newton line search failed")
    }
    beta <- candidate
  }
  final <- evaluate(beta,need_hessian=FALSE)
  if (!converged && max(abs(final$gradient)) > 1e-6*max(1,E))
    stop("Ridge Newton solver did not converge: residual score = ",
         signif(max(abs(final$gradient)),4))
  alpha <- log(E)-final$logtotal
  if (!is.finite(alpha) || any(!is.finite(beta)))
    stop("Non-finite ridge hazard coefficients")
  list(alpha=as.numeric(alpha),bX=as.numeric(beta[1L]),
       bN2=if(useN2) as.numeric(beta[2L]) else 0,
       fallback=FALSE)
}

semi_predict_rate <- function(fit,x,n2=0) {
  eta <- fit$alpha + fit$bX*x + fit$bN2*n2
  if (any(!is.finite(eta)) || any(eta>50) || any(eta< -700))
    stop("Extreme fitted log hazard; inspect this replication")
  pmax(exp(eta),1e-12)
}

semi_fit_models <- function(d, oracle=FALSE, par=NULL,c0=NULL) {
  if (oracle) {
    if (is.null(par) || is.null(c0)) stop("Oracle needs par,c0")
    return(lapply(0:1,function(a) list(
      C=list(alpha=log(c0),bX=par$gammaX,bN2=par$gammaN2,
             fallback=FALSE),
      T1=list(alpha=log(par$lambda10)+par$beta1A*a,
              bX=par$beta1X,bN2=par$beta1N2,fallback=FALSE),
      T2=list(alpha=log(par$lambda20)+par$beta2A*a,
              bX=par$beta2X,bN2=0,fallback=FALSE)
    )))
  }
  ans <- vector("list",2)
  for (a in 0:1) {
    da <- d[d$A==a,,drop=FALSE]
    pedC <- semi_ped(da,"C")
    pedT1 <- semi_ped(da,"T1")
    pedT2 <- data.frame(
      exposure=pmax(da$T2tilde,1e-12),event=da$Delta2,X=da$X)
    ans[[a+1L]] <- list(
      C=semi_fit_rate(pedC),T1=semi_fit_rate(pedT1),
      T2=semi_fit_rate(pedT2,FALSE))
  }
  ans
}

semi_observed_t2 <- function(d) ifelse(d$Delta2==1,d$T2tilde,Inf)
semi_G <- function(t,x,t2obs,fitC) {
  lam0<-semi_predict_rate(fitC,x,0)
  lam1<-semi_predict_rate(fitC,x,1)
  exp(-lam0*pmin(t,t2obs)-lam1*pmax(t-t2obs,0))
}

## Fast weighted cdf; weights can total <1 for the T1=tau subset.
semi_wcdf <- function(z,w,tol=1e-10) {
  ord<-order(z); z<-z[ord]; w<-w[ord]
  sums<-c(0,cumsum(w)); total<-sum(w)
  less<-function(y) sums[findInterval(y-tol,z)+1L]
  leq <-function(y) sums[findInterval(y+tol,z)+1L]
  list(less=less,leq=leq,total=total)
}

## Training-set normalized IPCW marginal distributions, coherent h and psi.
semi_make_projection <- function(train,models,tau,g_floor=1e-8) {
  dlist<-vector("list",2)
  for (a in 0:1) {
    z<-train[train$A==a & train$R==1,,drop=FALSE]
    if (nrow(z)==0) stop("No complete outcomes in one arm")
    t2obs<-semi_observed_t2(z)
    G<-semi_G(z$obs_time,z$X,t2obs,models[[a+1]]$C)
    if (any(!is.finite(G)) || any(G<=g_floor))
      stop("Fitted censoring survival below numerical floor")
    w<-1/G
    dlist[[a+1]]<-list(Yp=z$Yp,Ys=z$Ys,X=z$X,w=w/sum(w),G=G)
  }
  r1<-dlist[[2]]; r0<-dlist[[1]]
  F1p<-semi_wcdf(r1$Yp,r1$w)
  F0p<-semi_wcdf(r0$Yp,r0$w)
  i1<-which(abs(r1$Yp-tau)<1e-10)
  i0<-which(abs(r0$Yp-tau)<1e-10)
  F1s<-semi_wcdf(r1$Ys[i1],r1$w[i1])
  F0s<-semi_wcdf(r0$Ys[i0],r0$w[i0])
  h1<-function(yp,ys) {
    out<-numeric(length(yp)); early<-yp<tau-1e-10
    if (any(early)) out[early]<-F0p$less(yp[early])
    if (any(!early)) {
      ss<-ys[!early]
      out[!early]<-F0p$less(tau) + F0s$less(ss) +
        .5*(F0s$leq(ss)-F0s$less(ss))
    }
    pmin(pmax(out,0),1)
  }
  h0<-function(yp,ys) {
    out<-numeric(length(yp)); early<-yp<tau-1e-10
    if (any(early)) out[early]<-1-F1p$leq(yp[early])
    if (any(!early)) {
      ss<-ys[!early]
      out[!early]<-F1s$total-F1s$leq(ss)+
        .5*(F1s$leq(ss)-F1s$less(ss))
    }
    pmin(pmax(out,0),1)
  }
  psi<-sum(r1$w*h1(r1$Yp,r1$Ys))
  list(psi=psi,h1=h1,h0=h0)
}

semi_fit_mu <- function(train,models,proj,g_floor=1e-8) {
  ans<-vector("list",2)
  for (a in 0:1) {
    z<-train[train$A==a & train$R==1,,drop=FALSE]
    if (nrow(z)<6) {ans[[a+1]]<-list(fit=NULL,fallback=.5);next}
    h<-if (a==1) proj$h1(z$Yp,z$Ys) else proj$h0(z$Yp,z$Ys)
    g<-pmax(semi_G(z$obs_time,z$X,semi_observed_t2(z),
                   models[[a+1]]$C),g_floor)
    fit<-suppressWarnings(tryCatch(lm(h~X+I(X^2),
                  weights=1/g,data=data.frame(h=h,X=z$X)),
                  error=function(e) NULL))
    ans[[a+1]]<-list(fit=fit,fallback=weighted.mean(h,1/g))
  }
  ans
}
semi_mu <- function(f,x) {
  out<-if (is.null(f$fit)) rep(f$fallback,length(x))
       else suppressWarnings(as.numeric(predict(f$fit,
                                   newdata=data.frame(X=x))))
  if (any(!is.finite(out))) out[!is.finite(out)]<-f$fallback
  pmin(pmax(out,0),1)
}

## Deterministic conditional-quantile evaluation of q(u,H_{u-}).
## If secondary event has occurred, its EXACT observed time enters Ys.
## Before T2, simulate future T2 and T1 (with T1 hazard switching at T2).
semi_q <- function(a,u,x,t2obs,fit,proj,mu,tau,M_q=100) {
  if (u>tau) return(0)
  L<-max(4L,as.integer(ceiling(sqrt(M_q))))
  U<-(seq_len(L)-.5)/L
  rate2<-semi_predict_rate(fit$T2,x)
  rate10<-semi_predict_rate(fit$T1,x,0)
  rate11<-semi_predict_rate(fit$T1,x,1)
  if (is.finite(t2obs) && t2obs < u+1e-10) {
    e1<--log1p(-U)/rate11
    yp<-pmin(u+e1,tau)
    ys<-rep(pmin(t2obs,tau),L)
  } else {
    e2<-rep(-log1p(-U)/rate2,each=L)
    E1<-rep(-log1p(-U),times=L)
    t2<-u+e2
    before<-E1<=rate10*e2
    t1<-ifelse(before,u+E1/rate10,
               t2+(E1-rate10*e2)/rate11)
    yp<-pmin(t1,tau)
    ys<-pmin(t2,t1,tau)
  }
  h<-if (a==1) proj$h1(yp,ys) else proj$h0(yp,ys)
  mean(h)-mu
}

## Eight-point Gauss-Legendre integration, separately before/after T2.
semi_GLx<-c(-.9602898564975363,-.7966664774136267,
 -.5255324099163290,-.1834346424956498,.1834346424956498,
 .5255324099163290,.7966664774136267,.9602898564975363)
semi_GLw<-c(.1012285362903763,.2223810344533745,
 .3137066458778873,.3626837833783620,.3626837833783620,
 .3137066458778873,.2223810344533745,.1012285362903763)

semi_augment <- function(a,x,stop,cens,t2obs,fit,proj,mu,tau,
                         M_q=100,g_floor=1e-8) {
  if (stop<=0) return(0)
  q_eval<-function(u) semi_q(a,u,x,t2obs,fit,proj,mu,tau,M_q)
  G_eval<-function(u) {
    G<-semi_G(u,x,t2obs,fit$C)
    if (any(!is.finite(G)) || any(G<=g_floor))
      stop("Fitted censoring survival below numerical floor")
    G
  }
  integrate_part<-function(left,right,n2) {
    if (right<=left) return(0)
    u<-(left+right)/2+(right-left)/2*semi_GLx
    q<-vapply(u,q_eval,numeric(1))
    lamC<-semi_predict_rate(fit$C,x,n2)
    (right-left)/2*sum(semi_GLw*(q/G_eval(u))*lamC)
  }
  cut<-min(stop,t2obs)
  drift<-integrate_part(0,cut,0)
  if (is.finite(t2obs) && t2obs<stop)
    drift<-drift+integrate_part(t2obs,stop,1)
  jump<-if (cens==1) q_eval(stop)/G_eval(stop) else 0
  as.numeric(jump-drift)
}

## Cross-fitted IPCW comparator using ONLY out-of-fold fitted weights.
semi_ipcw_crossfit <- function(dat,w) {
  i1 <- which(dat$A==1L & dat$R==1L & is.finite(w) & w>0)
  i0 <- which(dat$A==0L & dat$R==1L & is.finite(w) & w>0)
  if (!length(i1) || !length(i0)) stop("No revealed outcomes for IPCW")
  w1 <- w[i1]/sum(w[i1]); w0 <- w[i0]/sum(w[i0])
  sum(vapply(seq_along(i1),function(j) {
    k <- i1[j]
    w1[j]*sum(w0*semi_h(rep(dat$Yp[k],length(i0)),
                            rep(dat$Ys[k],length(i0)),
                            dat$Yp[i0],dat$Ys[i0]))
  },numeric(1)))
}


## Censoring-model score/information for the parametric IPCW variance.
## Within each arm, log lambda_C(u) = b0 + bX X + bN2 I(T2 < u).
## Values below are evaluated on held-out subjects at their OUT-OF-FOLD
## fitted censoring coefficient vector, not at the simulation truth.
semi_censor_score_information <- function(d, fitC) {
  if (isTRUE(fitC$fallback))
    stop("Censoring fit used fallback; parametric IPCW SE unavailable")
  eta0 <- fitC$alpha + fitC$bX*d$X
  eta1 <- eta0 + fitC$bN2
  if (any(!is.finite(c(eta0,eta1))) ||
      any(c(eta0,eta1) > 50 | c(eta0,eta1) < -700))
    stop("Fitted censoring hazard reached a numerical bound; IPCW SE unavailable")
  lc0 <- exp(eta0)
  lc1 <- exp(eta1)
  t2 <- semi_observed_t2(d)
  u0 <- pmin(d$obs_time,t2)
  u1 <- pmax(d$obs_time-t2,0)
  H0 <- lc0*u0
  H1 <- lc1*u1
  H <- H0+H1
  z_jump <- as.numeric(is.finite(t2) & t2 < d$obs_time)
  S <- cbind(d$cens_event-H,
             d$X*(d$cens_event-H),
             d$cens_event*z_jump-H1)
  colnames(S) <- c("intercept","X","N2")
  dH <- cbind(H,d$X*H,H1)
  colnames(dH) <- colnames(S)
  J <- array(0,dim=c(3L,3L,nrow(d)))
  for (i in seq_len(nrow(d))) {
    z0 <- c(1,d$X[i],0)
    z1 <- c(1,d$X[i],1)
    J[,,i] <- H0[i]*tcrossprod(z0) + H1[i]*tcrossprod(z1)
  }
  list(score=S,dH=dH,info=J)
}

## Complete-case normalized IPCW U-statistic and its direct IF SE.
## This implementation accounts for estimation of the piecewise-exponential
## censoring hazards in each arm. It assumes KNOWN treatment probabilities.
##
## Weighted IF:    sum_a W_ai / mean(W_a) * [h_a(Y_i) - psi].
## Hazard IF:      sum_a c_fold(i) d_a' J_a^{-1} S_{C,a,i},
## where d_a = E[W_a (h_a-psi) grad_b Lambda_C(K)] / E(W_a),
##       J_a = E[I(A=a) int Y_C lambda_C z z' du], and
##       c_fold(i) = sum_{v != fold(i)} n_v/(n-n_v).
## c_fold equals 1 with equal-sized validation folds.
## Cross-fitted normalized IPCW + direct penalized-score IF SE.
## Fold-specific derivative D_{a,v} is evaluated on held-out
## complete subjects, with out-of-fold normalized IPCW weights.
## The nuisance perturbation for training subjects uses
## (sum_training info_i + diag(0,kappa,kappa))^{-1} score_i.
## This accounts for overlap of training folds and ridge curvature.
semi_ipcw_inference <- function(dat,w,folds,cfits,kappa=RIDGE_KAPPA) {
  n <- nrow(dat); V <- max(folds)
  i1 <- which(dat$A==1L & dat$R==1L & is.finite(w) & w>0)
  i0 <- which(dat$A==0L & dat$R==1L & is.finite(w) & w>0)
  if (!length(i1) || !length(i0)) stop("No complete cases in an arm")
  w1 <- w[i1]/sum(w[i1]); w0 <- w[i0]/sum(w[i0])
  tau <- attr(dat,"tau")
  F1p <- semi_wcdf(dat$Yp[i1],w1)
  F0p <- semi_wcdf(dat$Yp[i0],w0)
  t1 <- abs(dat$Yp[i1]-tau)<1e-10
  t0 <- abs(dat$Yp[i0]-tau)<1e-10
  F1s <- semi_wcdf(dat$Ys[i1][t1],w1[t1])
  F0s <- semi_wcdf(dat$Ys[i0][t0],w0[t0])
  h1 <- function(yp,ys) {
    ans <- numeric(length(yp)); early <- yp<tau-1e-10
    if (any(early)) ans[early] <- F0p$less(yp[early])
    if (any(!early)) {
      ss <- ys[!early]
      ans[!early] <- F0p$less(tau)+F0s$less(ss)+
        .5*(F0s$leq(ss)-F0s$less(ss))
    }
    ans
  }
  h0 <- function(yp,ys) {
    ans <- numeric(length(yp)); early <- yp<tau-1e-10
    if (any(early)) ans[early] <- 1-F1p$leq(yp[early])
    if (any(!early)) {
      ss <- ys[!early]
      ans[!early] <- F1s$total-F1s$leq(ss)+
        .5*(F1s$leq(ss)-F1s$less(ss))
    }
    ans
  }
  hv1 <- h1(dat$Yp[i1],dat$Ys[i1])
  hv0 <- h0(dat$Yp[i0],dat$Ys[i0])
  h <- numeric(n); h[i1] <- hv1; h[i0] <- hv0
  wn <- numeric(n); wn[i1] <- w1; wn[i0] <- w0
  psi <- sum(w1*hv1)
  if(abs(psi-sum(w0*hv0))>1e-9)
    stop("Kernel projection mismatch in IPCW")
  weighted_IF <- numeric(n)
  ii <- c(i1,i0)
  weighted_IF[ii] <- n*wn[ii]*(h[ii]-psi)
  nuisance_IF <- numeric(n)
  penalty_matrix <- diag(c(0,kappa,kappa),3L)
  for (v in seq_len(V)) {
    for (a in 0:1) {
      tr <- which(folds!=v & dat$A==a)
      va <- which(folds==v & dat$A==a & dat$R==1L)
      fitC <- cfits[[v]][[a+1L]]
      train <- semi_censor_score_information(dat[tr,,drop=FALSE],fitC)
      H <- apply(train$info,c(1,2),sum)+penalty_matrix
      if (any(!is.finite(H)) || rcond(H)<1e-10)
        stop("Near-singular penalized censor information")
      D <- numeric(3)
      if (length(va)) {
        valid <- semi_censor_score_information(dat[va,,drop=FALSE],fitC)
        D <- colSums(valid$dH * as.numeric(wn[va]*(h[va]-psi)))
      }
      coeff <- as.numeric(solve(H,D))
      nuisance_IF[tr] <- nuisance_IF[tr] +
        n*as.numeric(train$score %*% coeff)
    }
  }
  IF <- weighted_IF+nuisance_IF
  if (!all(is.finite(IF))) stop("Nonfinite IPCW IF")
  list(psi=psi,se=sqrt(mean((IF-mean(IF))^2)/n),
       phi=IF,weighted_if=weighted_IF,nuisance_if=nuisance_IF)
}

## Fold-specific full and q=0 one-step; return per-observation IF values.
## Default mu_method="model" integrates h under the fitted event model,
## matching the construction used for q. "quadratic" reproduces the
## working-regression mu from the original standalone scripts.
semi_one_step <- function(dat,V=5,M_q=64,g_floor=1e-8,seed=123,
                          piA=.5,mu_method=c("model","quadratic"),
                          oracle=FALSE, par=NULL,c0=NULL) {
  mu_method <- match.arg(mu_method)
  n <- nrow(dat); tau <- attr(dat,"tau")
  if (is.null(tau)) stop("Provide attr(dat,'tau')")
  if (V<2 || V>n) stop("V must be between 2 and n")
  if (piA<=0 || piA>=1) stop("piA must be in (0,1)")
  set.seed(seed)
  folds <- sample(rep(seq_len(V),length.out=n))
  phi_full <- numeric(n)
  psi_pl <- numeric(V); minG <- numeric(V); w_ipcw <- numeric(n)
  cfits <- vector("list",V)
  for (v in seq_len(V)) {
    train <- dat[folds!=v,,drop=FALSE]
    valid <- dat[folds==v,,drop=FALSE]
    idx <- which(folds==v)
    models <- semi_fit_models(train,oracle=oracle,par=par,c0=c0)
    proj <- semi_make_projection(train,models,tau,g_floor)
    mufit <- if (mu_method=="quadratic")
      semi_fit_mu(train,models,proj,g_floor) else NULL
    psi_pl[v] <- proj$psi
    minG[v] <- 1
    cfits[[v]] <- lapply(models,function(m) m$C)
    for (j in seq_along(idx)) {
      i <- idx[j]; a <- valid$A[j]; x <- valid$X[j]
      if (mu_method=="model") {
        mu1 <- semi_q(1L,0,x,Inf,models[[2]],proj,0,tau,M_q)
        mu0 <- semi_q(0L,0,x,Inf,models[[1]],proj,0,tau,M_q)
      } else {
        mu1 <- semi_mu(mufit[[2]],x)
        mu0 <- semi_mu(mufit[[1]],x)
      }
      mu1 <- pmin(pmax(mu1,0),1)
      mu0 <- pmin(pmax(mu0,0),1)
      mu_a <- if (a==1L) mu1 else mu0
      fit <- models[[a+1L]]
      t2obs <- if (valid$Delta2[j]==1L) valid$T2tilde[j] else Inf
      stop_i <- valid$obs_time[j]
      Gstop <- semi_G(stop_i,x,t2obs,fit$C)
      if (!is.finite(Gstop) || Gstop<=g_floor)
        stop("Fitted censoring survival below numerical floor")
      minG[v] <- min(minG[v],Gstop)
      ipcw <- 0
      if (valid$R[j]==1L) {
        h <- if (a==1L) proj$h1(valid$Yp[j],valid$Ys[j])
             else proj$h0(valid$Yp[j],valid$Ys[j])
        ipcw <- (h-mu_a)/Gstop
        pa <- if (a==1L) piA else 1-piA
        w_ipcw[i] <- 1/(pa*Gstop)
      }
      aug <- semi_augment(a,x,stop_i,valid$cens_event[j],t2obs,
                          fit,proj,mu_a,tau,M_q,g_floor)
      common <- mu1+mu0-2*psi_pl[v]
      w <- 1/(if(a==1L) piA else 1-piA)
      phi_full[i] <- common+w*(ipcw+aug)
    }
  }
  psi_pl_mean <- sum(tabulate(folds,nbins=V)*psi_pl)/n
  psi_full <- psi_pl_mean+mean(phi_full)
  se_full <- sqrt(mean((phi_full-mean(phi_full))^2)/n)
  if (oracle) {
    return(list(full=psi_full,se_full=se_full,plugin=psi_pl_mean,
                minG=min(minG),phi_full=phi_full,folds=folds))
  }
  ipcw_point <- semi_ipcw_crossfit(dat,w_ipcw)
  ipcw_fit <- tryCatch(
    semi_ipcw_inference(dat,w_ipcw,folds,cfits),
    error=function(e) list(se=NA_real_,phi=NULL,
                           error=conditionMessage(e)))
  if (!is.null(ipcw_fit$psi) && abs(ipcw_fit$psi-ipcw_point)>1e-8)
    stop("IPCW point estimator and IF projection disagree")
  list(full=psi_full,se_full=se_full,
       ipcw=ipcw_point,se_ipcw=ipcw_fit$se,
       ipcw_se_error=if(is.null(ipcw_fit$error)) NA_character_ else ipcw_fit$error,
       phi_ipcw=ipcw_fit$phi,
       plugin=psi_pl_mean,
       minG=min(minG),phi_full=phi_full,
       folds=folds,mu_method=mu_method)
}

## Simulation summary records estimation performance and failures.
## counts replications with a finite SE (N_SE).
semi_summarize <- function(x,se=NULL,truth) {
  ok <- is.finite(x); N <- sum(ok)
  z <- x[ok]
  seok <- ok & !is.null(se) & is.finite(se) & se>0
  s <- if (is.null(se)) numeric(0) else se[seok]
  data.frame(N=N,truth=truth,
    bias=if(N) mean(z-truth) else NA_real_,
    empirical_SD=if(N>1) sd(z) else NA_real_,
    RMSE=if(N) sqrt(mean((z-truth)^2)) else NA_real_,
    mean_SE=if(length(s)) mean(s) else NA_real_,
    RMS_SE=if(length(s)) sqrt(mean(s^2)) else NA_real_,
    coverage95=if(length(s)) mean(abs(x[seok]-truth)<=qnorm(.975)*s)
      else NA_real_,
    variance_ratio=if(length(s)>1) mean(s^2)/var(x[seok])
      else NA_real_,
    N_SE=length(s))
}

run_semicomp_simulation <- function(scenario="moderate",B=SIM_B,
    n=SIM_N,V=SIM_V,M_q=SIM_M_Q,M_ref=SIM_M_REF,
    g_floor=1e-8,seed=SIM_SEED,oracle=SIM_ORACLE,
    mu_method="model",verbose=TRUE) {
  par <- semi_parameters(scenario)
  truth_out <- semi_reference_truth(par,M=M_ref)
  truth <- as.numeric(truth_out)
  truth_mcse <- attr(truth_out,"mcse")
  c0 <- semi_calibrate_censoring(par)
  true_min_G <- semi_true_min_G(par,c0)
  if (scenario %in% c("stronger", "strongest") && true_min_G < 0.01)
    stop("DGP fails the prespecified true-G bound: ", true_min_G)
  if(verbose) cat("Semi-competing:",scenario,"n=",n,"B=",B,
                  "kappa=",RIDGE_KAPPA,"M_q=",M_q,
                  "truth=",signif(truth,7),"truth_MCSE=",
                  signif(truth_mcse,3),"c0=",signif(c0,5),
                  "min_true_G=",signif(true_min_G,5),"\n")
  out <- data.frame(replication=seq_len(B),full=NA_real_,
    se_full=NA_real_,
    ipcw=NA_real_,se_ipcw=NA_real_,plugin=NA_real_,
    oracle=NA_real_,se_oracle=NA_real_,
    censor_fraction=NA_real_,min_G=NA_real_,
    ipcw_se_error=rep("",B),error=rep("",B),
    stringsAsFactors=FALSE)
  set.seed(seed+100L)
  rep_seeds <- sample.int(.Machine$integer.max,B)
  for(b in seq_len(B)) {
    set.seed(rep_seeds[b]); dat <- semi_generate(n,par,c0)
    attr(dat,"tau") <- par$tau
    out$censor_fraction[b] <- mean(dat$R==0L)
    fold_seed <- as.integer(rep_seeds[b] %% 1000000L) + 101L
    fit <- tryCatch(semi_one_step(dat,V=V,M_q=M_q,g_floor=g_floor,
      seed=fold_seed,piA=par$piA,mu_method=mu_method),
      error=function(e)e)
    if(inherits(fit,"error")) {
      out$error[b] <- conditionMessage(fit)
      if(verbose) cat("Replication",b,"failed:",out$error[b],"\n")
    } else {
      out$full[b] <- fit$full;out$se_full[b] <- fit$se_full
      out$ipcw[b] <- fit$ipcw;out$se_ipcw[b] <- fit$se_ipcw
      out$plugin[b] <- fit$plugin;out$min_G[b] <- fit$minG
      out$ipcw_se_error[b] <- if(is.null(fit$ipcw_se_error) ||
        is.na(fit$ipcw_se_error)) "" else fit$ipcw_se_error
    }
    if(oracle) {
      ofit <- tryCatch(semi_one_step(dat,V=V,M_q=M_q,g_floor=g_floor,
        seed=fold_seed,piA=par$piA,mu_method=mu_method,
        oracle=TRUE,par=par,c0=c0),error=function(e)e)
      if(!inherits(ofit,"error")) {
        out$oracle[b] <- ofit$full
        out$se_oracle[b] <- ofit$se_full
      }
    }
    if(verbose && (b%%10L==0L || b==B))
      cat("Completed",b,"/",B,"[",scenario,",n=",n,"]\n")
  }
  summary <- rbind(
    FULL_EIF=semi_summarize(out$full,out$se_full,truth),
    IPCW=semi_summarize(out$ipcw,out$se_ipcw,truth),
    PLUGIN=semi_summarize(out$plugin,NULL,truth))
  if(oracle) summary <- rbind(summary,
    ORACLE_HAZARDS=semi_summarize(out$oracle,out$se_oracle,truth))
  ok <- is.finite(out$ipcw) & is.finite(out$full)
  RE <- if(sum(ok)>1) var(out$ipcw[ok])/var(out$full[ok])
    else NA_real_
  list(summary=summary,replicates=out,truth=truth,
    truth_mcse=truth_mcse,c0=c0,true_min_G=true_min_G,scenario=scenario,
    n=n,B=B,V=V,M_q=M_q,seed=seed,kappa=RIDGE_KAPPA,
    relative_efficiency=RE,
    achieved_censoring=mean(out$censor_fraction,na.rm=TRUE),
    full_failures=sum(!is.finite(out$full)),
    ipcw_se_failures=sum(is.finite(out$ipcw)&!is.finite(out$se_ipcw)))
}

semi_export <- function(res,prefix) {
  s <- data.frame(method=rownames(res$summary),res$summary,
    row.names=NULL,scenario=res$scenario,n=res$n,
    achieved_censoring=res$achieved_censoring,
    true_min_G=res$true_min_G,
    truth_mcse=res$truth_mcse,
    relative_efficiency=res$relative_efficiency)
  write.csv(s,paste0(prefix,"_summary.csv"),row.names=FALSE)
  write.csv(res$replicates,paste0(prefix,"_replicates.csv"),row.names=FALSE)
  saveRDS(res,paste0(prefix,"_results.rds"))
}

semi_sanity_check <- function(seed=2026) {
  set.seed(seed)
  p1 <- c(1,1,5,5,3);s1 <- c(.5,.5,1,5,2)
  p0 <- c(2,1,5,5,4);s0 <- c(1,.5,2,5,1)
  stopifnot(max(abs(semi_h(p1,s1,p0,s0)+
                    semi_h(p0,s0,p1,s1)-1))<1e-12)
  par <- semi_parameters("moderate")
  c0 <- semi_calibrate_censoring(par,M=10000L)
  d <- semi_generate(200L,par,c0)
  attr(d,"tau") <- par$tau
  stopifnot(all(is.na(d$Yp[d$R==0])),all(is.na(d$Ys[d$R==0])))
  mods <- semi_fit_models(d)
  # Check unpenalized intercept score and penalized slope scores.
  for(a in 0:1) {
    da <- d[d$A==a,,drop=FALSE]
    ped <- semi_ped(da,"C")
    fc <- mods[[a+1]]$C
    Z <- cbind(1,ped$X,ped$N2)
    lam <- semi_predict_rate(fc,ped$X,ped$N2)
    score <- colSums(Z*as.numeric(ped$event-ped$exposure*lam)) -
      c(0,RIDGE_KAPPA*fc$bX,RIDGE_KAPPA*fc$bN2)
    if(!fc$fallback) stopifnot(max(abs(score))<1e-4)
  }
  proj <- semi_make_projection(d,mods,par$tau)
  stopifnot(is.finite(proj$psi),proj$psi>=0,proj$psi<=1)
  result <- semi_one_step(d,V=5,M_q=64,piA=par$piA,seed=seed+1)
  stopifnot(is.finite(result$full),is.finite(result$se_full),
            is.finite(result$ipcw),is.finite(result$se_ipcw))
  cat("Ridge semi-competing sanity checks passed.\n")
  invisible(TRUE)
}

# RUN ALL THREE SCENARIOS; SAVE CSV + RDS.
semi_sanity_check()
results <- list()
for(sc in RUN_SCENARIOS) {
  res <- run_semicomp_simulation(scenario=sc)
  results[[sc]] <- res
  cat("\nScenario:",sc,"; ridge kappa=",RIDGE_KAPPA,"\n")
  print(res$summary)
  cat("True minimum G over support:",res$true_min_G,"\n")
  cat("Relative efficiency IPCW/FULL_EIF:",
      res$relative_efficiency,"\n")
  cat("FULL_EIF failures:",res$full_failures,
      "; IPCW SE failures:",res$ipcw_se_failures,"\n")
  if(any(is.finite(res$replicates$se_full)))
    print(quantile(res$replicates$se_full,
      c(.5,.9,.95,.99,1),na.rm=TRUE))
  good_G <- res$replicates$min_G[is.finite(res$replicates$min_G)]
  cat("Min fitted G:",if(length(good_G)) min(good_G) else NA_real_,"\n")
  prefix <- paste0("semicomp_ridge_",sc,"_n",SIM_N,
                   "_B",SIM_B,"_seed",SIM_SEED)
  semi_export(res,prefix)
}

