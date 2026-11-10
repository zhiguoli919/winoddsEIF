# ========================================================================
# SCALAR WIN-PROBABILITY SIMULATIONS: POOLED CENSORING RIDGE MODEL
# Complete base-R implementation. No bootstrap.
#
# CHANGE relative to scalar_ridge_all_scenarios_archived.R:
#   * event hazards remain ridge-regularized and arm-specific;
#   * censoring hazard fitted by ONE ridge Poisson model with A + X,
#     where treatment and X slopes are shared across the pooled sample;
#   * IPCW score/Hessian correction uses the pooled 3-coefficient model.
# All event/censoring fits are cross-fitted on the training folds only.
#
# All three scenarios share kappa=2; NO hazard/weight truncation.
# DGP censoring has no treatment effect, but fitted model allows one.
#
# REVISED MAIN STRONGEST SCENARIO (Nov 2026):
#   betaX=1.75 remains unchanged; gammaX is reduced from 1.25 to 0.80;
#   target censoring fraction remains 0.40. This yields calibrated
#   c0 about 0.16585549 and true min G about 0.01645 on X in [-2,2],
#   t in [0,5]; the target win probability is unchanged (~0.543225981).
#   Previous gammaX=1.25 scenario was a near-positivity stress test.
#
# All three scenarios and all estimators remain available in this ONE
# submission-ready script. RUN_SCENARIOS defaults to 'strongest' ONLY
# so that prior moderate/stronger experiments are not rerun.
# To run all three: RUN_SCENARIOS <- c('moderate','stronger','strongest')
#
# Default n=300, B=1000. Use B=5000 for a more informative RE check.
# SEED_STREAM_B=5000 preserves the old replication seed stream for B<=5000,
# assuming unchanged RNG settings; strongest DGP itself has changed.
# ========================================================================
SIM_N           <- 800L   # total n, approximately 150 per arm
SIM_B           <- 5000L  # change to 5000L for full comparison
SIM_V           <- 5L
SIM_SEED        <- 2027L
SEED_STREAM_B   <- 5000L  # same master seed stream as old B=5000 run
SIM_ORACLE      <- TRUE
RIDGE_KAPPA     <- 2.0
RUN_SCENARIOS <- c("moderate", "stronger", "strongest")
RUN_MAIN        <- TRUE   # set FALSE to source functions without running

SCENARIOS <- list(
  moderate = list(tau=5,lambda0=.10,betaA=-.40,betaX=.75,
                  piA=.50,gammaX=.65,target=.30),
  stronger = list(tau=5,lambda0=.10,betaA=-.40,betaX=1.00,
                  piA=.50,gammaX=.80,target=.40),
  strongest = list(tau=5,lambda0=.10,betaA=-.40,betaX=1.75,
                   piA=.50,gammaX=.80,target=.40)  # REVISED strongest
)
MIN_RATE <- 1e-12   # numerical safeguard ONLY; not a substantive rate cap

rX <- function(n) qnorm(runif(n,pnorm(-2),pnorm(2)))
rateT <- function(a,x,p) p$lambda0*exp(p$betaA*a+p$betaX*x)
rateC <- function(x,c0,p) c0*exp(p$gammaX*x)

gl_nodes <- function(n=80L) {
  k <- seq_len(n-1L)
  J <- matrix(0,n,n)
  J[cbind(k,k+1L)] <- k/sqrt(4*k*k-1)
  J[cbind(k+1L,k)] <- k/sqrt(4*k*k-1)
  e <- eigen(J,symmetric=TRUE)
  list(x=qnorm(pnorm(-2)+(pnorm(2)-pnorm(-2))*(e$values+1)/2),
       w=e$vectors[1L,]^2)
}

truth_psi <- function(p,n_nodes=80L) {
  q <- gl_nodes(n_nodes)
  l1 <- rateT(1,q$x,p); l0 <- rateT(0,q$x,p)
  L <- outer(l1,l0,"+")
  h <- outer(rep(1,length(l1)),l0,"*")/L*(-expm1(-p$tau*L)) +
       .5*exp(-p$tau*L)
  sum(outer(q$w,q$w,"*")*h)
}

censor_fraction <- function(c0,p,q=gl_nodes()) {
  lc <- rateC(q$x,c0,p)
  f <- function(a) {
    lt <- rateT(a,q$x,p)
    sum(q$w*lc/(lt+lc)*(-expm1(-(lt+lc)*p$tau)))
  }
  (1-p$piA)*f(0)+p$piA*f(1)
}

calibrate_censor <- function(p) {
  f <- function(c0) censor_fraction(c0,p)-p$target
  hi <- .5
  while(f(hi)<0) hi <- 2*hi
  uniroot(f,c(0,hi),tol=1e-12)$root
}

# For gammaX >= 0 and X in [-2,2], t in [0,tau], the smallest
# TRUE censoring survival occurs at X=2, t=tau. This is a DGP
# diagnostic, NOT truncation of the fitted censoring distribution.
true_min_G <- function(p,c0) {
  xmax <- if(p$gammaX >= 0) 2 else -2
  exp(-p$tau*rateC(xmax,c0,p))
}

# Prespecified practical-positivity threshold for MAIN scalar settings.
MAIN_G_THRESHOLD <- 0.01

generate_data <- function(n,p,c0) {
  X <- rX(n); A <- rbinom(n,1,p$piA)
  T <- rexp(n,rate=rateT(A,X,p))
  C <- rexp(n,rate=rateC(X,c0,p))
  K <- pmin(T,p$tau)
  data.frame(X=X,A=A,s=pmin(K,C),R=as.integer(C>=K),
             Delta=as.integer(T<=C & T<=p$tau))
}

# The penalized counting-process log-likelihood (up to constants):
# ell(alpha,beta) = sum(event)*(alpha) + sum(event*X)*beta
#                   - sum(s*exp(alpha+beta*X)) - kappa*beta^2/2.
# For each beta, alpha is explicitly profiled out. The profile score
# decreases strictly for kappa>0, so a single uniroot() is sufficient.
# This avoids external packages and unstable two-dimensional fitting.
fit_rate <- function(d,event_name) {
  event <- d[[event_name]]
  x <- d$X
  E <- sum(event)
  if(nrow(d)<5L || !is.finite(E) || E<1L ||
     any(!is.finite(d$s)) || any(d$s<=0)) {
    return(list(alpha=log(max(E/sum(d$s),MIN_RATE)),
                beta=0,pooled=max(E/sum(d$s),MIN_RATE),
                fallback=TRUE))
  }
  kappa <- RIDGE_KAPPA
  stopifnot(is.finite(kappa),kappa>0)
  SX <- sum(event*x)

  exposure_weights <- function(beta) {
    z <- beta*x
    em <- max(z)
    w <- d$s * exp(z-em)
    list(w=w, maxz=em)
  }
  profile_score <- function(beta) {
    ew <- exposure_weights(beta)
    SX - E*sum(ew$w*x)/sum(ew$w) - kappa*beta
  }
  bound <- max(8,4*E/kappa+1)
  beta <- uniroot(profile_score,c(-bound,bound),tol=1e-10)$root
  ew <- exposure_weights(beta)
  alpha <- log(E)-log(sum(ew$w))-ew$maxz
  if(!is.finite(alpha) || !is.finite(beta))
    stop("Non-finite penalized hazard coefficient")
  list(alpha=alpha,beta=beta,pooled=NA_real_,fallback=FALSE)
}

# Unlike the previous scripts, we do not cap the fitted hazard at 0.5,
# 1, 4, etc. The fixed ridge penalty stabilizes the coefficient estimate.
predict_rate <- function(f,x,clipped=TRUE) {
  if(f$fallback) out <- rep(f$pooled,length(x)) else {
    eta <- f$alpha + f$beta*as.numeric(x)
    if(any(!is.finite(eta)) || any(eta>50))
      stop("Extreme fitted log-hazard: inspect replication")
    out <- exp(eta)
  }
  if(any(!is.finite(out)) || any(out<=0))
    stop("Invalid penalized hazard prediction")
  pmax(out,MIN_RATE)
}

# The POOLED censoring hazard has log rate alpha + beta_A*A + beta_X*X.
# Penalize BOTH slopes but not the intercept, with the same kappa=2.
# Profile out alpha and maximize the strictly concave two-slope profile
# log-likelihood by Newton's method with backtracking. No packages used.
fit_censor_pooled <- function(train) {
  ev <- as.numeric(1L-train$R)
  E <- sum(ev)
  S <- sum(train$s)
  if(nrow(train)<5L || !is.finite(S) || S<=0 ||
     any(!is.finite(train$s)) || any(train$s<=0) || E<1) {
    raw <- if(is.finite(S) && S>0) max(E/S,MIN_RATE) else MIN_RATE
    return(list(alpha=log(raw), beta=c(A=0,X=0),
                pooled=raw, fallback=TRUE))
  }
  Z <- cbind(A=as.numeric(train$A), X=as.numeric(train$X))
  SX <- as.vector(crossprod(Z,ev))
  kappa <- RIDGE_KAPPA
  stopifnot(is.finite(kappa),kappa>0)

  evaluate <- function(beta) {
    eta <- as.vector(Z %*% beta)
    m <- max(eta)
    w <- train$s*exp(eta-m)
    total <- sum(w)
    if(!is.finite(total) || total<=0)
      stop("Invalid exposure sum in pooled censoring model")
    prob <- w/total
    zbar <- as.vector(crossprod(Z,prob))
    centered <- sweep(Z,2L,zbar,"-")
    covariance <- crossprod(centered,centered*as.vector(prob))
    score <- SX-E*zbar-kappa*beta
    information <- E*covariance+diag(kappa,2L)
    obj <- sum(SX*beta)-E*(log(total)+m)-kappa*sum(beta^2)/2
    alpha <- log(E)-log(total)-m
    list(obj=obj, score=score, information=information, alpha=alpha)
  }

  beta <- c(0,0)
  for(iter in seq_len(100L)) {
    z <- evaluate(beta)
    if(max(abs(z$score)) < 1e-7) break
    step <- as.vector(solve(z$information,z$score))
    if(any(!is.finite(step))) stop("Invalid pooled censoring Newton step")
    scale <- 1
    accepted <- FALSE
    for(j in seq_len(40L)) {
      proposal <- beta+scale*step
      new_z <- evaluate(proposal)
      # Permit a tiny objective rounding error close to the optimum.
      sufficient_increase <- new_z$obj >=
        z$obj + 1e-4*scale*sum(z$score*step)
      rounding_safe <- (max(abs(new_z$score)) < max(abs(z$score)) &&
                        new_z$obj >= z$obj-1e-12*(1+abs(z$obj)))
      if(is.finite(new_z$obj) &&
         (sufficient_increase || rounding_safe)) {
        beta <- proposal
        accepted <- TRUE
        break
      }
      scale <- scale/2
    }
    if(!accepted) stop("Pooled censoring line search failed")
  }
  z <- evaluate(beta)
  if(max(abs(z$score)) > 1e-6 ||
     any(!is.finite(c(z$alpha,beta))))
    stop("Pooled censoring model did not converge")
  names(beta) <- c("A","X")
  list(alpha=z$alpha,beta=beta,pooled=NA_real_,fallback=FALSE)
}

predict_censor_pooled <- function(f,a,x) {
  nn <- max(length(a),length(x))
  if(f$fallback) return(rep(f$pooled,nn))
  eta <- f$alpha + f$beta["A"]*a + f$beta["X"]*x
  if(any(!is.finite(eta)) || any(eta>50))
    stop("Extreme pooled censoring log-hazard: inspect replication")
  out <- exp(eta)
  if(any(!is.finite(out)) || any(out<=0))
    stop("Invalid pooled censoring hazard")
  pmax(out,MIN_RATE)
}

fit_models <- function(train,p,c0,oracle=FALSE) {
  if(oracle) return(list(
    T=function(a,x) rateT(a,x,p),
    C=function(a,x) rateC(x,c0,p),
    cfit=NULL,fallbacks=0L))
  tf <- lapply(0:1,function(a)
    fit_rate(train[train$A==a,,drop=FALSE],"Delta"))
  cf <- fit_censor_pooled(train)
  list(T=function(a,x) predict_rate(tf[[a+1L]],x),
       C=function(a,x) predict_censor_pooled(cf,a,x),
       cfit=cf,
       fallbacks=sum(vapply(tf,function(f) f$fallback,logical(1)))+
                 as.integer(cf$fallback))
}

# Normalized IPCW law of revealed outcomes for one treatment arm.
outcome_law <- function(d,a,fit,tau) {
  z <- d[d$A==a & d$R==1L,,drop=FALSE]
  if(nrow(z)<2L) stop("Too few complete outcomes")
  lw <- fit$C(a,z$X)*z$s
  w <- exp(lw-max(lw)); w <- w/sum(w)
  list(y=z$s,w=w,is_tau=(z$s==tau),
       ess=1/sum(w^2))
}

# The only outcome atom is at tau, so comparisons use exact equality.
h_projection <- function(a,y,F) {
  if(a==1L) return(sum(F$w*(as.numeric(y>F$y)+.5*(y==F$y))))
  sum(F$w*(as.numeric(y<F$y)+.5*(y==F$y)))
}

psi_plugin <- function(F1,F0)
  sum(F1$w*vapply(F1$y,function(y) h_projection(1L,y,F0),numeric(1)))

# E[h_a(min(T,tau)) | T>=u, X=x,A=a] under fitted exponential event model.
conditional_h <- function(a,u,lt,F,tau) {
  y <- F$y; s <- exp(-lt*pmax(y-u,0))
  if(a==1L) {
    val <- s
    val[F$is_tau] <- .5*exp(-lt*(tau-u))
  } else {
    val <- ifelse(y>u,-expm1(-lt*(y-u)),0)
    val[F$is_tau] <- 1-.5*exp(-lt*(tau-u))
  }
  sum(F$w*val)
}

# Exact integral int_0^s (q(u)/G(u))lambda_C du.
compensator <- function(a,s,lt,lc,F,tau,mu) {
  if(s<=0) return(0)
  y <- F$y; t <- pmin(y,s)
  J <- lc/(lc+lt)*exp(lc*t-lt*(y-t))*(-expm1(-(lc+lt)*t))
  if(a==1L) {
    I <- J
    ii <- which(s>y)
    if(length(ii)) I[ii] <- I[ii]+exp(lc*s)-exp(lc*y[ii])
    I[F$is_tau] <- .5*J[F$is_tau]
  } else {
    I <- expm1(lc*t)-J
    I[F$is_tau] <- expm1(lc*s)-.5*J[F$is_tau]
  }
  sum(F$w*I)-mu*expm1(lc*s)
}

# First-order direct SE of the normalized IPCW U-statistic.
# For the POOLED censoring regression the nuisance score/Hessian are
# three-dimensional (intercept, A, X). Because each fold's censoring fit
# is trained on ALL non-validation subjects, its score correction is
# distributed across ALL those subjects, not separately by arm.
#
# Weight normalization and variation of BOTH empirical outcome laws enter
# through the projected kernel h_i-psi. Fixed kappa produces O(1/n) bias;
# this remains a first-order asymptotic variance approximation and is
# evaluated by empirical coverage in the simulations.
ipcw_direct_se <- function(d,fold,cfits,weights,h1,h0,psi,ix1,ix0) {
  n <- nrow(d); V <- max(fold)
  IF <- numeric(n); h <- numeric(n)
  h[ix1] <- h1; h[ix0] <- h0
  ii <- c(ix1,ix0)
  IF[ii] <- n*weights[ii]*(h[ii]-psi)
  nuisance <- numeric(n)
  for(v in seq_len(V)) {
    fit <- cfits[[v]]
    if(is.null(fit) || fit$fallback)
      stop("Pooled censor model fallback: no valid IPCW SE")
    tr <- which(fold!=v)
    va <- which(fold==v & d$R==1L)
    ntr <- length(tr)

    D <- numeric(3L)
    if(length(va)) {
      zva <- cbind(1,d$A[va],d$X[va])
      lam_va <- predict_censor_pooled(fit,d$A[va],d$X[va])
      D <- as.vector(crossprod(
        zva,weights[va]*(h[va]-psi)*d$s[va]*lam_va))
    }

    ztr <- cbind(1,d$A[tr],d$X[tr])
    lam_tr <- predict_censor_pooled(fit,d$A[tr],d$X[tr])
    exposure <- lam_tr*d$s[tr]
    H <- (crossprod(ztr,ztr*as.vector(exposure)) +
          diag(c(0,RIDGE_KAPPA,RIDGE_KAPPA),3L))/ntr
    if(any(!is.finite(H)) || rcond(H)<1e-10)
      stop("Pooled censoring information matrix singular")
    score <- ztr*as.vector((1L-d$R[tr])-exposure)
    nuisance[tr] <- nuisance[tr]+(n/ntr)*
      as.vector(score %*% solve(H,D))
  }
  IF <- IF+nuisance
  if(any(!is.finite(IF))) stop("Non-finite IPCW IF")
  list(se=sqrt(mean((IF-mean(IF))^2)/n))
}

estimate_once <- function(d,p,c0,fold,oracle=FALSE) {
  n <- nrow(d); V <- max(fold); tau <- p$tau
  phi <- numeric(n); phi0 <- numeric(n)
  plug <- numeric(V); ww <- numeric(n)
  G <- numeric(n); ess1 <- ess0 <- numeric(V)
  cfits <- vector("list",V); nfallback <- 0L
  for(v in seq_len(V)) {
    train <- d[fold!=v,,drop=FALSE]
    model <- fit_models(train,p,c0,oracle)
    cfits[[v]] <- model$cfit
    nfallback <- nfallback+model$fallbacks
    F1 <- outcome_law(train,1L,model,tau)
    F0 <- outcome_law(train,0L,model,tau)
    ess1[v] <- F1$ess; ess0[v] <- F0$ess
    plug[v] <- psi_plugin(F1,F0)
    for(i in which(fold==v)) {
      a <- d$A[i]; x <- d$X[i]; s <- d$s[i]
      lt <- model$T(a,x); lc <- model$C(a,x)
      mu1 <- conditional_h(1L,0,model$T(1L,x),F0,tau)
      mu0 <- conditional_h(0L,0,model$T(0L,x),F1,tau)
      mu <- if(a==1L) mu1 else mu0
      F <- if(a==1L) F0 else F1
      pa <- if(a==1L) p$piA else 1-p$piA
      invG <- exp(lc*s)
      if(!is.finite(invG)) stop("Inverse censoring overflow")
      ipart <- if(d$R[i]==1L)
        invG*(h_projection(a,s,F)-mu) else 0
      qjump <- if(d$R[i]==0L)
        invG*(conditional_h(a,s,lt,F,tau)-mu) else 0
      aug <- qjump-compensator(a,s,lt,lc,F,tau,mu)
      base <- mu1+mu0-2*plug[v]
      phi[i] <- base+(ipart+aug)/pa
      phi0[i] <- base+ipart/pa
      G[i] <- 1/invG
      if(d$R[i]==1L) ww[i] <- invG/pa
    }
  }
  plug_avg <- sum(tabulate(fold,nbins=V)*plug)/n
  out <- list(
    full=plug_avg+mean(phi),
    se_full=sqrt(mean((phi-mean(phi))^2)/n),
    q0=plug_avg+mean(phi0),
    se_q0=sqrt(mean((phi0-mean(phi0))^2)/n),
    plugin=plug_avg,
    min_G=min(G),
    min_ess_treated=min(ess1),
    min_ess_control=min(ess0),
    fallback=nfallback)
  if(oracle) return(out)

  ix1 <- which(d$A==1L & d$R==1L)
  ix0 <- which(d$A==0L & d$R==1L)
  if(!length(ix1) || !length(ix0)) stop("No complete outcomes in one arm")
  ww[ix1] <- ww[ix1]/sum(ww[ix1])
  ww[ix0] <- ww[ix0]/sum(ww[ix0])
  F0 <- list(y=d$s[ix0],w=ww[ix0])
  F1 <- list(y=d$s[ix1],w=ww[ix1])
  h1 <- vapply(d$s[ix1],function(y) h_projection(1L,y,F0),numeric(1))
  h0 <- vapply(d$s[ix0],function(y) h_projection(0L,y,F1),numeric(1))
  out$ipcw <- sum(ww[ix1]*h1)
  out$se_ipcw <- NA_real_
  out$ipcw_error <- ""
  ipse <- tryCatch(ipcw_direct_se(d,fold,cfits,ww,h1,h0,
                                  out$ipcw,ix1,ix0),error=function(e)e)
  if(inherits(ipse,"error")) out$ipcw_error <- conditionMessage(ipse)
  else {
    out$se_ipcw <- ipse$se
  }
  out
}

summarize_one <- function(x,se,truth) {
  good <- is.finite(x)
  z <- x[good]; N <- length(z)
  seok <- good & is.finite(se) & se>0
  s <- se[seok]
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

run_simulation <- function(scenario,n=SIM_N,B=SIM_B,V=SIM_V,
                           seed=SIM_SEED,oracle=SIM_ORACLE,progress=50L) {
  p <- SCENARIOS[[scenario]]
  if(is.null(p)) stop("Unknown scenario")
  if(n<40 || V<2 || V>n/4 || B<2) stop("Invalid n, B, or V")
  truth <- truth_psi(p); c0 <- calibrate_censor(p)
  g0 <- true_min_G(p,c0)
  if(g0 < MAIN_G_THRESHOLD)
    stop(sprintf("True min G=%.8g violates specified main-scenario threshold %.3f",
                 g0, MAIN_G_THRESHOLD))
  cat(sprintf(paste0("%s: n=%d, B=%d; psi=%.9f; c0=%.9f; ",
                     "target censor=%.2f; min true G=%.6f\n"),
              scenario,n,B,truth,c0,p$target,g0))
  r <- data.frame(replication=seq_len(B),full=NA_real_,se_full=NA_real_,
      q0=NA_real_,se_q0=NA_real_,ipcw=NA_real_,se_ipcw=NA_real_,
      plugin=NA_real_,oracle=NA_real_,se_oracle=NA_real_,
      min_G=NA_real_,min_ess_treated=NA_real_,min_ess_control=NA_real_,
      fallback=NA_integer_,censor_fraction=NA_real_,
      ipcw_se_error=rep("",B),error=rep("",B),stringsAsFactors=FALSE)
  set.seed(seed+100L)
  rep_seeds <- sample.int(.Machine$integer.max,max(B,SEED_STREAM_B))
  rep_seeds <- rep_seeds[seq_len(B)]
  for(b in seq_len(B)) {
    set.seed(rep_seeds[b])
    d <- generate_data(n,p,c0)
    r$censor_fraction[b] <- mean(d$R==0L)
    set.seed(as.integer(rep_seeds[b] %% 1000000L)+101L)
    fold <- sample(rep(seq_len(V),length.out=n))
    z <- tryCatch(estimate_once(d,p,c0,fold),error=function(e)e)
    if(inherits(z,"error")) {
      r$error[b] <- conditionMessage(z)
    } else {
      r$full[b] <- z$full; r$se_full[b] <- z$se_full
      r$q0[b] <- z$q0; r$se_q0[b] <- z$se_q0
      r$ipcw[b] <- z$ipcw; r$se_ipcw[b] <- z$se_ipcw
      r$plugin[b] <- z$plugin; r$min_G[b] <- z$min_G
      r$min_ess_treated[b] <- z$min_ess_treated
      r$min_ess_control[b] <- z$min_ess_control
      r$fallback[b] <- z$fallback
      r$ipcw_se_error[b] <- z$ipcw_error
    }
    if(oracle) {
      o <- tryCatch(estimate_once(d,p,c0,fold,oracle=TRUE),
                    error=function(e)e)
      if(!inherits(o,"error")) {
        r$oracle[b] <- o$full; r$se_oracle[b] <- o$se_full
      }
    }
    if(progress>0 && (b%%progress==0 || b==B))
      cat(sprintf("Completed %d/%d [%s,n=%d]\n",b,B,scenario,n))
  }
  summary <- rbind(
    FULL_EIF=summarize_one(r$full,r$se_full,truth),
    Q0_DIAGNOSTIC=summarize_one(r$q0,r$se_q0,truth),
    IPCW=summarize_one(r$ipcw,r$se_ipcw,truth),
    PLUGIN=summarize_one(r$plugin,rep(NA_real_,B),truth))
  if(oracle) summary <- rbind(summary,
    ORACLE_HAZARDS=summarize_one(r$oracle,r$se_oracle,truth))
  good <- is.finite(r$full) & is.finite(r$ipcw)
  RE <- if(sum(good)>1L) var(r$ipcw[good])/var(r$full[good])
        else NA_real_
  list(summary=summary,replicates=r,truth=truth,c0=c0,
       true_min_G=g0,scenario_parameters=p,scenario=scenario,
       n=n,B=B,V=V,seed=seed,seed_stream_B=SEED_STREAM_B,
       ridge_kappa=RIDGE_KAPPA,censor_model="pooled_A_plus_X",
       relative_efficiency=RE,
       full_failures=sum(!is.finite(r$full)),
       ipcw_se_failures=sum(is.finite(r$ipcw)&!is.finite(r$se_ipcw)))
}

# Numerical check of the exact censoring compensator integral.
sanity_check <- function() {
  # Verify all three DGPs satisfy the main-scenario minimum TRUE G.
  for(nm in c("moderate","stronger","strongest")) {
    pp <- SCENARIOS[[nm]]
    cc <- calibrate_censor(pp)
    stopifnot(true_min_G(pp,cc) >= MAIN_G_THRESHOLD,
              abs(censor_fraction(cc,pp)-pp$target) < 1e-8)
  }
  # These reference numbers pertain to the REVISED strongest DGP.
  p_new <- SCENARIOS$strongest
  c_new <- calibrate_censor(p_new)
  stopifnot(abs(c_new-0.16585549) < 2e-6,
            abs(true_min_G(p_new,c_new)-0.016450) < 1e-4,
            abs(truth_psi(p_new)-0.543225981) < 1e-6)
  p <- SCENARIOS$moderate
  stopifnot(abs(truth_psi(p,80L)-truth_psi(p,100L))<1e-8)
  c0 <- calibrate_censor(p)
  stopifnot(abs(censor_fraction(c0,p)-p$target)<1e-8)
  F <- list(y=c(.6,1.7,p$tau),w=c(.2,.3,.5),
            is_tau=c(FALSE,FALSE,TRUE))
  for(a in 0:1) {
    lt <- .12; lc <- .15; s <- 2.1
    mu <- conditional_h(a,0,lt,F,p$tau)
    approx <- integrate(function(u) vapply(u,function(v)
      (conditional_h(a,v,lt,F,p$tau)-mu)*lc*exp(lc*v),numeric(1)),
      0,s,rel.tol=1e-10)$value
    exact <- compensator(a,s,lt,lc,F,p$tau,mu)
    stopifnot(abs(approx-exact)<1e-7)
  }
  # Penalized score equations for the arm-specific event rate.
  toy <- data.frame(X=c(-1.4,-.9,-.4,0,.4,.7,1.1,1.5),
                    s=c(1,2,1.5,3,1,.7,2,1.3),
                    ev=c(0,1,0,0,1,0,1,0),
                    A=c(0,0,1,0,1,1,0,1),
                    R=c(1,0,1,1,0,1,0,1))
  f <- fit_rate(toy,"ev")
  stopifnot(!f$fallback)
  Z <- cbind(1,toy$X)
  rate <- predict_rate(f,toy$X)
  score <- colSums(Z*as.numeric(toy$ev-rate*toy$s)) -
           c(0,RIDGE_KAPPA*f$beta)
  stopifnot(max(abs(score))<1e-7)

  # Test the pooled ridge censoring score AND the 3x3 information matrix.
  cf <- fit_censor_pooled(toy)
  stopifnot(!cf$fallback)
  Zc <- cbind(1,toy$A,toy$X)
  lc <- predict_censor_pooled(cf,toy$A,toy$X)
  beta_c <- c(cf$alpha,cf$beta)
  score_c <- as.vector(crossprod(Zc,(1-toy$R)-lc*toy$s)) -
    c(0,RIDGE_KAPPA*cf$beta)
  stopifnot(max(abs(score_c))<1e-6,
            all(is.finite(beta_c)))
  # Crossfit model: event remains arm-specific, censor is pooled.
  model <- fit_models(toy,p,c0,oracle=FALSE)
  stopifnot(length(model$cfit$beta)==2L,
            identical(names(model$cfit$beta),c("A","X")))
  cat("Sanity checks passed (including ridge score equations).\n")
}

# ========================================================================
# RECONSTRUCT ANY REPLICATION with identical previous B=5000 seed stream.
# For example (after setting RUN_MAIN=FALSE above and sourcing this file):
#   z <- reproduce_replication("stronger",1798L,n=300L)
#   z$estimate[c("full","se_full","ipcw","min_G")]
# ========================================================================
reproduce_replication <- function(scenario="stronger",b=1798L,
                                  n=SIM_N,seed=SIM_SEED,V=SIM_V) {
  p <- SCENARIOS[[scenario]]
  if(is.null(p)) stop("Invalid scenario")
  stopifnot(b>=1L,b==as.integer(b),n>=40,V>=2L)
  c0 <- calibrate_censor(p)
  set.seed(seed+100L)
  seeds <- sample.int(.Machine$integer.max,max(b,SEED_STREAM_B))
  rep_seed <- seeds[b]
  set.seed(rep_seed)
  d <- generate_data(n,p,c0)
  set.seed(as.integer(rep_seed %% 1000000L)+101L)
  fold <- sample(rep(seq_len(V),length.out=n))
  list(estimate=estimate_once(d,p,c0,fold),
       data=d,fold=fold,scenario=scenario,replication=b,seed=seed)
}

# ========================================================================
# MAIN RUN: all 3 scenarios AVAILABLE, strongest ONLY by default.
# For submission, set RUN_SCENARIOS to all three when desired.
# For n=600 change SIM_N; for B=5000 change SIM_B.
# If running the script in an interactive session, results is populated.
# ========================================================================
sanity_check()
if(RUN_MAIN) {
  results <- list()
  for(sc in RUN_SCENARIOS) {
    res <- run_simulation(scenario=sc)
    results[[sc]] <- res
    cat("\nScenario: ",sc,"; ridge kappa=",RIDGE_KAPPA,
        "; pooled censoring A+X\n",sep="")
    print(res$summary)
    cat("Relative efficiency IPCW/FULL_EIF:",res$relative_efficiency,"\n")
    cat("FULL_EIF failures:",res$full_failures,
        "; IPCW SE failures:",res$ipcw_se_failures,"\n")
    print(quantile(res$replicates$se_full,
                   c(.5,.90,.95,.99,1),na.rm=TRUE))
    cat("Min TRUE G on covariate/time support:",res$true_min_G,"\n")
    cat("Min fitted G:",min(res$replicates$min_G,na.rm=TRUE),"\n")
    cat("Hazard fallbacks:",sum(res$replicates$fallback,na.rm=TRUE),"\n")
    prefix <- paste0("scalar_pooled_ridge_revised_",sc,"_n",SIM_N,"_B",SIM_B,
                     "_seed",SIM_SEED,"_kappa",
                     gsub("\\.","p",as.character(RIDGE_KAPPA)))
    write.csv(res$summary,paste0(prefix,"_summary.csv"))
    write.csv(res$replicates,paste0(prefix,"_replicates.csv"),row.names=FALSE)
    saveRDS(res,paste0(prefix,"_results.rds"))
  }
  combined <- do.call(rbind,lapply(names(results),function(sc) {
    tab <- results[[sc]]$summary
    data.frame(scenario=sc,method=rownames(tab),tab,
               row.names=NULL,check.names=FALSE)
  }))
  write.csv(combined,
      paste0("scalar_pooled_ridge_revised_n",SIM_N,"_B",SIM_B,
             "_seed",SIM_SEED,"_combined_summary.csv"),row.names=FALSE)
}
