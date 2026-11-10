# Scalar estimator machinery adapted from the supplied simulation script.
# These helpers are internal. Fixed ridge penalty matches the study code.

.wo_scalar_fit_rate <- function(d,event_name) {
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
.wo_scalar_predict_rate <- function(f,x,clipped=TRUE) {
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
.wo_scalar_fit_censor_pooled <- function(train) {
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

.wo_scalar_predict_censor_pooled <- function(f,a,x) {
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

.wo_scalar_fit_models <- function(train) {
  tf <- lapply(0:1,function(a)
    .wo_scalar_fit_rate(train[train$A==a,,drop=FALSE],"Delta"))
  cf <- .wo_scalar_fit_censor_pooled(train)
  list(T=function(a,x) .wo_scalar_predict_rate(tf[[a+1L]],x),
       C=function(a,x) .wo_scalar_predict_censor_pooled(cf,a,x),
       cfit=cf,
       fallbacks=sum(vapply(tf,function(f) f$fallback,logical(1)))+
                 as.integer(cf$fallback))
}

# Normalized IPCW law of revealed outcomes for one treatment arm.
.wo_scalar_outcome_law <- function(d,a,fit,tau) {
  z <- d[d$A==a & d$R==1L,,drop=FALSE]
  if(nrow(z)<2L) stop("Too few complete outcomes")
  lw <- fit$C(a,z$X)*z$s
  w <- exp(lw-max(lw)); w <- w/sum(w)
  list(y=z$s,w=w,is_tau=(z$s==tau),
       ess=1/sum(w^2))
}

# The only outcome atom is at tau, so comparisons use exact equality.
.wo_scalar_h_projection <- function(a,y,F) {
  if(a==1L) return(sum(F$w*(as.numeric(y>F$y)+.5*(y==F$y))))
  sum(F$w*(as.numeric(y<F$y)+.5*(y==F$y)))
}

.wo_scalar_psi_plugin <- function(F1,F0)
  sum(F1$w*vapply(F1$y,function(y) .wo_scalar_h_projection(1L,y,F0),numeric(1)))

# E[h_a(min(T,tau)) | T>=u, X=x,A=a] under fitted exponential event model.
.wo_scalar_conditional_h <- function(a,u,lt,F,tau) {
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
.wo_scalar_compensator <- function(a,s,lt,lc,F,tau,mu) {
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
.wo_scalar_ipcw_direct_se <- function(d,fold,cfits,weights,h1,h0,psi,ix1,ix0) {
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
      lam_va <- .wo_scalar_predict_censor_pooled(fit,d$A[va],d$X[va])
      D <- as.vector(crossprod(
        zva,weights[va]*(h[va]-psi)*d$s[va]*lam_va))
    }

    ztr <- cbind(1,d$A[tr],d$X[tr])
    lam_tr <- .wo_scalar_predict_censor_pooled(fit,d$A[tr],d$X[tr])
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

.wo_scalar_estimate_once <- function(d,p,fold) {
  n <- nrow(d); V <- max(fold); tau <- p$tau
  phi <- numeric(n)
  plug <- numeric(V); ww <- numeric(n)
  G <- numeric(n); ess1 <- ess0 <- numeric(V)
  cfits <- vector("list",V); nfallback <- 0L
  for(v in seq_len(V)) {
    train <- d[fold!=v,,drop=FALSE]
    model <- .wo_scalar_fit_models(train)
    cfits[[v]] <- model$cfit
    nfallback <- nfallback+model$fallbacks
    F1 <- .wo_scalar_outcome_law(train,1L,model,tau)
    F0 <- .wo_scalar_outcome_law(train,0L,model,tau)
    ess1[v] <- F1$ess; ess0[v] <- F0$ess
    plug[v] <- .wo_scalar_psi_plugin(F1,F0)
    for(i in which(fold==v)) {
      a <- d$A[i]; x <- d$X[i]; s <- d$s[i]
      lt <- model$T(a,x); lc <- model$C(a,x)
      mu1 <- .wo_scalar_conditional_h(1L,0,model$T(1L,x),F0,tau)
      mu0 <- .wo_scalar_conditional_h(0L,0,model$T(0L,x),F1,tau)
      mu <- if(a==1L) mu1 else mu0
      F <- if(a==1L) F0 else F1
      pa <- if(a==1L) p$piA else 1-p$piA
      invG <- exp(lc*s)
      if(!is.finite(invG)) stop("Inverse censoring overflow")
      ipart <- if(d$R[i]==1L)
        invG*(.wo_scalar_h_projection(a,s,F)-mu) else 0
      qjump <- if(d$R[i]==0L)
        invG*(.wo_scalar_conditional_h(a,s,lt,F,tau)-mu) else 0
      aug <- qjump-.wo_scalar_compensator(a,s,lt,lc,F,tau,mu)
      base <- mu1+mu0-2*plug[v]
      phi[i] <- base+(ipart+aug)/pa
      G[i] <- 1/invG
      if(d$R[i]==1L) ww[i] <- invG/pa
    }
  }
  plug_avg <- sum(tabulate(fold,nbins=V)*plug)/n
  out <- list(
    full=plug_avg+mean(phi),
    se_full=sqrt(mean((phi-mean(phi))^2)/n),
    plugin=plug_avg,
    min_G=min(G),
    min_ess_treated=min(ess1),
    min_ess_control=min(ess0),
    fallback=nfallback)
  ix1 <- which(d$A==1L & d$R==1L)
  ix0 <- which(d$A==0L & d$R==1L)
  if(!length(ix1) || !length(ix0)) stop("No complete outcomes in one arm")
  ww[ix1] <- ww[ix1]/sum(ww[ix1])
  ww[ix0] <- ww[ix0]/sum(ww[ix0])
  F0 <- list(y=d$s[ix0],w=ww[ix0])
  F1 <- list(y=d$s[ix1],w=ww[ix1])
  h1 <- vapply(d$s[ix1],function(y) .wo_scalar_h_projection(1L,y,F0),numeric(1))
  h0 <- vapply(d$s[ix0],function(y) .wo_scalar_h_projection(0L,y,F1),numeric(1))
  out$ipcw <- sum(ww[ix1]*h1)
  out$se_ipcw <- NA_real_
  out$ipcw_error <- ""
  ipse <- tryCatch(.wo_scalar_ipcw_direct_se(d,fold,cfits,ww,h1,h0,
                                  out$ipcw,ix1,ix0),error=function(e)e)
  if(inherits(ipse,"error")) out$ipcw_error <- conditionMessage(ipse)
  else {
    out$se_ipcw <- ipse$se
  }
  out
}

