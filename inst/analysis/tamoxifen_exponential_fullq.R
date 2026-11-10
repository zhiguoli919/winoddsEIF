# Tamoxifen trial: exponential working-hazard analysis (sensitivity)
# -------------------------------------------------------------------
# Treatment A=1: tamoxifen plus radiotherapy; A=0: tamoxifen alone. Primary endpoint: death; secondary endpoint: first local, axillary, or
# distant relapse strictly before death; same-day ties are death-first.
## Nuisance models: exponential working hazards, with history change after relapse for death and censoring. Arm-specific event models, pooled
# censoring model; ridge penalty kappa=2, five-fold cross-fitting, full q.
# Covariates: age10, log_size, hormone-receptor status.
# One-step SE is influence-function based; IPCW SE uses a first-order censoring-score correction.
# This is a sensitivity analysis, not the manuscript's main spline model.
# -------------------------------------------------------------------

RUN_HORIZONS <- c(4, 5)
N_FOLDS <- 5L
FOLD_SEED <- 2027L
RIDGE_KAPPA <- 2.0
M_Q <- 256L   # 16 x 16 conditional quantile grid
G_NUMERIC_FLOOR <- 1e-8  # numerical safety check, NOT weight truncation
OUTPUT_PREFIX <- "tamoxifen_fullq"

BASE_NAMES <- c("age10", "log_size", "hr")

# ------------------------ Prepare and verify patient data --------------------
prepare_tamoxifen <- function() {
  if (!requireNamespace("collett", quietly=TRUE))
    stop("Please install the CRAN package 'collett' first: install.packages('collett')")
  raw <- collett::tamoxifen
  required <- c("id","treat","age","size","hr","hist","hb",
                "lsurv","ls","asurv","as","dsurv","ds","tsurv","ts")
  if (!all(required %in% names(raw)))
    stop("Unexpected tamoxifen column names. Check collett version.")
  if (anyNA(raw[,required])) stop("Missing values in original dataset")
  if (anyDuplicated(raw$id)) stop("Duplicate patient IDs")
  if (any(!raw$treat %in% 0:1)) stop("Unexpected treatment coding")
  if (any(!raw$hr %in% 0:1)) stop("Unexpected hormone-receptor coding")
  if (any(raw$size<0) || any(raw$tsurv<=0)) stop("Invalid follow-up/size")

  relapse_times <- as.matrix(raw[,c("lsurv","asurv","dsurv")])
  relapse_events <- as.matrix(raw[,c("ls","as","ds")])
  if (any(relapse_times<=0) || any(!relapse_events %in% 0:1))
    stop("Invalid relapse follow-up times or indicators")

  # First observed relapse of any type, or first loss of observation of any
  # relapse type, whichever comes first. Thus we DO NOT search only among
  # event times (which might select an event after an earlier censoring time).
  first_time <- apply(relapse_times,1L,min)
  event_at_first <- rowSums(
    sweep(relapse_times,1L,first_time,`==`) & (relapse_events==1L)
  ) > 0L
  any_recorded <- rowSums(relapse_events==1L)>0L
  ambiguous <- sum(any_recorded & !event_at_first)
  if (ambiguous>0L)
    stop("Some relapses occur after another relapse type's follow-up ends: ",
         ambiguous, ". Composite outcome requires further handling.")

  td <- as.numeric(raw$tsurv)
  dd <- as.integer(raw$ts)
  same_day <- event_at_first & first_time==td & dd==1L
  # Treatment code in CRAN documentation:
  # treat=0 tamoxifen + radiotherapy; treat=1 tamoxifen alone.
  # A=1 is the combination treatment.
  d <- data.frame(id=raw$id, A=as.integer(raw$treat==0L),
                  td=td,dd=dd,tr=pmin(first_time,td),
                  dr=as.integer(event_at_first & first_time<td),
                  age=as.numeric(raw$age),size=as.numeric(raw$size),
                  hr=as.numeric(raw$hr),hist=as.integer(raw$hist),
                  hb=as.numeric(raw$hb))
  d <- d[order(d$id),,drop=FALSE]
  rownames(d) <- NULL

  if (any(d$dr==1L & d$tr>=d$td))
    stop("A relapse was counted at/after death or end of observation")
  if (any(d$dr==0L & d$tr<d$td))
    stop("A relapse process was censored earlier than death follow-up")
  if (any(d$tr<=0 | d$td<=0)) stop("Zero or negative observed time")
  if (any(!dd %in% 0:1)) stop("Unexpected death status")

  # Predeclared transformations; no data-derived score or imputation is
  # learned using the held-out fold. These are fixed mappings.
  d$age10 <- (d$age-60)/10
  d$log_size <- log1p(d$size)
  if (any(!is.finite(as.matrix(d[,BASE_NAMES,drop=FALSE]))))
    stop("Non-finite model predictor")

  cat("Patients:",nrow(d),"; treatment arms:",
      paste(names(table(d$A)),table(d$A),collapse="; "),"\n")
  cat("Same-day relapse/death handled as death-first:",sum(same_day),"\n")
  cat("Any first relapse strictly before death:",sum(d$dr),"\n")
  cat("Relapse records inconsistent with first-censor construction:",ambiguous,"\n")
  d
}

restrict_tamoxifen <- function(d,years) {
  stopifnot(length(years)==1L,is.finite(years),years>0)
  tau <- as.numeric(years)
  td <- d$td/365.25
  tr <- d$tr/365.25
  Delta1 <- as.integer(d$dd==1L & td<=tau)
  # The lower-priority endpoint is known until death or tau because of the
  # verified common observation horizon for all three relapse processes.
  R <- as.integer(Delta1==1L | td>=tau)
  obs_time <- pmin(td,tau)
  Delta2 <- as.integer(d$dr==1L & tr<obs_time+1e-12)
  T2tilde <- ifelse(Delta2==1L,tr,obs_time)
  out <- data.frame(d[,c("id","A",BASE_NAMES),drop=FALSE],
                    obs_time=obs_time,R=R,cens_event=1L-R,
                    Delta1=Delta1,Delta2=Delta2,
                    T1tilde=obs_time,T2tilde=T2tilde,
                    Yp=ifelse(R==1L,obs_time,NA_real_),
                    Ys=ifelse(R==1L,T2tilde,NA_real_))
  attr(out,"tau") <- tau
  if (any(out$Delta2==1L & out$T2tilde>=out$obs_time+1e-10) ||
      any(out$T2tilde>out$T1tilde) || any(out$obs_time<=0))
    stop("Inconsistent restricted outcomes")
  cat(sprintf("Horizon %.1f y: %d patients; %.1f%% censored; %.1f%% death; %.1f%% first relapse\n",
      tau,nrow(out),100*mean(out$R==0),100*mean(out$Delta1),
      100*mean(out$Delta2)))
  out
}

# ---------------------- Exponential nuisance models -----------------------
# The fitting and prediction interface uses a constant set of named numeric
# predictors. It does not construct any scores from the validation outcomes.

make_ridge_fit <- function(exposure,event,Z,kappa=RIDGE_KAPPA) {
  exposure <- as.numeric(exposure)
  event <- as.numeric(event)
  Z <- as.matrix(Z)
  if (!length(exposure) || nrow(Z)!=length(exposure) ||
      any(!is.finite(exposure)) || any(exposure<=0) ||
      any(!is.finite(Z)) || any(!event %in% c(0,1)) ||
      !is.finite(kappa) || kappa<=0)
    stop("Invalid hazard-fitting data")
  E <- sum(event)
  predictors <- colnames(Z)
  if (is.null(predictors) || anyDuplicated(predictors))
    stop("Ridge predictors must be uniquely named")
  if (E==0) {
    return(list(alpha=log(1e-8), beta=setNames(rep(0,ncol(Z)),predictors),
                fallback=TRUE))
  }
  zz <- log(exposure)
  sE <- as.numeric(crossprod(Z,event))
  tolerance <- 1e-8*max(1,E)
  evaluate <- function(beta,hessian=TRUE) {
    eta <- drop(Z%*%beta)+zz
    shift <- max(eta)
    ex <- exp(eta-shift)
    prob <- ex/sum(ex)
    logtotal <- shift+log(sum(ex))
    zm <- as.numeric(crossprod(Z,prob))
    value <- E*logtotal-sum(sE*beta)+.5*kappa*sum(beta^2)
    grad <- E*zm-sE+kappa*beta
    if (!hessian) return(list(value=value,gradient=grad,logtotal=logtotal))
    Zc <- sweep(Z,2L,zm,"-")
    H <- E*crossprod(Zc,Zc*as.numeric(prob)) + diag(kappa,ncol(Z))
    list(value=value,gradient=grad,hessian=H,logtotal=logtotal)
  }
  beta <- rep(0,ncol(Z))
  for (iter in seq_len(100L)) {
    cur <- evaluate(beta)
    if (max(abs(cur$gradient)) < tolerance) break
    step <- as.numeric(solve(cur$hessian,cur$gradient))
    direction <- sum(step*cur$gradient)
    if (!is.finite(direction) || direction<=0) stop("Invalid ridge direction")
    accepted <- FALSE
    scale <- 1.0
    for (bt in seq_len(50L)) {
      candidate <- beta-scale*step
      nxt <- evaluate(candidate,FALSE)
      if (is.finite(nxt$value) &&
          nxt$value <= cur$value-1e-4*scale*direction) {
        accepted <- TRUE
        break
      }
      scale <- scale/2
    }
    if (!accepted) {
      if (max(abs(cur$gradient)) < 1e-6*max(1,E)) break
      stop("Ridge line search did not converge")
    }
    beta <- candidate
  }
  last <- evaluate(beta,FALSE)
  if (max(abs(last$gradient))>1e-6*max(1,E))
    stop("Ridge solver has excessive score residual")
  alpha <- log(E)-last$logtotal
  if (!is.finite(alpha) || any(!is.finite(beta))) stop("Invalid ridge coefficients")
  list(alpha=as.numeric(alpha),beta=setNames(beta,predictors),fallback=FALSE)
}

semi_t2obs <- function(d) ifelse(d$Delta2==1L,d$T2tilde,Inf)

# Exponential hazards within the pre/post-recurrence states.
# The baseline time hazard is constant within each recurrence state.
# These are working models, not nonparametric baseline hazard estimators.
semi_ped <- function(d, event_type, predictors) {
  stopifnot(event_type %in% c("C","T1","T2"))
  if (event_type=="T2") {
    Z <- as.matrix(d[,predictors,drop=FALSE])
    return(list(exposure=pmax(d$T2tilde,1e-12),
                event=d$Delta2,Z=Z))
  }
  end <- d$obs_time
  event <- if(event_type=="C") d$cens_event else d$Delta1
  t2 <- semi_t2obs(d)
  pre <- pmin(end,t2)
  post <- pmax(end-t2,0)
  z0 <- d; z0$N2 <- 0
  z1 <- d; z1$N2 <- 1
  use0 <- which(pre>0)
  use1 <- which(post>0)
  if (!length(use0) && !length(use1)) stop("No at-risk exposure")
  X0 <- as.matrix(z0[use0,predictors,drop=FALSE])
  X1 <- as.matrix(z1[use1,predictors,drop=FALSE])
  list(exposure=c(pre[use0],post[use1]),
       event=as.numeric(c(event[use0]*(end[use0]<=t2[use0]),
                          event[use1]*(end[use1]>t2[use1]))),
       Z=rbind(X0,X1))
}

semi_fit_models <- function(train) {
  # Pooled censoring: allow treatment-specific intercept, shared covariate
  # and recurrence-history slopes. Arm-specific event-hazard regressions.
  cvars <- c("A",BASE_NAMES,"N2")
  cped <- semi_ped(train,"C",cvars)
  fitC <- make_ridge_fit(cped$exposure,cped$event,cped$Z)
  out <- vector("list",2L)
  for(a in 0:1) {
    da <- train[train$A==a,,drop=FALSE]
    if (nrow(da)<20L) stop("Too few training patients in an arm")
    t1vars <- c(BASE_NAMES,"N2")
    t1ped <- semi_ped(da,"T1",t1vars)
    t2ped <- semi_ped(da,"T2",BASE_NAMES)
    out[[a+1L]] <- list(C=fitC,
      T1=make_ridge_fit(t1ped$exposure,t1ped$event,t1ped$Z),
      T2=make_ridge_fit(t2ped$exposure,t2ped$event,t2ped$Z))
  }
  out
}

# Prediction for a single baseline row x (named vector), arm a and N2 status.
semi_design <- function(fit,x,a,n2) {
  z <- c(A=as.numeric(a),x,N2=as.numeric(n2))
  if (!all(names(fit$beta) %in% names(z))) stop("Missing hazard predictor")
  c("(Intercept)"=1, z[names(fit$beta)])
}

semi_rate <- function(fit,x,a,n2=0) {
  z <- semi_design(fit,x,a,n2)
  eta <- fit$alpha + sum(fit$beta*z[-1])
  if (!is.finite(eta) || eta>50 || eta< -700)
    stop("Extreme fitted log hazard: inspect censoring model")
  max(exp(eta),1e-12)
}

semi_G <- function(t,x,a,t2obs,fitC) {
  l0 <- semi_rate(fitC,x,a,0)
  l1 <- semi_rate(fitC,x,a,1)
  exp(-l0*pmin(t,t2obs)-l1*pmax(t-t2obs,0))
}

get_x <- function(d,i) {
  setNames(as.numeric(unlist(d[i,BASE_NAMES,drop=FALSE],use.names=FALSE)),BASE_NAMES)
}

# ------------------------ Weighted outcome projections ----------------------
# All calculations of h use strict ordering for death and recurrence times,
# with half credit for exact ties (including tied death times before tau).

semi_wcdf <- function(z,w,tol=1e-10) {
  ord <- order(z)
  z <- z[ord]; w <- w[ord]
  sums <- c(0,cumsum(w))
  less <- function(y) sums[findInterval(y-tol,z)+1L]
  leq  <- function(y) sums[findInterval(y+tol,z)+1L]
  list(less=less,leq=leq,total=sum(w))
}

semi_make_h <- function(F1,F0,tau,tol=1e-10) {
  F1p <- semi_wcdf(F1$Yp,F1$w,tol)
  F0p <- semi_wcdf(F0$Yp,F0$w,tol)
  jj1 <- abs(F1$Yp-tau)<=tol
  jj0 <- abs(F0$Yp-tau)<=tol
  F1s <- semi_wcdf(F1$Ys[jj1],F1$w[jj1],tol)
  F0s <- semi_wcdf(F0$Ys[jj0],F0$w[jj0],tol)

  h1 <- function(yp,ys) {
    out <- F0p$less(yp)
    mass <- F0p$leq(yp)-F0p$less(yp)
    ii <- which(mass > 1e-12)
    if (length(ii)) for(j in ii) {
      if (abs(yp[j]-tau)<=tol) {
        out[j] <- out[j]+F0s$less(ys[j])+
          .5*(F0s$leq(ys[j])-F0s$less(ys[j]))
      } else {
        tie <- which(abs(F0$Yp-yp[j])<=tol)
        out[j] <- out[j]+sum(F0$w[tie]*
           (as.numeric(ys[j]>F0$Ys[tie]+tol)+
             .5*as.numeric(abs(ys[j]-F0$Ys[tie])<=tol)))
      }
    }
    pmin(pmax(out,0),1)
  }
  h0 <- function(yp,ys) {
    out <- 1-F1p$leq(yp)
    mass <- F1p$leq(yp)-F1p$less(yp)
    ii <- which(mass > 1e-12)
    if (length(ii)) for(j in ii) {
      if (abs(yp[j]-tau)<=tol) {
        out[j] <- out[j]+F1s$total-F1s$leq(ys[j])+
          .5*(F1s$leq(ys[j])-F1s$less(ys[j]))
      } else {
        tie <- which(abs(F1$Yp-yp[j])<=tol)
        out[j] <- out[j]+sum(F1$w[tie]*
          (as.numeric(F1$Ys[tie]>ys[j]+tol)+
            .5*as.numeric(abs(F1$Ys[tie]-ys[j])<=tol)))
      }
    }
    pmin(pmax(out,0),1)
  }
  list(h1=h1,h0=h0)
}

semi_make_projection <- function(train,models,tau,g_floor=G_NUMERIC_FLOOR) {
  laws <- vector("list",2L)
  for(a in 0:1) {
    z <- train[train$A==a & train$R==1L,,drop=FALSE]
    if(nrow(z)<5L) stop("Too few completely revealed outcomes")
    g <- vapply(seq_len(nrow(z)),function(i) {
      x <- get_x(z,i)
      semi_G(z$obs_time[i],x,a,
             if(z$Delta2[i]==1L) z$T2tilde[i] else Inf,
             models[[a+1L]]$C)
    }, numeric(1))
    if(any(!is.finite(g)) || any(g<=g_floor))
      stop("Fitted G below numerical floor in outcome projection")
    w <- 1/g
    w <- w/sum(w)
    laws[[a+1L]] <- list(Yp=z$Yp,Ys=z$Ys,w=w)
  }
  hs <- semi_make_h(laws[[2L]],laws[[1L]],tau)
  psi <- sum(laws[[2L]]$w * hs$h1(laws[[2L]]$Yp,laws[[2L]]$Ys))
  rev_psi <- sum(laws[[1L]]$w * hs$h0(laws[[1L]]$Yp,laws[[1L]]$Ys))
  if (abs(psi-rev_psi)>1e-8)
    stop(sprintf("h1/h0 projection mismatch: %.9g vs %.9g",psi,rev_psi))
  list(psi=psi,h1=hs$h1,h0=hs$h0)
}

# ---------------- Conditional residual q and censoring augmentation ---------
# This is conditional expectation under the fitted event-history models.
# M_q points are used in a deterministic L x L conditional-quantile grid.

semi_q <- function(a,u,x,t2obs,fit,proj,mu,tau,M_q=M_Q) {
  if (u>tau) return(0)
  L <- max(4L,as.integer(ceiling(sqrt(M_q))))
  U <- (seq_len(L)-.5)/L
  rate2 <- semi_rate(fit$T2,x,a,0)
  rate10 <- semi_rate(fit$T1,x,a,0)
  rate11 <- semi_rate(fit$T1,x,a,1)
  if (is.finite(t2obs) && t2obs<u+1e-10) {
    e1 <- -log1p(-U)/rate11
    yp <- pmin(u+e1,tau)
    ys <- rep(pmin(t2obs,tau),L)
  } else {
    e2 <- rep(-log1p(-U)/rate2,each=L)
    E1 <- rep(-log1p(-U),times=L)
    t2 <- u+e2
    before <- E1<=rate10*e2
    t1 <- ifelse(before,u+E1/rate10,
                 t2+(E1-rate10*e2)/rate11)
    yp <- pmin(t1,tau)
    ys <- pmin(t2,t1,tau)
  }
  h <- if(a==1L) proj$h1(yp,ys) else proj$h0(yp,ys)
  mean(h)-mu
}

GLX <- c(-.9602898564975363,-.7966664774136267,
         -.5255324099163290,-.1834346424956498,
          .1834346424956498,.5255324099163290,
          .7966664774136267,.9602898564975363)
GLW <- c(.1012285362903763,.2223810344533745,
         .3137066458778873,.3626837833783620,
         .3626837833783620,.3137066458778873,
         .2223810344533745,.1012285362903763)

semi_augment <- function(a,x,stop,cens,t2obs,fit,proj,mu,tau,
                         M_q=M_Q,g_floor=G_NUMERIC_FLOOR) {
  if (stop<=0) return(0)
  q_eval <- function(u) semi_q(a,u,x,t2obs,fit,proj,mu,tau,M_q)
  G_eval <- function(u) {
    g <- semi_G(u,x,a,t2obs,fit$C)
    if (any(!is.finite(g)) || any(g<=g_floor))
      stop("Fitted G below numerical floor in augmentation")
    g
  }
  integrate_segment <- function(l,r,n2) {
    if(r<=l) return(0)
    uu <- (l+r)/2+(r-l)/2*GLX
    q <- vapply(uu,q_eval,numeric(1))
    lc <- semi_rate(fit$C,x,a,n2)
    (r-l)/2*sum(GLW*q/G_eval(uu)*lc)
  }
  cut <- min(stop,t2obs)
  drift <- integrate_segment(0,cut,0)
  if(is.finite(t2obs) && t2obs<stop)
    drift <- drift+integrate_segment(t2obs,stop,1)
  jump <- if(cens==1L) q_eval(stop)/G_eval(stop) else 0
  as.numeric(jump-drift)
}

# ------------------- IPCW direct SE: pooled censoring -----------------------

semi_censor_score_info <- function(d,fitC) {
  if(isTRUE(fitC$fallback))
    stop("No censoring events in fold training set: IPCW SE unavailable")
  n <- nrow(d)
  p <- length(fitC$beta)+1L
  score <- matrix(0,n,p)
  dH <- matrix(0,n,p)
  info <- matrix(0,p,p)
  for(i in seq_len(n)) {
    x <- get_x(d,i); a <- d$A[i]
    t2 <- if(d$Delta2[i]==1L) d$T2tilde[i] else Inf
    H0 <- semi_rate(fitC,x,a,0)*min(d$obs_time[i],t2)
    H1 <- semi_rate(fitC,x,a,1)*max(d$obs_time[i]-t2,0)
    z0 <- unname(semi_design(fitC,x,a,0))
    z1 <- unname(semi_design(fitC,x,a,1))
    zjump <- if(is.finite(t2) && t2<d$obs_time[i]) z1 else z0
    score[i,] <- d$cens_event[i]*zjump-H0*z0-H1*z1
    dH[i,] <- H0*z0+H1*z1
    info <- info + H0*tcrossprod(z0)+H1*tcrossprod(z1)
  }
  list(score=score,dH=dH,info=info)
}

semi_ipcw_inference <- function(d,w,folds,cfits,tau,kappa=RIDGE_KAPPA) {
  n <- nrow(d)
  i1 <- which(d$A==1L & d$R==1L & is.finite(w) & w>0)
  i0 <- which(d$A==0L & d$R==1L & is.finite(w) & w>0)
  if(!length(i1) || !length(i0)) stop("No revealed outcomes for IPCW")
  w1 <- w[i1]/sum(w[i1]); w0 <- w[i0]/sum(w[i0])
  F1 <- list(Yp=d$Yp[i1],Ys=d$Ys[i1],w=w1)
  F0 <- list(Yp=d$Yp[i0],Ys=d$Ys[i0],w=w0)
  hs <- semi_make_h(F1,F0,tau)
  h <- numeric(n)
  h[i1] <- hs$h1(d$Yp[i1],d$Ys[i1])
  h[i0] <- hs$h0(d$Yp[i0],d$Ys[i0])
  psi <- sum(w1*h[i1])
  if(abs(psi-sum(w0*h[i0]))>1e-8)
    stop("IPCW h1/h0 projection mismatch")
  wn <- numeric(n);wn[i1] <- w1;wn[i0] <- w0
  weighted_IF <- numeric(n)
  ii <- c(i1,i0)
  weighted_IF[ii] <- n*wn[ii]*(h[ii]-psi)
  nuisance_IF <- numeric(n)
  for(v in seq_len(max(folds))) {
    tr <- which(folds!=v)
    va <- which(folds==v & d$R==1L)
    fitC <- cfits[[v]]
    si <- semi_censor_score_info(d[tr,,drop=FALSE],fitC)
    p <- ncol(si$score)
    H <- si$info+diag(c(0,rep(kappa,p-1L)),p)
    if(rcond(H)<1e-10 || any(!is.finite(H)))
      stop("Censoring score information matrix near singular")
    D <- numeric(p)
    if(length(va)) {
      sv <- semi_censor_score_info(d[va,,drop=FALSE],fitC)
      D <- colSums(sv$dH*as.numeric(wn[va]*(h[va]-psi)))
    }
    coeff <- as.numeric(solve(H,D))
    nuisance_IF[tr] <- nuisance_IF[tr] +
      n*as.numeric(si$score%*%coeff)
  }
  phi <- weighted_IF+nuisance_IF
  se <- sqrt(mean((phi-mean(phi))^2)/n)
  list(psi=psi,se=se,phi=phi,
       weighted_IF=weighted_IF,nuisance_IF=nuisance_IF)
}

# ---------------------------- Cross-fitted estimator -------------------------
semi_one_step <- function(d,V=N_FOLDS,seed=FOLD_SEED,M_q=M_Q,
                          g_floor=G_NUMERIC_FLOOR,piA=.5) {
  n <- nrow(d); tau <- attr(d,"tau")
  if(is.null(tau)) stop("attr(d,'tau') is required")
  if(V<2 || V>=n || piA<=0 || piA>=1) stop("Invalid cross-fitting inputs")
  set.seed(seed)
  # Balance folds within treatment arms, preserving representation.
  folds <- integer(n)
  for(a in 0:1) {
    ix <- which(d$A==a)
    folds[ix] <- sample(rep(seq_len(V),length.out=length(ix)))
  }
  phi <- numeric(n)
  psi_pl <- numeric(V)
  w_ipcw <- numeric(n)
  cfits <- vector("list",V)
  fitted_G <- numeric(n)
  base_IF <- residual_IF <- augment_IF <- numeric(n)

  for(v in seq_len(V)) {
    train <- d[folds!=v,,drop=FALSE]
    ids <- which(folds==v)
    models <- semi_fit_models(train)
    proj <- semi_make_projection(train,models,tau,g_floor)
    psi_pl[v] <- proj$psi
    cfits[[v]] <- models[[1L]]$C
    cat(sprintf("  Fold %d: train %d, validation %d, plug-in %.5f\n",
                v,nrow(train),length(ids),psi_pl[v]))
    for(i in ids) {
      a <- d$A[i]; x <- get_x(d,i)
      mu1 <- semi_q(1L,0,x,Inf,models[[2L]],proj,0,tau,M_q)
      mu0 <- semi_q(0L,0,x,Inf,models[[1L]],proj,0,tau,M_q)
      mu1 <- min(max(mu1,0),1)
      mu0 <- min(max(mu0,0),1)
      mu_a <- if(a==1L) mu1 else mu0
      fit <- models[[a+1L]]
      t2obs <- if(d$Delta2[i]==1L) d$T2tilde[i] else Inf
      s <- d$obs_time[i]
      G <- semi_G(s,x,a,t2obs,fit$C)
      if(!is.finite(G) || G<=g_floor)
        stop(sprintf("Bad fitted G for patient id=%d fold=%d",d$id[i],v))
      fitted_G[i] <- G
      pa <- if(a==1L) piA else 1-piA
      rpart <- 0
      if(d$R[i]==1L) {
        h <- if(a==1L) proj$h1(d$Yp[i],d$Ys[i]) else
                            proj$h0(d$Yp[i],d$Ys[i])
        rpart <- (h-mu_a)/G
        w_ipcw[i] <- 1/(pa*G)
      }
      aug <- semi_augment(a,x,s,d$cens_event[i],t2obs,
                          fit,proj,mu_a,tau,M_q,g_floor)
      base_IF[i] <- mu1+mu0-2*psi_pl[v]
      residual_IF[i] <- rpart/pa
      augment_IF[i] <- aug/pa
      phi[i] <- base_IF[i]+residual_IF[i]+augment_IF[i]
    }
  }
  pl <- sum(tabulate(folds,nbins=V)*psi_pl)/n
  psi_full <- pl+mean(phi)
  se_full <- sqrt(mean((phi-mean(phi))^2)/n)
  ipc <- semi_ipcw_inference(d,w_ipcw,folds,cfits,tau)
  # Useful diagnostics of the influence of extreme fitted censoring weights.
  ord <- head(order(abs(phi),decreasing=TRUE),10L)
  largest <- data.frame(id=d$id[ord],fold=folds[ord],A=d$A[ord],
                        Ghat=fitted_G[ord],
                        R=d$R[ord],IF=phi[ord],
                        baseline=base_IF[ord],
                        residual=residual_IF[ord],
                        augmentation=augment_IF[ord])
  list(psi_full=psi_full,se_full=se_full,
       psi_ipcw=ipc$psi,se_ipcw=ipc$se,
       plugin=pl,phi_full=phi,phi_ipcw=ipc$phi,
       baseline_IF=base_IF,residual_IF=residual_IF,
       augmentation_IF=augment_IF,
       folds=folds,Ghat=fitted_G,
       min_G=min(fitted_G),
       quantiles_G=quantile(fitted_G,c(0,.01,.05,.1,.5,1)),
       largest_IF=largest,
       censoring_fallback=any(vapply(cfits,function(z)isTRUE(z$fallback),logical(1))))
}

# ---------------------------- Summaries and exports -------------------------
format_estimate <- function(method,psi,se,tau,ratio=NA_real_) {
  if(!is.finite(psi) || !is.finite(se) || se<=0)
    stop("Invalid estimate or standard error")
  if(psi<=0 || psi>=1)
    stop("Win-probability estimate outside (0,1); inspect influence weights")
  z <- qnorm(.975)
  theta <- psi/(1-psi)
  se_log_odds <- se/(psi*(1-psi))
  ci_odds <- exp(log(theta)+c(-1,1)*z*se_log_odds)
  ci_psi <- psi+c(-1,1)*z*se
  data.frame(horizon_years=tau,method=method,
    psi=psi,se_psi=se,psi_ci_low=ci_psi[1],psi_ci_high=ci_psi[2],
    win_odds=theta,win_odds_ci_low=ci_odds[1],win_odds_ci_high=ci_odds[2],
    log_odds_p=2*pnorm(-abs(log(theta)/se_log_odds)),
    estimated_RE=ratio)
}

run_tamoxifen_analysis <- function() {
  raw <- prepare_tamoxifen()
  output <- list()
  for(tau in RUN_HORIZONS) {
    cat("\n",paste(rep("=",65),collapse=""),"\n",sep="")
    d <- restrict_tamoxifen(raw,tau)
    # Same fold assignment at both horizons (patient order unchanged).
    ans <- semi_one_step(d,V=N_FOLDS,seed=FOLD_SEED,M_q=M_Q)
    re <- ans$se_ipcw^2/ans$se_full^2
    tab <- rbind(
      format_estimate("IPCW",ans$psi_ipcw,ans$se_ipcw,tau,ratio=1),
      format_estimate("One-step (full q)",ans$psi_full,ans$se_full,tau,ratio=re))
    cat("\nResults (RE is ratio of squared ESTIMATED SEs):\n")
    print(tab,digits=5,row.names=FALSE)
    cat("\nEstimated censoring survival summary:\n")
    print(ans$quantiles_G,digits=5)
    cat("Largest absolute one-step IF contributions:\n")
    print(ans$largest_IF,digits=5,row.names=FALSE)
    cat("Censoring fallback in any fold:",ans$censoring_fallback,"\n")
    suffix <- gsub("\\.","p",sprintf("%.1f",tau))
    write.csv(tab,paste0(OUTPUT_PREFIX,"_",suffix,"y_estimates.csv"),
              row.names=FALSE)
    write.csv(ans$largest_IF,
              paste0(OUTPUT_PREFIX,"_",suffix,"y_largest_IF.csv"),
              row.names=FALSE)
    output[[suffix]] <- list(table=tab,detail=ans)
  }
  all_tables <- do.call(rbind,lapply(output,`[[`,"table"))
  rownames(all_tables) <- NULL
  write.csv(all_tables,paste0(OUTPUT_PREFIX,"_all_results.csv"),
            row.names=FALSE)
  saveRDS(output,paste0(OUTPUT_PREFIX,"_full_output.rds"))
  cat("\nAll outputs saved.\n")
  invisible(output)
}

# ------------------------ Lightweight deterministic checks ------------------
tamoxifen_sanity <- function() {
  h <- function(p1,s1,p0,s0,tol=1e-10) {
    gt <- p1>p0+tol
    eq <- abs(p1-p0)<=tol
    out <- as.numeric(gt)
    out[eq] <- as.numeric(s1[eq]>s0[eq]+tol)+
                 .5*as.numeric(abs(s1[eq]-s0[eq])<=tol)
    out
  }
  a <- c(1,1,6,6,3); b <- c(.5,.5,1,6,2)
  c <- c(2,1,6,6,4); d <- c(1,.5,2,6,1)
  stopifnot(max(abs(h(a,b,c,d)+h(c,d,a,b)-1))<1e-12)
  raw <- prepare_tamoxifen()
  stopifnot(nrow(raw)==641L,
            sum(raw$dr)==76L,
            sum(raw$A==1L)==320L,
            sum(raw$A==0L)==321L)
  for (t in RUN_HORIZONS) {
    z <- restrict_tamoxifen(raw,t)
    stopifnot(all(z$R %in% 0:1),
              all(is.na(z$Yp[z$R==0])),
              all(z$T2tilde<=z$obs_time))
  }
  # counts from this version of the public dataset; stop if changed.
  z4 <- restrict_tamoxifen(raw,4)
  z5 <- restrict_tamoxifen(raw,5)
  stopifnot(sum(z4$R==0L)==110L,sum(z5$R==0L)==218L,
            sum(z4$Delta1)==36L,sum(z5$Delta1)==38L,
            sum(z4$Delta2)==54L,sum(z5$Delta2)==57L)
  cat("Data and win-kernel checks passed.\n")
  invisible(TRUE)
}

# Verify input data and event ordering before analysis.
tamoxifen_sanity()
results <- run_tamoxifen_analysis()
