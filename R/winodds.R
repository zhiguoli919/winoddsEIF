#' Estimate marginal win probability and win odds: scalar survival outcomes
#'
#' Internal fitting uses the same ridge-penalized exponential event model,
#' pooled A+X censoring model, and cross-fitted one-step correction as the
#' supplied scalar simulation code. Treatment probability is constant and known.
#'
#' @param data A data frame with columns A, X, time, status.
#' @param tau Positive analysis horizon in the units of time.
#' @param treatment_prob Known allocation probability P(A=1).
#' @param folds Number of cross-fitting folds (default 5).
#' @param seed Random seed for fold assignment.
#' @return A winoddsEIF object with one-step estimate, SE, confidence intervals,
#'   and an IPCW comparator.
#' @export
winprob_scalar <- function(data, tau, treatment_prob = 0.5,
                           folds = 5L, seed = 2027L) {
  .wo_validate_basic(data, tau, treatment_prob, folds, seed,
                     c("A", "X", "time", "status"))
  .wo_validate_binary(data$status, "status")
  .wo_check_time(data$time, "time")
  if (any(data$time < tau & data$status == 1L & data$time <= 0))
    stop("Event times must be positive.", call. = FALSE)
  n <- nrow(data)
  d <- data.frame(
    X = as.numeric(data$X), A = as.integer(data$A),
    s = pmin(as.numeric(data$time), tau),
    Delta = as.integer(data$status == 1L & data$time <= tau))
  d$R <- as.integer(d$Delta == 1L | data$time >= tau)
  .wo_check_fold_viability(d$A, d$R, folds)
  set.seed(as.integer(seed))
  fold <- sample(rep(seq_len(as.integer(folds)), length.out = n))
  # Only tau and the constant randomization probability are required by the
  # fitted estimator; simulation parameters are never passed into this call.
  p <- list(tau = tau, piA = treatment_prob)
  result <- .wo_scalar_estimate_once(d, p, fold = fold)
  .wo_make_result(
    psi = result$full, se = result$se_full,
    psi_plugin = result$plugin,
    psi_ipcw = result$ipcw, se_ipcw = result$se_ipcw,
    type = "scalar", tau = tau, n = n, folds = folds,
    diagnostics = list(min_fitted_G = result$min_G,
                       min_effective_n_treated = result$min_ess_treated,
                       min_effective_n_control = result$min_ess_control,
                       hazard_fallbacks = result$fallback,
                       ipcw_se_message = result$ipcw_error))
}

#' Estimate marginal win probability: prioritized semi-competing risks
#'
#' The terminal event is prioritized over the nonterminal event. Among pairs
#' tied on the restricted terminal event time, restricted nonterminal times
#' are compared. Censoring may depend on observed nonterminal-event history.
#'
#' @param data A data frame with columns A, X, time1, status1, time2, status2.
#'   time1 is min(terminal time, censoring time, tau); time2 is
#'   min(nonterminal time, terminal time, censoring time, tau).
#' @param tau Positive analysis horizon.
#' @param treatment_prob Known allocation probability P(A=1).
#' @param folds Number of cross-fitting folds (default 5).
#' @param seed Random seed for fold assignment.
#' @param m_q Number of deterministic integration grid points (rounded up to
#'   an L by L grid), default 256.
#' @param mu_method Either "model" (default; event-model integration) or
#'   "quadratic" (weighted working regression).
#' @param g_floor Numerical minimum acceptable fitted censoring survival.
#'   Below this value, estimation fails rather than truncating weights.
#' @return A winoddsEIF object.
#' @examples
#' \dontrun{
#' set.seed(2027)
#' n <- 240L
#' tau <- 5
#' A <- rbinom(n, 1, 0.5)
#' X <- runif(n, -1, 1)
#' # Terminal event has higher priority than nonterminal event.
#' T1 <- rexp(n, rate = 0.10 * exp(-0.35*A + 0.7*X))
#' T2 <- rexp(n, rate = 0.20 * exp(-0.30*A + 0.4*X))
#' C  <- rexp(n, rate = 0.08 * exp(0.2*X))
#' d <- data.frame(
#'   A = A, X = X,
#'   time1 = pmin(T1, C, tau),
#'   status1 = as.integer(T1 <= C & T1 <= tau),
#'   time2 = pmin(T2, T1, C, tau),
#'   status2 = as.integer(T2 < T1 & T2 <= C & T2 <= tau)
#' )
#' fit <- winprob_semicomp(d, tau = tau, treatment_prob = 0.5,
#'                         folds = 5, seed = 2027, m_q = 16)
#' print(fit)
#' fit$diagnostics
#' }
#' @export
winprob_semicomp <- function(data, tau, treatment_prob = 0.5,
                             folds = 5L, seed = 2027L, m_q = 256L,
                             mu_method = c("model", "quadratic"),
                             g_floor = 1e-8) {
  .wo_validate_basic(data, tau, treatment_prob, folds, seed,
                     c("A", "X", "time1", "status1", "time2", "status2"))
  .wo_validate_binary(data$status1, "status1")
  .wo_validate_binary(data$status2, "status2")
  .wo_check_time(data$time1, "time1")
  .wo_check_time(data$time2, "time2")
  if (any(data$time2 > data$time1 + 1e-8))
    stop("time2 cannot exceed time1.", call. = FALSE)
  if (any(data$status2 == 1L & data$time2 >= data$time1 - 1e-8))
    stop("A nonterminal event must strictly precede the observed terminal/censoring end.",
         call. = FALSE)
  if (!is.numeric(m_q) || length(m_q) != 1L || !is.finite(m_q) || m_q < 16)
    stop("m_q must be at least 16.", call. = FALSE)
  if (!is.numeric(g_floor) || length(g_floor) != 1L ||
      !is.finite(g_floor) || g_floor <= 0 || g_floor >= 1)
    stop("g_floor must lie strictly between 0 and 1.", call. = FALSE)
  mu_method <- match.arg(mu_method)
  time1 <- pmin(as.numeric(data$time1), tau)
  time2 <- pmin(as.numeric(data$time2), time1)
  delta1 <- as.integer(data$status1 == 1L & data$time1 <= tau)
  delta2 <- as.integer(data$status2 == 1L & data$time2 <= tau)
  R <- as.integer(delta1 == 1L | data$time1 >= tau)
  d <- data.frame(
    X = as.numeric(data$X), A = as.integer(data$A),
    obs_time = time1, R = R, cens_event = 1L - R,
    T1tilde = time1, Delta1 = delta1,
    T2tilde = time2, Delta2 = delta2,
    Yp = ifelse(R == 1L, time1, NA_real_),
    Ys = ifelse(R == 1L, time2, NA_real_))
  attr(d, "tau") <- tau
  .wo_check_fold_viability(d$A, d$R, folds)
  result <- semi_one_step(d, V = as.integer(folds), M_q = as.integer(m_q),
                          g_floor = g_floor, seed = as.integer(seed),
                          piA = treatment_prob, mu_method = mu_method,
                          oracle = FALSE)
  .wo_make_result(
    psi = result$full, se = result$se_full,
    psi_plugin = result$plugin,
    psi_ipcw = result$ipcw, se_ipcw = result$se_ipcw,
    type = "semi-competing risks", tau = tau, n = nrow(data), folds = folds,
    diagnostics = list(min_fitted_G = result$minG,
                       ipcw_se_message = result$ipcw_se_error,
                       mu_method = mu_method, m_q = as.integer(m_q)))
}

.wo_validate_basic <- function(data, tau, treatment_prob, folds, seed, required) {
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (!all(required %in% names(data)))
    stop("Missing required columns: ",
         paste(setdiff(required, names(data)), collapse = ", "), call. = FALSE)
  if (nrow(data) < 40L) stop("At least 40 observations are required.", call. = FALSE)
  if (anyNA(data[, required, drop = FALSE]))
    stop("Required columns cannot contain missing values.", call. = FALSE)
  .wo_validate_binary(data$A, "A")
  if (!is.numeric(data$X) || any(!is.finite(data$X)))
    stop("X must be a finite numeric baseline covariate.", call. = FALSE)
  if (length(tau) != 1L || !is.numeric(tau) || !is.finite(tau) || tau <= 0)
    stop("tau must be a positive finite number.", call. = FALSE)
  if (length(treatment_prob) != 1L || !is.numeric(treatment_prob) ||
      !is.finite(treatment_prob) || treatment_prob <= 0 || treatment_prob >= 1)
    stop("treatment_prob must be between 0 and 1.", call. = FALSE)
  if (length(folds) != 1L || !is.numeric(folds) || !is.finite(folds) ||
      folds != as.integer(folds) || folds < 2 || folds > nrow(data) / 4)
    stop("folds must be an integer from 2 to nrow(data)/4.", call. = FALSE)
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed < 0 || seed != as.integer(seed))
    stop("seed must be a nonnegative integer.", call. = FALSE)
  invisible(TRUE)
}
.wo_validate_binary <- function(x, nm) {
  if (!is.numeric(x) && !is.integer(x))
    stop(nm, " must be numeric 0/1.", call. = FALSE)
  if (any(!is.finite(x)) || any(!(x %in% c(0, 1))))
    stop(nm, " must consist of 0/1 values.", call. = FALSE)
  if (nm == "A" && length(unique(x)) != 2L)
    stop("Both treatment groups must be represented.", call. = FALSE)
  invisible(TRUE)
}
.wo_check_time <- function(x, nm) {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x <= 0))
    stop(nm, " must contain positive finite times.", call. = FALSE)
  invisible(TRUE)
}
.wo_check_fold_viability <- function(A, R, folds) {
  for (a in 0:1) {
    if (sum(A == a) < 2 * folds)
      stop("Each arm needs at least 2*folds subjects for reliable fitting.", call. = FALSE)
    if (sum(A == a & R == 1L) < 2 * folds)
      stop("Each arm needs at least 2*folds revealed outcomes.", call. = FALSE)
  }
  invisible(TRUE)
}
.wo_make_result <- function(psi, se, psi_plugin, psi_ipcw, se_ipcw,
                            type, tau, n, folds, diagnostics) {
  if (!is.finite(psi) || !is.finite(se) || se <= 0)
    stop("The fitted one-step estimate or SE is not finite.", call. = FALSE)
  ci_psi <- as.numeric(psi + c(-1, 1) * qnorm(.975) * se)
  names(ci_psi) <- c("lower", "upper")
  valid_odds <- psi > 0 && psi < 1
  odds <- if (valid_odds) psi/(1-psi) else NA_real_
  ci_odds <- if (valid_odds) {
    delta <- qnorm(.975) * se / (psi*(1-psi))
    exp(log(odds)+c(-delta,delta))
  } else c(NA_real_,NA_real_)
  names(ci_odds) <- c("lower", "upper")
  pvalue <- if(valid_odds) 2*pnorm(-abs(log(odds)*psi*(1-psi)/se)) else NA_real_
  structure(list(
    setting = type, n = as.integer(n), tau = tau, folds = as.integer(folds),
    psi = psi, se_psi = se, ci_psi = ci_psi,
    win_odds = odds, ci_win_odds = ci_odds, pvalue_win_odds = pvalue,
    plugin = psi_plugin, ipcw = psi_ipcw, se_ipcw = se_ipcw,
    diagnostics = diagnostics, ridge_kappa = RIDGE_KAPPA),
    class = "winoddsEIF")
}

#' @export
print.winoddsEIF <- function(x, ...) {
  cat("Cross-fitted win probability (", x$setting, ")\n", sep = "")
  cat("n =", x$n, "; tau =", x$tau, "; folds =", x$folds, "\n")
  cat("One-step psi:", signif(x$psi, 6), "; SE:", signif(x$se_psi, 5), "\n")
  cat("95% CI for psi:", paste(signif(x$ci_psi, 5), collapse = " to "), "\n")
  cat("Win odds:", signif(x$win_odds, 6), "; 95% CI:",
      paste(signif(x$ci_win_odds, 5), collapse = " to "), "\n")
  cat("Normalized IPCW psi:", signif(x$ipcw, 6), "\n")
  invisible(x)
}
