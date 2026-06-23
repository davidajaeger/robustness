#' Range tests for equality and equivalence across specifications
#'
#' Computes range-based and Wald-based tests of equality and equivalence
#' across a set of specifications, from saved bootstrap draws of the
#' coefficient of interest. This is the method of Jaeger (2026), and the R
#' counterpart of the Stata \code{robustness} command: given the same
#' bootstrap draws, the two produce the same statistics (up to Monte Carlo
#' error when the draws are generated separately).
#'
#' @details
#' For \eqn{K} specifications with full-sample estimates and \eqn{B} bootstrap
#' draws of each, \code{robustness()} computes, for every requested comparison,
#' two readings of the \emph{same} bootstrap draws.
#'
#' \emph{Equivalence.} The minimum equivalence bound
#' \eqn{R^*_{1-\alpha}} is the \eqn{(1-\alpha)} quantile of the
#' \emph{uncentred} bootstrap range, taken as the inverse of the empirical
#' c.d.f.: the smallest \eqn{d} such that at least a \eqn{(1-\alpha)} share of
#' bootstrap ranges are \eqn{\leq d}. This is the type-1 sample quantile and
#' is exactly the definition in the paper. \eqn{W^*_{1-\alpha}} is the
#' analogous bound on the Wald (Mahalanobis) scale, the square root of the
#' \eqn{(1-\alpha)} quantile of the uncentred bootstrap Wald statistic.
#'
#' \emph{Equality.} The p-values \eqn{p_R} and \eqn{p_W} are the shares of the
#' \emph{recentred} bootstrap statistics at or above the observed statistic.
#' Recentring subtracts each estimate's deviation from the cross-specification
#' mean, imposing a common probability limit (\eqn{\Delta = 0}).
#'
#' The two readings come from the same resamples and differ only in the
#' recentring. The bootstrap draws must be uncentred (raw resample estimates);
#' all recentring is internal.
#'
#' \code{robustness()} consumes bootstrap draws; it does not generate them. The
#' same shared resample must be used for all specifications on each
#' replication. Drawing specifications independently silently destroys the
#' joint distribution and gives wrong p-values and bounds.
#'
#' @param theta Numeric vector of \eqn{K} full-sample point estimates.
#' @param draws Numeric matrix of uncentred bootstrap draws, \eqn{B} rows by
#'   \eqn{K} columns, column \eqn{k} matching \code{theta[k]}.
#' @param comparisons Either a single integer vector of column indices (one
#'   comparison), or a named list of integer vectors (several comparisons, the
#'   names labelling them). Defaults to all specifications as one comparison
#'   named \code{"all"}.
#' @param alpha Numeric vector of significance levels for the equivalence
#'   bounds. Defaults to \code{c(0.50, 0.05)}, giving \eqn{R^*_{.50}} (the
#'   median of the bootstrap range) and \eqn{R^*_{.95}} (the bound).
#' @param n Optional per-specification sample sizes: a \eqn{B} by \eqn{K}
#'   matrix or a length-\eqn{K} vector. Reported descriptively; does not enter
#'   any statistic.
#' @param max_drop Maximum proportion of incomplete bootstrap replications
#'   tolerated before the function stops. Defaults to 0.01.
#' @param keep_draws Logical. If \code{TRUE}, the per-replication bootstrap
#'   series (the uncentred and recentred range and Wald) are retained on each
#'   comparison's result and can be extracted with \code{\link{bootstrap_draws}}
#'   for plotting. Defaults to \code{FALSE}, since these are \eqn{B}-length
#'   vectors per comparison.
#'
#' @return An object of class \code{"robustness"}: a list of per-comparison
#'   results (each of class \code{"range_test"}), with \code{print},
#'   \code{summary}, and \code{as.data.frame} methods. With
#'   \code{keep_draws = TRUE}, \code{\link{bootstrap_draws}} returns the
#'   per-replication series in long form.
#'
#' @references
#' Jaeger, David A. (2026). Robustness? Range Tests for Equality and
#' Equivalence Across Specifications.
#'
#' @examples
#' set.seed(1)
#' theta <- c(0.20, 0.22, 0.19, 0.21)
#' draws <- sapply(theta, function(m) rnorm(9999, m, 0.03))
#' robustness(theta, draws)
#' robustness(theta, draws,
#'            comparisons = list(all = 1:4, first_two = 1:2, extremes = c(1, 4)))
#'
#' # Retain the bootstrap series and plot the range distribution:
#' r <- robustness(theta, draws, keep_draws = TRUE)
#' d <- bootstrap_draws(r)
#' hist(d$range_unc[d$comparison == "all"])
#'
#' @export
robustness <- function(theta, draws, comparisons = NULL,
                       alpha = c(0.50, 0.05), n = NULL, max_drop = 0.01,
                       keep_draws = FALSE) {

  theta <- as.numeric(theta)
  draws <- as.matrix(draws)
  if (ncol(draws) != length(theta)) {
    stop("draws must have one column per element of theta (",
         ncol(draws), " columns vs ", length(theta), " estimates).")
  }

  if (is.null(comparisons)) {
    comparisons <- list(all = seq_along(theta))
  } else if (!is.list(comparisons)) {
    comparisons <- list(comparison = as.integer(comparisons))
  }
  if (is.null(names(comparisons)) || any(names(comparisons) == "")) {
    nm <- names(comparisons)
    if (is.null(nm)) nm <- rep("", length(comparisons))
    blank <- nm == ""
    nm[blank] <- paste0("comparison", seq_along(nm))[blank]
    names(comparisons) <- nm
  }

  results <- lapply(names(comparisons), function(cn) {
    .range_test_one(theta, draws, cols = comparisons[[cn]],
                    alpha = alpha, n = n, max_drop = max_drop, label = cn,
                    keep_draws = keep_draws)
  })
  names(results) <- names(comparisons)

  out <- list(results = results, alpha = alpha,
              K_total = length(theta), B = nrow(draws))
  class(out) <- "robustness"
  out
}


# Internal: statistics for one comparison. The computational core.
#
# Inputs:  theta  full-sample estimates (length >= K)
#          draws  uncentred bootstrap draws, B x (>= K)
#          cols   the specifications in this comparison (column indices)
#
# Output (class "range_test"):
#          R, W       observed range and Wald statistic, from the estimates
#          p_R, p_W   equality p-values, from the recentred bootstrap
#          equivalence  data frame: alpha, Rstar, Wstar (one row per alpha)
#
# The same uncentred draws are read two ways. Uncentred -> equivalence bounds
# (Rstar, Wstar); recentred -> equality p-values (p_R, p_W). Nothing else
# differs between the two.
.range_test_one <- function(theta, draws, cols = NULL, alpha = c(0.50, 0.05),
                            n = NULL, max_drop = 0.01, label = NULL,
                            keep_draws = FALSE) {

  if (is.null(cols)) cols <- seq_along(theta)
  cols <- as.integer(cols)
  if (anyDuplicated(cols))                         # a spec may appear only once
    stop("Comparison '", label, "' lists a specification more than once.")
  if (any(cols < 1L | cols > length(theta)))
    stop("Comparison '", label, "' references columns out of range ",
         "(valid 1 to ", length(theta), ").")

  th <- theta[cols]
  D  <- draws[, cols, drop = FALSE]
  K  <- length(cols)
  if (K < 2L)
    stop("Comparison '", label, "' needs at least 2 specifications, got ", K, ".")

  # Drop replications with any missing draw among these specs; count them and
  # stop if the incomplete share exceeds max_drop.
  B_orig    <- nrow(D)
  ok        <- stats::complete.cases(D)
  D         <- D[ok, , drop = FALSE]
  B         <- nrow(D)
  B_dropped <- B_orig - B
  if (B < K)
    stop(sprintf("Comparison '%s': only %d complete replications (of %d), need %d.",
                 label, B, B_orig, K))
  if (B_dropped / B_orig > max_drop)
    stop(sprintf(paste0("Comparison '%s': %d of %d replications incomplete ",
                        "(%.1f%%), exceeding max_drop = %.1f%%."),
                 label, B_dropped, B_orig, 100 * B_dropped / B_orig, 100 * max_drop))

  # Grand-mean contrast, (K-1) x K: row i is e_i minus the mean over all K
  # coordinates. The K estimates are equal iff all K-1 contrasts are zero.
  Rmat <- diag(K)[-K, , drop = FALSE] - 1 / K
  RVR  <- Rmat %*% stats::var(D) %*% t(Rmat)

  # Range statistics never use the contrast covariance, so they are always
  # defined.
  R_obs <- max(th) - min(th)
  R_unc <- apply(D, 1L, function(r) max(r) - min(r))
  # Recentred draws: subtract each spec's deviation from the cross-spec mean,
  # imposing Delta = 0. The common shift cancels in the range and in the
  # contrast, so it never affects the p-value; written to match the paper.
  Dc    <- sweep(D, 2L, th - mean(th), "-")
  R_rc  <- apply(Dc, 1L, function(r) max(r) - min(r))
  # Monte Carlo p-value (1 + #)/(B + 1): the observed statistic joins its own
  # reference set, so it is bounded away from zero (minimum 1/(B+1)) and uniform
  # under the null by exchangeability (Davison and Hinkley 1997). B is the
  # number of complete replications.
  p_R   <- (1 + sum(R_rc >= R_obs)) / (B + 1)

  # Wald statistics require a full-rank contrast covariance. Duplicate columns
  # are rejected above, so a rank-deficient RVR means the specifications are
  # genuinely collinear in the draws. Detect with a relative eigenvalue
  # tolerance rather than a rank() call: a default rank tolerance can miss
  # structured collinearity (e.g. one spec a constant shift of another) whose
  # true zero eigenvalue rounding lifts to a tiny positive value. Rather than
  # invert with a generalized inverse (a degenerate Wald), return W, p_W, and
  # W* as NA and keep the range results.
  ev      <- eigen(RVR, symmetric = TRUE, only.values = TRUE)$values
  wald_ok <- max(ev) > 0 && min(ev) > 1e-12 * max(ev)
  if (wald_ok) {
    RVRinv <- solve(RVR)
    d_obs  <- Rmat %*% th
    W_obs  <- as.numeric(t(d_obs) %*% RVRinv %*% d_obs)
    Bd     <- Rmat %*% t(D)
    W_unc  <- colSums(Bd * (RVRinv %*% Bd))
    Bdc    <- Rmat %*% t(Dc)
    W_rc   <- colSums(Bdc * (RVRinv %*% Bdc))
    p_W    <- (1 + sum(W_rc >= W_obs)) / (B + 1)
  } else {
    warning(sprintf(paste0("Comparison '%s': contrast covariance is rank ",
                           "deficient (collinear specifications); Wald ",
                           "statistics (W, p_W, W*) returned as NA. Range ",
                           "statistics are unaffected."), label))
    W_obs <- NA_real_
    p_W   <- NA_real_
    W_unc <- rep(NA_real_, B)
    W_rc  <- rep(NA_real_, B)
  }

  # Equivalence bounds: (1-alpha) quantile of the uncentred distribution, type 1
  # (inverse empirical c.d.f.). Rstar is always defined; Wstar is NA when the
  # Wald is undefined.
  eq <- data.frame(alpha = alpha, Rstar = NA_real_, Wstar = NA_real_)
  for (a in seq_along(alpha)) {
    eq$Rstar[a] <- stats::quantile(R_unc, 1 - alpha[a], type = 1, names = FALSE)
    if (!anyNA(W_unc))
      eq$Wstar[a] <- sqrt(stats::quantile(W_unc, 1 - alpha[a], type = 1, names = FALSE))
  }

  avg_n <- NULL
  if (!is.null(n)) {
    if (is.matrix(n)) avg_n <- colMeans(n[ok, cols, drop = FALSE], na.rm = TRUE)
    else              avg_n <- as.numeric(n)[cols]
  }

  out <- list(label = label, K = K, B = B, B_dropped = B_dropped,
              R = R_obs, W = W_obs, p_R = p_R, p_W = p_W,
              equivalence = eq, avg_n = avg_n)

  # The four bootstrap series, retained only on request (each is B-length).
  # These are the distributions the summaries above collapse: the (1-alpha)
  # quantile of range_unc is Rstar, the share of range_rc at or above R_obs is
  # p_R. Column names match the Stata command's saving() output.
  if (keep_draws) {
    out$draws <- data.frame(draw = seq_len(B),
                            range_unc = R_unc, range_rc = R_rc,
                            wald_unc = W_unc, wald_rc = W_rc)
  }
  class(out) <- "range_test"
  out
}
