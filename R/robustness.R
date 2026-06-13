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
#' draws of each, \code{robustness()} computes, for every requested comparison:
#' the observed range and Wald statistic; equivalence bounds
#' \eqn{\delta^*_R(\alpha)} and \eqn{\delta^*_W(\alpha)} from the uncentred
#' bootstrap distribution (\eqn{\delta^*_R(\alpha)} is the observed range plus
#' the \eqn{(1-\alpha)} quantile of the centred bootstrap range, and
#' \eqn{\delta^*_R(.50)} its median); and equality-test p-values \eqn{p_R} and
#' \eqn{p_W} from the recentred bootstrap distribution.
#'
#' The bootstrap draws must be uncentred (raw resample estimates). All
#' recentring is internal. Replications with any missing draw among the
#' selected specifications are dropped and counted; the function stops if more
#' than \code{max_drop} are incomplete.
#'
#' \code{robustness()} consumes bootstrap draws; it does not generate them. The
#' same shared resample must be used for all specifications on each
#' replication. Drawing specifications independently silently destroys the
#' joint distribution and gives wrong p-values and bounds. See the vignette for
#' a worked generation example.
#'
#' @param theta Numeric vector of \eqn{K} full-sample point estimates.
#' @param draws Numeric matrix of uncentred bootstrap draws, \eqn{B} rows by
#'   \eqn{K} columns, column \eqn{k} matching \code{theta[k]}.
#' @param comparisons Either a single integer vector of column indices (one
#'   comparison), or a named list of integer vectors (several comparisons, the
#'   names labelling them). Defaults to all specifications as one comparison
#'   named \code{"all"}.
#' @param alpha Numeric vector of significance levels for equivalence bounds.
#'   Defaults to \code{c(0.50, 0.05)}.
#' @param n Optional per-specification sample sizes: a \eqn{B} by \eqn{K}
#'   matrix or a length-\eqn{K} vector. Reported descriptively; does not enter
#'   any statistic.
#' @param max_drop Maximum proportion of incomplete bootstrap replications
#'   tolerated before the function stops. Defaults to 0.01.
#'
#' @return An object of class \code{"robustness"}: a list of per-comparison
#'   results (each of class \code{"range_test"}), with \code{print},
#'   \code{summary}, and \code{as.data.frame} methods.
#'
#' @references
#' Jaeger, David A. (2026). Robustness? Range Tests for Equality and
#' Equivalence Across Specifications. To cite this method, please cite that
#' paper; see \code{citation("robustness")}.
#'
#' @examples
#' set.seed(1)
#' theta <- c(0.20, 0.22, 0.19, 0.21)
#' draws <- sapply(theta, function(m) rnorm(999, m, 0.03))
#' robustness(theta, draws)
#' robustness(theta, draws,
#'            comparisons = list(all = 1:4, first_two = 1:2, extremes = c(1, 4)))
#'
#' @export
robustness <- function(theta, draws, comparisons = NULL,
                       alpha = c(0.50, 0.05), n = NULL, max_drop = 0.01) {

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
                    alpha = alpha, n = n, max_drop = max_drop, label = cn)
  })
  names(results) <- names(comparisons)

  out <- list(results = results, alpha = alpha,
              K_total = length(theta), B = nrow(draws))
  class(out) <- "robustness"
  out
}


# Internal: statistics for one comparison. The verified computational core.
.range_test_one <- function(theta, draws, cols = NULL, alpha = c(0.50, 0.05),
                            n = NULL, max_drop = 0.01, label = NULL) {

  if (is.null(cols)) cols <- seq_along(theta)
  cols <- as.integer(cols)
  if (any(cols < 1L | cols > length(theta))) {
    stop("Comparison '", label, "' references columns out of range ",
         "(valid 1 to ", length(theta), ").")
  }

  th <- theta[cols]
  D  <- draws[, cols, drop = FALSE]
  K  <- length(cols)
  if (K < 2L) {
    stop("Comparison '", label, "' needs at least 2 specifications, got ", K, ".")
  }

  B_orig    <- nrow(D)
  ok        <- stats::complete.cases(D)
  D         <- D[ok, , drop = FALSE]
  B         <- nrow(D)
  B_dropped <- B_orig - B
  if (B < K) {
    stop(sprintf("Comparison '%s': only %d complete replications (of %d), need %d.",
                 label, B, B_orig, K))
  }
  if (B_dropped / B_orig > max_drop) {
    stop(sprintf(paste0("Comparison '%s': %d of %d replications incomplete ",
                        "(%.1f%%), exceeding max_drop = %.1f%%."),
                 label, B_dropped, B_orig,
                 100 * B_dropped / B_orig, 100 * max_drop))
  }

  Rmat <- matrix(-1 / K, nrow = K - 1L, ncol = K)
  for (i in seq_len(K - 1L)) Rmat[i, i] <- Rmat[i, i] + 1

  Vhat   <- stats::var(D)
  RVR    <- Rmat %*% Vhat %*% t(Rmat)
  RVRinv <- solve(RVR)

  d_obs <- Rmat %*% th
  W_obs <- as.numeric(t(d_obs) %*% RVRinv %*% d_obs)
  R_obs <- max(th) - min(th)

  Db    <- Rmat %*% t(D)
  W_unc <- colSums(Db * (RVRinv %*% Db))
  R_unc <- apply(D, 1L, function(r) max(r) - min(r))

  thbar <- mean(th)
  Dc    <- sweep(D, 2L, th, "-") + thbar
  Dbc   <- Rmat %*% t(Dc)
  W_rc  <- colSums(Dbc * (RVRinv %*% Dbc))
  R_rc  <- apply(Dc, 1L, function(r) max(r) - min(r))

  p_W <- mean(W_rc >= W_obs)
  p_R <- mean(R_rc >= R_obs)

  eq <- data.frame(alpha = alpha, delta_R = NA_real_, delta_W = NA_real_)
  for (a in seq_along(alpha)) {
    q_R <- stats::quantile(R_unc - R_obs, probs = 1 - alpha[a],
                           names = FALSE, type = 7)
    q_W <- stats::quantile(W_unc - W_obs, probs = 1 - alpha[a],
                           names = FALSE, type = 7)
    eq$delta_R[a] <- R_obs + q_R
    eq$delta_W[a] <- sqrt(max(W_obs + q_W, 0))
  }

  avg_n <- NULL
  if (!is.null(n)) {
    if (is.matrix(n)) {
      nc    <- n[ok, cols, drop = FALSE]
      avg_n <- colMeans(nc, na.rm = TRUE)
    } else {
      avg_n <- as.numeric(n)[cols]
    }
  }

  out <- list(label = label, K = K, B = B, B_dropped = B_dropped,
              R = R_obs, W = W_obs, p_R = p_R, p_W = p_W,
              equivalence = eq, avg_n = avg_n)
  class(out) <- "range_test"
  out
}
