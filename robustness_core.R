#!/usr/bin/env Rscript
#===========================================================================
#  robustness_core.R
#  Shared engine for the Monte Carlo simulations in "Robustness? Range
#  Tests for Equality and Equivalence Across Specifications" -- Jaeger (2026)
#
#  Single source of truth. Every per-design runner sources this file, so
#  the quantile convention and the test statistics are defined exactly once.
#  Sourced by: run_D1.R ... run_D7.R, and the D8 module (which reuses
#  contrast_matrix, safe_quantile, and robustness_stats).
#
#  This file is design-independent. It contains no run parameters
#  (n_obs, MC, seed, n_cores, output paths) and no DGPs. Those live in the
#  runners and in dgps.R respectively.
#===========================================================================

library(parallel)
library(MASS)  # for ginv if needed

#===========================================================================
# Convention constants (defined once here, not in any runner)
#===========================================================================

# B is chosen so that (1 - alpha)(B + 1) = 9,500 is an integer, so the
# (1 - alpha) equivalence bound lands on an exact order statistic and
# safe_quantile (type = 1) returns it without interpolation. See eq. (529).
B     <- 9999   # bootstrap replications
alpha <- 0.05   # significance level

#===========================================================================
# Core functions (extracted verbatim from robustness_simulations.R)
#===========================================================================

#---------------------------------------------------------------------------
# Contrast matrix: grand-mean contrasts (K-1 x K)
#---------------------------------------------------------------------------
contrast_matrix <- function(K) {
  R <- matrix(-1/K, nrow = K-1, ncol = K)
  for (j in 1:(K-1)) R[j, j] <- R[j, j] + 1
  R
}

#---------------------------------------------------------------------------
# Quantile helper (handles NAs)
#---------------------------------------------------------------------------
safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  quantile(x, p, type = 1)  # type=1: order-statistic quantile (eq. 529 convention), no interpolation [was type=7]
}

#---------------------------------------------------------------------------
# robustness_stats: compute all test statistics from bootstrap draws
#
# theta:  K-vector of full-sample estimates
# DRAWS:  B x K matrix of bootstrap coefficient draws
# alpha:  significance level
#
# Returns a list with all statistics
#---------------------------------------------------------------------------
robustness_stats <- function(theta, DRAWS, alpha = 0.05) {

  K <- length(theta)
  B <- nrow(DRAWS)

  # Drop rows with any NA
  ok <- complete.cases(DRAWS)
  DRAWS <- DRAWS[ok, , drop = FALSE]
  B <- nrow(DRAWS)
  if (B < K) return(NULL)

  # Bootstrap V-hat
  Vhat <- cov(DRAWS)

  # Contrast matrix
  R <- contrast_matrix(K)
  RVR <- R %*% Vhat %*% t(R)

  # Try to invert; use ginv if singular
  RVR_inv <- tryCatch(solve(RVR), error = function(e) ginv(RVR))

  # Observed statistics
  d_obs <- R %*% theta
  W_obs <- as.numeric(t(d_obs) %*% RVR_inv %*% d_obs)
  R_obs <- max(theta) - min(theta)

  # ── Equivalence: uncentred bootstrap ──────────────────────────────
  W_boot <- numeric(B)
  R_boot <- numeric(B)
  for (b in 1:B) {
    d_b <- R %*% DRAWS[b, ]
    W_boot[b] <- as.numeric(t(d_b) %*% RVR_inv %*% d_b)
    R_boot[b] <- max(DRAWS[b, ]) - min(DRAWS[b, ])
  }

  q_W <- safe_quantile(W_boot - W_obs, 1 - alpha)
  delta_W <- sqrt(max(W_obs + q_W, 0))
  delta_R_05 <- safe_quantile(R_boot, 0.95)
  delta_R_10 <- safe_quantile(R_boot, 0.90)
  delta_R_50 <- safe_quantile(R_boot, 0.50)

  # ── Equality: recentred bootstrap ─────────────────────────────────
  theta_bar <- mean(theta)
  W_boot_rc <- numeric(B)
  R_boot_rc <- numeric(B)
  for (b in 1:B) {
    theta_rc <- DRAWS[b, ] - theta + theta_bar
    d_rc <- R %*% theta_rc
    W_boot_rc[b] <- as.numeric(t(d_rc) %*% RVR_inv %*% d_rc)
    R_boot_rc[b] <- max(theta_rc) - min(theta_rc)
  }

  p_W <- mean(W_boot_rc >= W_obs)
  p_R <- mean(R_boot_rc >= R_obs)

  # ── Analytical Wald: chi2(K-1) reference distribution ─────────────
  p_W_chi2 <- pchisq(W_obs, df = K - 1, lower.tail = FALSE)

  # ── Eyeballing criteria (on original estimates) ───────────────────
  bs_se <- sqrt(diag(Vhat))
  z_crit <- qnorm(1 - alpha/2)

  # Criterion 1: all same sign and individually significant
  all_same_sign <- all(theta > 0) | all(theta < 0)
  all_indiv_sig <- all(abs(theta) / bs_se > z_crit)
  eyeball_sign_sig <- all_same_sign & all_indiv_sig

  # Criterion 2: all pairwise CIs overlap
  ci_lo <- theta - z_crit * bs_se
  ci_hi <- theta + z_crit * bs_se
  all_overlap <- TRUE
  for (j in 1:(K-1)) {
    for (l in (j+1):K) {
      if (ci_hi[j] < ci_lo[l] || ci_hi[l] < ci_lo[j]) {
        all_overlap <- FALSE
        break
      }
    }
    if (!all_overlap) break
  }
  eyeball_ci_overlap <- all_overlap

  list(
    K = K, B = B,
    W_obs = W_obs, R_obs = R_obs,
    delta_W = delta_W, delta_R_05 = delta_R_05, delta_R_10 = delta_R_10, delta_R_50 = delta_R_50,
    p_W = p_W, p_R = p_R, p_W_chi2 = p_W_chi2,
    eyeball_sign_sig = eyeball_sign_sig,
    eyeball_ci_overlap = eyeball_ci_overlap,
    range_rejects = (p_R < alpha)
  )
}


#===========================================================================
# Estimation: OLS coefficient on D (first column of spec matrix)
#===========================================================================

estimate_specs <- function(dgp) {
  K <- dgp$K
  theta <- numeric(K)
  multi_Y <- isTRUE(dgp$multi_Y)

  for (k in 1:K) {
    Xmat <- dgp$specs[[k]]
    if (multi_Y) {
      Y <- dgp$Y_list[[k]]
    } else {
      Y <- dgp$Y
    }
    # Add intercept
    Xmat_int <- cbind(1, Xmat)
    fit <- lm.fit(Xmat_int, Y)
    # D is second column (after intercept)
    theta[k] <- fit$coefficients[2]
  }
  theta
}

#===========================================================================
# Population estimands: run the DGP once at very large n to obtain the
# plim of each spec's coefficient on D.  Used to validate the
# upper-bound coverage P(Delta <= delta*_R(p)) >= 1 - p in the simulation.
# A fixed seed (independent of the MC seed) gives a deterministic
# "truth" against which empirical coverage is computed.
#===========================================================================

compute_pop_estimands <- function(dgp_fn, n_pop = 500000L, seed_pop = 99999L) {
  set.seed(seed_pop)
  dgp <- dgp_fn(n_pop)
  estimate_specs(dgp)
}

#---------------------------------------------------------------------------
# Analytical Delta for each design (belt-and-braces).
#
# These are the population values of Delta = max_k tau_k - min_k tau_k
# derived analytically from the DGP.  The simulations use these
# hard-coded values; verify_Delta() runs a large-N Monte Carlo as a
# guard against drift between the analytical derivation and the
# implemented DGP.
#
# Derivations (all by FWL on the population regression):
#   D1: D indep X, all five specs unbiased for tau = 0.5.       Delta = 0.
#   D2: spec 1 has OVB from X1 (plim 0.62), specs 2-5 partial
#       X1 out and recover 0.50.                                 Delta = 0.12.
#   D3: D indep X.  Specs 1-4 (no M) target the total effect
#       0.3 + 0.4*0.6 = 0.54; spec 5 (with M) partials out
#       the mediator and identifies the direct effect 0.30.      Delta = 0.24.
#   D4: each spec omits one X_k with gamma_k = 0.05*(-1)^k.
#       Var(D_res) = 0.09 + 1 = 1.09.  Cov(D_res, Y_res) =
#       0.545 + 0.3*gamma_k.  Three specs at 0.530/1.09,
#       three at 0.560/1.09.                                     Delta = 0.030/1.09.
#   D5: tau values (0.3, 0.5, 0.7, 0.4, 0.6).                    Delta = 0.4.
#   D6: D indep X, all 10 specs target tau = 0.5.                Delta = 0.
#   D7: same DGP family as D1 (D indep X, clusters added).       Delta = 0.
#   D8: handled separately because Delta depends on pi
#       (see analytical_Delta_d8()).
#---------------------------------------------------------------------------

analytical_Delta <- list(
  D1 = 0,
  D2 = 0.12,
  D3 = 0.24,
  D4 = 0.030 / 1.09,
  D5 = 0.40,
  D6 = 0,
  D7 = 0
)

# Tolerance for the large-N verification.  0.005 catches anything beyond
# the second decimal place, well below the precision Delta is reported
# at in the paper.
delta_tol <- 0.005

verify_Delta <- function(dname, dgp_fn, n_pop = 500000L) {
  Delta_th <- analytical_Delta[[dname]]
  if (is.null(Delta_th)) {
    stop(sprintf("No analytical Delta defined for %s", dname))
  }
  pop_theta <- compute_pop_estimands(dgp_fn, n_pop = n_pop)
  Delta_emp <- max(pop_theta) - min(pop_theta)
  diff <- abs(Delta_emp - Delta_th)
  cat(sprintf("  Population estimands: %s\n",
              paste(sprintf("%.4f", pop_theta), collapse = ", ")))
  cat(sprintf("  Delta: analytical = %.5f, large-N = %.5f, diff = %.5f\n",
              Delta_th, Delta_emp, diff))
  if (diff > delta_tol) {
    stop(sprintf("Delta mismatch in %s: analytical=%.5f, large-N=%.5f (tol=%.4f)",
                 dname, Delta_th, Delta_emp, delta_tol))
  }
  Delta_th
}

#===========================================================================
# Single MC replication
#===========================================================================

one_mc_rep <- function(m, design_fn, n_obs, B, alpha) {

  # Generate data
  dgp <- design_fn(n_obs)

  # Full-sample estimates
  theta <- estimate_specs(dgp)

  # Bootstrap: resample individuals with replacement
  K <- dgp$K
  DRAWS <- matrix(NA_real_, nrow = B, ncol = K)

  for (b in 1:B) {
    idx <- sample.int(n_obs, n_obs, replace = TRUE)

    # Resample DGP
    dgp_b <- dgp
    if (isTRUE(dgp$multi_Y)) {
      dgp_b$Y_list <- lapply(dgp$Y_list, function(Y) Y[idx])
    } else {
      dgp_b$Y <- dgp$Y[idx]
    }
    dgp_b$specs <- lapply(dgp$specs, function(X) X[idx, , drop = FALSE])

    DRAWS[b, ] <- estimate_specs(dgp_b)
  }

  # Compute all statistics
  stats <- robustness_stats(theta, DRAWS, alpha)
  if (is.null(stats)) return(NULL)

  # Return as data.frame row
  data.frame(
    rep = m,
    K = stats$K,
    W_obs = stats$W_obs,
    R_obs = stats$R_obs,
    delta_W = stats$delta_W,
    delta_R_05 = stats$delta_R_05,
    delta_R_10 = stats$delta_R_10,
    delta_R_50 = stats$delta_R_50,
    p_W = stats$p_W,
    p_W_chi2 = stats$p_W_chi2,
    p_R = stats$p_R,
    eyeball_sign_sig = stats$eyeball_sign_sig,
    eyeball_ci_overlap = stats$eyeball_ci_overlap,
    range_rejects = stats$range_rejects
  )
}
