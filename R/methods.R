#' @export
print.range_test <- function(x, ...) {
  cat(strrep("-", 70), "\n", sep = "")
  cat(sprintf("  %s  (K=%d, B=%d)\n", x$label %||% "comparison", x$K, x$B))
  cat(strrep("-", 70), "\n", sep = "")
  if (x$B_dropped > 0L) {
    cat(sprintf("  (%d incomplete replication(s) dropped)\n", x$B_dropped))
  }
  cat(sprintf("  theta_bar     = %9.5f\n", x$theta_bar))
  cat(sprintf("  Equality      R = %9.5f   p_R = %7.4f\n", x$R, x$p_R))
  cat(sprintf("                W = %9.5f   p_W = %7.4f\n", x$W, x$p_W))
  cat("  Equivalence\n")
  for (i in seq_len(nrow(x$equivalence))) {
    cat(sprintf("    alpha=%.2f   R* = %9.5f   W* = %9.5f\n",
                x$equivalence$alpha[i], x$equivalence$Rstar[i],
                x$equivalence$Wstar[i]))
  }
  if (is.na(x$ratio)) {
    cat("  Robustness ratio = .         (theta_bar = 0, undefined)\n")
  } else {
    cat(sprintf("  Robustness ratio = %9.5f   (R*(.95) / |theta_bar|)\n",
                x$ratio))
  }
  if (!is.null(x$p_tau) && !is.na(x$p_tau)) {
    # Equivalence at the pre-specified tolerance tau. By the duality of
    # Jaeger (2026), H0: Delta >= tau rejects at level alpha exactly when
    # R*(1-alpha) <= tau.
    cat(sprintf("  Equivalence   tau = %9.5f   p_tau = %7.4f\n",
                x$tau, x$p_tau))
  }
  if (isFALSE(x$wald_ok)) {
    cat("  Wald undefined: contrast covariance rank deficient.\n")
    cat("  Range statistics above are valid.\n")
  }
  if (!is.null(x$avg_n)) {
    nmin <- min(x$avg_n); nmax <- max(x$avg_n)
    if (isTRUE(all.equal(nmin, nmax))) {
      cat(sprintf("    avg n = %g (equal across specs)\n", round(nmin)))
    } else {
      cat(sprintf("    avg n per spec: %g to %g (samples differ across specs)\n",
                  round(nmin), round(nmax)))
    }
  }
  invisible(x)
}

#' Print a robustness result
#'
#' Two-panel layout matching the Stata counterpart's reporting. Panel A
#' (specification-level estimates) is printed only when the result carries one
#' (\code{se} was supplied to \code{robustness}). Panel B (per-comparison
#' statistics) is always printed. After the panels, each comparison's full
#' detail (Wald, equivalence bounds, ratio) is shown.
#'
#' @param x A \code{"robustness"} object from \code{\link{robustness}}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.robustness <- function(x, ...) {
  cat(strrep("=", 70), "\n", sep = "")
  cat("  Robustness: range tests for equality and equivalence\n")
  cat("  across specifications (Jaeger 2026)\n")
  cat(sprintf("  B = %d (supplied)   comparisons = %d   alpha = %s\n",
              x$B, length(x$results),
              paste(formatC(x$alpha, format = "g"), collapse = ", ")))
  cat(strrep("=", 70), "\n", sep = "")

  # Panel A: specification-level estimates, only when available.
  if (!is.null(x$panel_a)) {
    cat("\n")
    cat(strrep("-", 70), "\n", sep = "")
    cat("  Panel A: Specification-level estimates\n")
    cat(strrep("-", 70), "\n", sep = "")
    pa <- x$panel_a
    show_n_full <- !all(is.na(pa$n_full))
    show_n_boot <- !all(is.na(pa$n_boot))
    .print_panel_a_rows(pa, show_n_full, show_n_boot)
  }

  # Panel B: comparison-level statistics, always.
  cat("\n")
  cat(strrep("-", 70), "\n", sep = "")
  cat("  Panel B: Comparison-set statistics\n")
  cat(strrep("-", 70), "\n", sep = "")
  .print_panel_b_rows(x$results)
  cat("  Note: the robustness ratio is R*(.95) / |theta_bar|. When |theta_bar|\n")
  cat("        is close to zero, interpret it with caution.\n")

  # Equivalence at a pre-specified tolerance tau, printed only when tau was
  # supplied. Mirrors the Stata command's block below Panel B: R*(.95), the
  # equivalence p-value p_tau, and whether H0: Delta >= tau rejects at the
  # .05 level (equivalently, by the duality lemma, R*(.95) <= tau).
  tau_val <- x$results[[1]]$tau
  if (!is.null(tau_val) && !is.na(tau_val)) {
    cat("\n")
    cat(strrep("-", 70), "\n", sep = "")
    cat(sprintf("  Equivalence at pre-specified tolerance  tau = %.5f\n", tau_val))
    cat(strrep("-", 70), "\n", sep = "")
    .print_tau_rows(x$results, tau_val)
    cat("  A comparison is equivalent at tau (rejects H0: Delta >= tau at the\n")
    cat("  .05 level) exactly when R*(.95) <= tau, equivalently p_tau <= .05.\n")
  }

  # Per-comparison detail (Wald, equivalence bounds, full breakdown).
  cat("\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat("  Per-comparison detail\n")
  cat(strrep("=", 70), "\n", sep = "")
  for (r in x$results) print(r)
  cat(strrep("=", 70), "\n", sep = "")
  invisible(x)
}

# Internal: print Panel A as right-aligned columns.
.print_panel_a_rows <- function(pa, show_n_full, show_n_boot) {
  spec_w <- max(nchar(pa$spec), 4L) + 2L
  hdr <- paste0(
    "  ", formatC("Spec", width = -spec_w),
    formatC("theta_hat", width = 12),
    " ", formatC("SE", width = 8)
  )
  if (show_n_full) hdr <- paste0(hdr, " ", formatC("n_full",     width = 10))
  if (show_n_boot) hdr <- paste0(hdr, " ", formatC("avg_n_boot", width = 10))
  cat(hdr, "\n", sep = "")
  for (i in seq_len(nrow(pa))) {
    row <- paste0(
      "  ", formatC(pa$spec[i], width = -spec_w),
      formatC(pa$theta[i], width = 12, digits = 5, format = "f"),
      " ", formatC(pa$se[i],   width = 8,  digits = 5, format = "f")
    )
    if (show_n_full) {
      v <- if (is.na(pa$n_full[i])) formatC(".", width = 10) else
        formatC(round(pa$n_full[i]), width = 10, format = "f", digits = 0)
      row <- paste0(row, " ", v)
    }
    if (show_n_boot) {
      v <- if (is.na(pa$n_boot[i])) formatC(".", width = 10) else
        formatC(round(pa$n_boot[i]), width = 10, format = "f", digits = 0)
      row <- paste0(row, " ", v)
    }
    cat(row, "\n", sep = "")
  }
}

# Internal: print Panel B as right-aligned columns.
.print_panel_b_rows <- function(results) {
  cnames <- vapply(results, function(r) r$label %||% "comparison", character(1))
  cn_w <- max(nchar(cnames), nchar("Comparison set")) + 2L
  hdr <- paste0(
    "  ", formatC("Comparison set", width = -cn_w),
    formatC("K",          width = 3),
    " ", formatC("theta_bar",  width = 10),
    " ", formatC("R(theta)",   width = 10),
    " ", formatC("R*(.50)",    width = 10),
    " ", formatC("R*(.95)",    width = 10),
    " ", formatC("p_R",        width = 8),
    " ", formatC("Rob. ratio", width = 10)
  )
  cat(hdr, "\n", sep = "")
  for (r in results) {
    # R*(.95): always reported as the bootstrap .95 quantile of the uncentred
    # range, computed by .range_test_one independent of the user's alpha
    # argument (so it has a fixed meaning and is comparable across calls).
    # R*(.50): pulled from the equivalence frame if alpha = 0.50 is in the
    # user's request; NA otherwise.
    rstar_95 <- r$Rstar_95
    eq50_idx <- which(r$equivalence$alpha == 0.50)
    rstar_50 <- if (length(eq50_idx)) r$equivalence$Rstar[eq50_idx[1]] else NA_real_
    # Ratio: right-justified period in width 10 when undefined, matching
    # the numeric format used otherwise.
    ratio_s  <- if (is.na(r$ratio)) "         ." else
      formatC(r$ratio, width = 10, digits = 4, format = "f")
    cat(
      "  ", formatC(r$label %||% "comparison", width = -cn_w),
      formatC(r$K,           width = 3),
      " ", formatC(r$theta_bar, width = 10, digits = 5, format = "f"),
      " ", formatC(r$R,         width = 10, digits = 5, format = "f"),
      " ", formatC(rstar_50,    width = 10, digits = 5, format = "f"),
      " ", formatC(rstar_95,    width = 10, digits = 5, format = "f"),
      " ", formatC(r$p_R,       width = 8,  digits = 4, format = "f"),
      " ", ratio_s,
      "\n", sep = ""
    )
  }
}

# Internal: print the equivalence-at-tolerance rows. One row per comparison
# with R*(.95), p_tau, and the reject/accept decision at the .05 level.
.print_tau_rows <- function(results, tau_val) {
  cnames <- vapply(results, function(r) r$label %||% "comparison", character(1))
  cn_w <- max(nchar(cnames), nchar("Comparison set")) + 2L
  hdr <- paste0(
    "  ", formatC("Comparison set", width = -cn_w),
    formatC("R*(.95)", width = 10),
    " ", formatC("p_tau", width = 8),
    " ", formatC("equiv. at tau", width = 14)
  )
  cat(hdr, "\n", sep = "")
  for (r in results) {
    equiv <- if (isTRUE(r$Rstar_95 <= tau_val)) "yes" else "no"
    cat(
      "  ", formatC(r$label %||% "comparison", width = -cn_w),
      formatC(r$Rstar_95, width = 10, digits = 5, format = "f"),
      " ", formatC(r$p_tau, width = 8, digits = 4, format = "f"),
      " ", formatC(equiv, width = 14),
      "\n", sep = ""
    )
  }
}

#' Coerce a robustness result to a data frame
#'
#' One row per comparison-by-alpha, with the range and Wald statistics,
#' equality p-values, equivalence bounds, the comparison mean and robustness
#' ratio, the \code{wald_ok} flag (\code{FALSE} when the Wald is undefined
#' through rank deficiency), and (when \code{tau} was supplied to
#' \code{robustness}) the tolerance \code{tau} and its equivalence p-value
#' \code{p_tau}. Convenient for assembling tables.
#'
#' @param x A \code{"robustness"} object.
#' @param row.names,optional Ignored, present for S3 consistency.
#' @param ... Unused.
#' @return A data frame.
#' @export
as.data.frame.robustness <- function(x, row.names = NULL, optional = FALSE,
                                     ...) {
  rows <- lapply(x$results, function(r) {
    data.frame(
      comparison = r$label,
      K          = r$K,
      B          = r$B,
      theta_bar  = r$theta_bar,
      R          = r$R,
      W          = r$W,
      p_R        = r$p_R,
      p_W        = r$p_W,
      wald_ok    = r$wald_ok,
      alpha      = r$equivalence$alpha,
      Rstar      = r$equivalence$Rstar,
      Wstar      = r$equivalence$Wstar,
      ratio      = r$ratio,
      tau        = if (is.null(r$tau)) NA_real_ else r$tau,
      p_tau      = if (is.null(r$p_tau)) NA_real_ else r$p_tau,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}

#' Summarise a robustness result
#'
#' @param object A \code{"robustness"} object.
#' @param ... Unused.
#' @return A data frame, as from \code{as.data.frame}, returned invisibly
#'   after printing.
#' @export
summary.robustness <- function(object, ...) {
  df <- as.data.frame(object)
  print(df, row.names = FALSE)
  invisible(df)
}

#' Extract the per-replication bootstrap statistics
#'
#' Returns the bootstrap series that the reported summaries collapse to
#' scalars: for each comparison, the uncentred range and Wald (whose
#' \code{1 - alpha} quantiles are \code{R*} and \code{W*}) and the recentred
#' range and Wald (whose tails at or above the observed statistic give
#' \code{p_R} and \code{p_W}). Long form, one row per comparison-by-draw,
#' matching the layout the Stata command's \code{saving()} option writes.
#' Requires the object to have been produced with \code{keep_draws = TRUE}.
#'
#' @param x A \code{"robustness"} object created with \code{keep_draws = TRUE}.
#' @return A data frame with columns \code{comparison}, \code{draw},
#'   \code{range_unc}, \code{range_rc}, \code{wald_unc}, \code{wald_rc}.
#' @export
bootstrap_draws <- function(x) {
  if (!inherits(x, "robustness"))
    stop("x must be a 'robustness' object.")
  have <- vapply(x$results, function(r) !is.null(r$draws), logical(1))
  if (!all(have))
    stop("No stored bootstrap draws. Re-run robustness(..., keep_draws = TRUE).")
  rows <- lapply(x$results, function(r)
    data.frame(comparison = r$label, r$draws, stringsAsFactors = FALSE))
  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}

#' Extract Panel A from a robustness result
#'
#' Returns the specification-level data frame (one row per specification) when
#' the result was produced with full-sample SEs supplied to
#' \code{\link{robustness}}, and \code{NULL} otherwise.
#'
#' @param x A \code{"robustness"} object.
#' @return A data frame with columns \code{spec}, \code{theta}, \code{se},
#'   \code{n_full}, \code{n_boot}, or \code{NULL}.
#' @export
panel_a <- function(x) {
  if (!inherits(x, "robustness"))
    stop("x must be a 'robustness' object.")
  x$panel_a
}

# Null-coalescing helper.
`%||%` <- function(a, b) if (is.null(a)) b else a
