#' @export
print.range_test <- function(x, ...) {
  cat(strrep("-", 66), "\n", sep = "")
  cat(sprintf("  %s  (K=%d, B=%d)\n", x$label %||% "comparison", x$K, x$B))
  cat(strrep("-", 66), "\n", sep = "")
  if (x$B_dropped > 0L) {
    cat(sprintf("  (%d incomplete replication(s) dropped)\n", x$B_dropped))
  }
  cat(sprintf("  Equality      R = %9.4f   p_R = %7.4f\n", x$R, x$p_R))
  cat(sprintf("                W = %9.4f   p_W = %7.4f\n", x$W, x$p_W))
  cat("  Equivalence\n")
  for (i in seq_len(nrow(x$equivalence))) {
    cat(sprintf("    alpha=%.2f   delta*_R = %9.4f   delta*_W = %9.4f\n",
                x$equivalence$alpha[i], x$equivalence$delta_R[i],
                x$equivalence$delta_W[i]))
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
#' @param x A \code{"robustness"} object from \code{\link{robustness}}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.robustness <- function(x, ...) {
  cat(strrep("=", 66), "\n", sep = "")
  cat("  Robustness: range tests for equality and equivalence\n")
  cat("  across specifications (Jaeger 2026)\n")
  cat(sprintf("  B = %d   comparisons = %d   alpha = %s\n",
              x$B, length(x$results),
              paste(formatC(x$alpha, format = "g"), collapse = ", ")))
  cat(strrep("=", 66), "\n", sep = "")
  for (r in x$results) print(r)
  cat(strrep("=", 66), "\n", sep = "")
  invisible(x)
}

#' Coerce a robustness result to a data frame
#'
#' One row per comparison-by-alpha, with the range and Wald statistics,
#' equality p-values, and equivalence bounds. Convenient for assembling
#' tables.
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
      R          = r$R,
      W          = r$W,
      p_R        = r$p_R,
      p_W        = r$p_W,
      alpha      = r$equivalence$alpha,
      delta_R    = r$equivalence$delta_R,
      delta_W    = r$equivalence$delta_W,
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

# Null-coalescing helper.
`%||%` <- function(a, b) if (is.null(a)) b else a
