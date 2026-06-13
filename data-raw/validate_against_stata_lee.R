# Validate the R `robustness` package against the Stata package, using the
# SAME Lee bootstrap draws. Because the draws are identical (not regenerated),
# the two implementations should agree to many decimal places, not merely up
# to Monte Carlo error. This is the definitive cross-language check.
#
# Usage: edit `draws_dir` to point at the Stata package's output/draws folder,
# then run from the package root:
#     Rscript data-raw/validate_against_stata_lee.R
#
# Requires the package to be installed or sourced, plus haven.

library(haven)
# If the package is installed: library(robustness)
# Otherwise source the core directly:
if (!requireNamespace("robustness", quietly = TRUE)) {
  source("R/robustness.R"); source("R/methods.R")
}

# ---- Point this at your Stata draws ----
draws_dir <- "~/Dropbox/Research/robustness/replication-package/output/draws"
stem      <- file.path(path.expand(draws_dir), "lee_draws")

# ---- Read the Stata three-file draws ----
meta  <- as.data.frame(read_dta(paste0(stem, "_meta.dta")))
draws <- as.data.frame(read_dta(paste0(stem, ".dta")))

theta <- meta$theta
K     <- length(theta)
D     <- as.matrix(draws[, paste0("coef", seq_len(K))])

# Per-spec n if present
ncols <- paste0("n", seq_len(K))
nmat  <- if (all(ncols %in% names(draws))) as.matrix(draws[, ncols]) else NULL

# ---- The Lee comparisons (paper order) ----
comps <- list(
  all            = 1:7,
  main           = 1:5,
  main_resid     = 1:6,
  main_firstdiff = c(1, 2, 3, 4, 5, 7)
)

res <- robustness(theta, D, comparisons = comps, alpha = c(0.50, 0.05), n = nmat)
df  <- as.data.frame(res)
cat("\n=== R package results on the Stata Lee draws ===\n")
print(df, row.names = FALSE)

# ---- Expected Stata values (from stats_lee.log) ----
# columns: comparison, delta_R(.50), delta_R(.05), p_R
expected <- data.frame(
  comparison = c("all", "main", "main_resid", "main_firstdiff"),
  dR50 = c(0.0128, 0.0032, 0.0088, 0.0076),
  dR05 = c(0.0270, 0.0081, 0.0242, 0.0214),
  p_R  = c(0.9208, 0.8706, 0.7690, 0.9191),
  stringsAsFactors = FALSE
)

cat("\n=== Comparison to Stata (stats_lee.log) ===\n")
tol <- 5e-4   # draws are identical; differences should be ~0, this is generous
all_ok <- TRUE
for (i in seq_len(nrow(expected))) {
  cn <- expected$comparison[i]
  r  <- res$results[[cn]]
  got_dR50 <- r$equivalence$delta_R[r$equivalence$alpha == 0.50]
  got_dR05 <- r$equivalence$delta_R[r$equivalence$alpha == 0.05]
  got_pR   <- r$p_R
  d1 <- abs(got_dR50 - expected$dR50[i])
  d2 <- abs(got_dR05 - expected$dR05[i])
  d3 <- abs(got_pR   - expected$p_R[i])
  ok <- (d1 < tol) && (d2 < tol) && (d3 < tol)
  all_ok <- all_ok && ok
  cat(sprintf("  %-15s  dR50 %.4f vs %.4f (d=%.5f)  dR05 %.4f vs %.4f (d=%.5f)  p_R %.4f vs %.4f (d=%.5f)  [%s]\n",
              cn, got_dR50, expected$dR50[i], d1,
              got_dR05, expected$dR05[i], d2,
              got_pR, expected$p_R[i], d3, if (ok) "PASS" else "FAIL"))
}
cat(if (all_ok) "\nALL MATCH: R reproduces Stata on the same draws.\n"
    else "\nMISMATCH: investigate.\n")
cat("\nNote: the Stata values above are rounded to 4 dp (read from the log),\n",
    "so differences up to ~5e-5 are just rounding. Larger gaps in p_R or\n",
    "delta* on the SAME draws would indicate a real discrepancy to chase.\n", sep = "")
