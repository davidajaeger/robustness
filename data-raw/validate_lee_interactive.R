# Interactive validation of the R `robustness` package against the Stata
# package, using the SAME Lee draws. Paste into the RStudio console after
# loading the package. Because the draws are identical (not regenerated), the
# two should agree to rounding, not merely up to Monte Carlo error.
#
# Setup (run once):
#   devtools::load_all()            # from the package project, OR
#   library(robustness)             # if installed
#   install.packages("haven")       # if not already installed

library(haven)

# ---- EDIT this to your Stata draws folder ----
stem <- "~/Dropbox/Research/robustness/replication-package/output/draws/lee_draws"
stem <- path.expand(stem)

# ---- Read the Stata three-file draws ----
meta  <- as.data.frame(read_dta(paste0(stem, "_meta.dta")))
draws <- as.data.frame(read_dta(paste0(stem, ".dta")))

theta <- meta$theta
K     <- length(theta)
D     <- as.matrix(draws[, paste0("coef", seq_len(K))])

# Per-spec sample sizes, if the draws carry them
ncols <- paste0("n", seq_len(K))
nmat  <- if (all(ncols %in% names(draws))) as.matrix(draws[, ncols]) else NULL

# ---- Run the R package ----
res <- robustness(
  theta, D,
  comparisons = list(all            = 1:7,
                     main           = 1:5,
                     main_resid     = 1:6,
                     main_firstdiff = c(1, 2, 3, 4, 5, 7)),
  alpha = c(0.50, 0.05),
  n     = nmat
)

# Look at the results
res                    # printed, Stata-log style
as.data.frame(res)     # tidy table to compare against stats_lee.log

# ---- Side-by-side with the Stata numbers (from stats_lee.log) ----
stata <- data.frame(
  comparison = c("all", "main", "main_resid", "main_firstdiff"),
  s_dR50 = c(0.0128, 0.0032, 0.0088, 0.0076),
  s_dR05 = c(0.0270, 0.0081, 0.0242, 0.0214),
  s_pR   = c(0.9208, 0.8706, 0.7690, 0.9191)
)
r_df <- as.data.frame(res)
compare <- merge(
  stata,
  data.frame(
    comparison = r_df$comparison[r_df$alpha == 0.05],
    r_dR05     = r_df$delta_R[r_df$alpha == 0.05],
    r_pR       = r_df$p_R[r_df$alpha == 0.05]
  ),
  by = "comparison"
)
dR50 <- r_df$delta_R[r_df$alpha == 0.50]
names(dR50) <- r_df$comparison[r_df$alpha == 0.50]
compare$r_dR50 <- dR50[compare$comparison]
compare$d_dR50 <- abs(compare$r_dR50 - compare$s_dR50)
compare$d_dR05 <- abs(compare$r_dR05 - compare$s_dR05)
compare$d_pR   <- abs(compare$r_pR   - compare$s_pR)
print(compare[, c("comparison", "s_dR50", "r_dR50", "d_dR50",
                  "s_dR05", "r_dR05", "d_dR05",
                  "s_pR", "r_pR", "d_pR")], row.names = FALSE)

cat("\nDifferences (d_*) should be ~1e-4 or smaller (Stata values are rounded\n",
    "to 4 dp in the log). Same draws means same statistics, not just close.\n", sep = "")
