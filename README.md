# robustness (R)

Range tests for equality and equivalence across specifications, from saved bootstrap draws. This is the R implementation of the method in Jaeger (2026), "Robustness? Range Tests for Equality and Equivalence Across Specifications." It is the counterpart of the Stata [`robustness`](https://github.com/davidajaeger/robustness-stata) command: given the same bootstrap draws, the two produce the same statistics.

## What it does

Applied work routinely presents several specifications of the same coefficient and declares the results "robust" when the estimates look similar. The range of the estimates is the implicit object of that claim, but its joint sampling distribution is rarely reported. This package computes, for each comparison:

- **`R*(.95)` (the minimum equivalence bound):** the smallest tolerance within which the specifications can be certified equivalent at the 5 percent level, in the units of the coefficient. It is the `(1 - alpha)` quantile of the uncentred bootstrap range.
- **`R*(.50)`:** the median of the bootstrap range, a point estimate of the range.
- **`p_R` (the range-based equality test):** the bootstrap p-value for the null that every specification shares a common probability limit.
- the robustness ratio `R*(.95) / |theta_bar|`, where `theta_bar` is the mean estimate over the comparison's specs.

The Wald statistic (`W`), its p-value (`p_W`), and the Wald bound (`W*`) are also computed, with a `wald_ok` flag for collinearity-induced rank deficiency (the Wald is returned as `NA` with a warning; the range statistics are unaffected). With `keep_draws = TRUE` the per-replication bootstrap series are retained for plotting.

The package **computes statistics from draws; it does not generate the draws.** You supply the full-sample estimates and a matrix of uncentred (raw) bootstrap draws in which the same resampled units were used for every specification on each replication. Resampling specifications independently destroys the joint distribution and silently produces wrong results.

## Installation

```r
# install.packages("remotes")
remotes::install_github("davidajaeger/robustness")
```

The package depends only on base R (>= 3.5.0) and `stats`. The optional `read_robustness_dta()` reader for Stata .dta files uses `haven`; install it if you plan to consume Stata pipeline outputs (`install.packages("haven")`).

## Quick start

A synthetic example shaped like the Lee (2008) application ships with the package:

```r
library(robustness)
robustness(lee_example$theta, lee_example$draws,
           comparisons = lee_example$comparisons)
```

A minimal example from scratch:

```r
set.seed(1)
theta <- c(0.20, 0.22, 0.19, 0.21)
draws <- sapply(theta, function(m) rnorm(9999, m, 0.03))   # 9999 x 4 uncentred draws

robustness(theta, draws)

# Several comparisons at once
robustness(theta, draws,
           comparisons = list(all = 1:4, first_two = 1:2, extremes = c(1, 4)))
```

## Two workflows

The package supports two ways of supplying inputs.

### Workflow 1: R-native draws

Build `theta` (length-K numeric vector) and `draws` (B × K numeric matrix) from your own R bootstrap. The column order of `draws` must match `theta`. Optionally supply `se` (length-K vector of full-sample standard errors), `labels`, `n_full`, and `n_boot` for Panel A reporting.

```r
result <- robustness(theta, draws,
                     comparisons = list(all = 1:K),
                     se     = full_sample_se,
                     labels = c("Baseline", "+X1", "+X2", "+X1+X2"),
                     n_full = c(2000, 2000, 2000, 2000))
print(result)
```

If `se` is supplied, the print method emits a two-panel layout (Panel A: per-spec; Panel B: per-comparison). If not, only Panel B and the per-comparison detail blocks are shown.

### Workflow 2: Stata pipeline files

If you run the Stata `robustness` generation pipeline and want to consume its three canonical `.dta` files from R, use `read_robustness_dta()`:

```r
x <- read_robustness_dta(
  meta  = "output/draws/lee_draws_meta.dta",
  draws = "output/draws/lee_draws.dta",
  comps = "output/draws/lee_draws_comps.dta"
)
result <- robustness(x$theta, x$draws, x$comparisons,
                     se = x$se, labels = x$labels,
                     n_full = x$n_full, n_boot = x$n_boot)
print(result)
```

`read_robustness_dta()` returns a plain list whose components map directly onto the arguments of `robustness()`. The reader uses `haven`.

## Output

`robustness()` returns an object of class `"robustness"` with `print`, `summary`, and `as.data.frame` methods. The `as.data.frame` method returns one row per comparison-by-alpha, with columns `comparison`, `K`, `B`, `theta_bar`, `R`, `W`, `p_R`, `p_W`, `wald_ok`, `alpha`, `Rstar`, `Wstar`, `ratio`. `panel_a()` returns the Panel A data frame when available, `NULL` otherwise.

## Plotting the bootstrap distribution

The range distribution is the object the "robustness" claim implicitly invokes and rarely shows. With `keep_draws = TRUE` you can plot it:

```r
r <- robustness(theta, draws, keep_draws = TRUE)

# Per-replication series: comparison, draw, range_unc, range_rc, wald_unc, wald_rc
d <- bootstrap_draws(r)

# Range distribution for one comparison, with the equivalence bound marked.
# Read R* from the result rather than recomputing it: R* is a type-1 order
# statistic, so quantile() with its default (type 7) would interpolate and
# the line would sit a hair off. (If you do recompute, pass type = 1.)
tab   <- as.data.frame(r)
bound <- tab$Rstar[tab$comparison == "all" & tab$alpha == 0.05]   # R*(.95)
hist(d$range_unc[d$comparison == "all"], breaks = 40,
     main = "Bootstrap range distribution")
abline(v = bound, lty = 2)
```

The uncentred series (`range_unc`, `wald_unc`) are the ones whose `(1 - alpha)` quantiles are `R*` and `W*`; the recentred series (`range_rc`, `wald_rc`) are the ones whose tails at or above the observed statistic give `p_R` and `p_W`. This matches the layout the Stata command's `saving()` option writes, so either package gives the same series.

## Generating the draws

The package consumes draws; it does not produce them. The vignette (`vignette("robustness")`) gives a worked generation example, and the [Stata replication package](https://github.com/davidajaeger/robustness-replications) for Jaeger (2026) contains the generation scripts used for each application in the paper. The single requirement is that every specification on a given replication is estimated on the **same** resampled units.

## Citing

The software implements the method; please cite the paper.

```r
citation("robustness")
```

> Jaeger, David A. (2026). Robustness? Range Tests for Equality and Equivalence Across Specifications.

## License

MIT.
