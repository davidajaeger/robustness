# robustness (R)

Range tests for equality and equivalence across specifications, from saved bootstrap draws. This is the R implementation of the method in Jaeger (2026), "Robustness? Range Tests for Equality and Equivalence Across Specifications." It is the counterpart of the Stata [`robustness`](https://github.com/davidajaeger/robustness-stata) command: given the same bootstrap draws, the two produce the same statistics.

## What it does

Applied work routinely presents several specifications of the same coefficient and declares the results "robust" when the estimates look similar. The range of the estimates is the implicit object of that claim, but its joint sampling distribution is rarely reported. This package computes two statistics from bootstrap draws of the coefficient:

- **`R*` (the minimum equivalence bound):** the smallest tolerance within which the specifications can be certified equivalent at a chosen level. Reported in the units of the coefficient.
- **`p_R` (the range-based equality test):** the bootstrap p-value for the null that every specification shares a common probability limit.

The two answer different questions. A large `p_R` means the estimates cannot be told apart; a small `R*` means they are affirmatively close. Imprecise estimates do not earn robustness by default.

The package **computes statistics from draws; it does not generate the draws.** You supply the full-sample estimates and a matrix of uncentred (raw) bootstrap draws in which the same resampled units were used for every specification on each replication. Resampling specifications independently destroys the joint distribution and silently produces wrong results.

## Installation

```r
# install.packages("remotes")
remotes::install_github("davidajaeger/robustness")
```

The package depends only on base R (>= 3.5.0) and `stats`.

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
theta <- c(0.20, 0.22, 0.19, 0.21)            # 4 full-sample estimates
draws <- sapply(theta, function(m) rnorm(999, m, 0.03))  # 999 x 4 uncentred draws

robustness(theta, draws)

# several comparison sets at once
robustness(theta, draws,
           comparisons = list(all = 1:4, first_two = 1:2, extremes = c(1, 4)))
```

## Interface

```r
robustness(theta, draws, comparisons = NULL,
           alpha = c(0.50, 0.05), n = NULL, max_drop = 0.01)
```

- **`theta`** — numeric vector of `K` full-sample point estimates.
- **`draws`** — `B` x `K` numeric matrix of uncentred bootstrap draws; column `k` matches `theta[k]`.
- **`comparisons`** — a single integer vector of column indices, or a named list of integer vectors for several comparisons. Defaults to all specifications as one comparison.
- **`alpha`** — significance levels for the equivalence bounds. Defaults to `c(0.50, 0.05)`: the median bound (a point estimate of the range) and the 95th-percentile upper bound.
- **`n`** — optional per-specification sample sizes, reported descriptively.
- **`max_drop`** — maximum proportion of incomplete replications tolerated before the function stops. Defaults to `0.01`.

Returns an object of class `robustness` with `print`, `summary`, and `as.data.frame` methods, so results drop straight into a data frame for tables.

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
