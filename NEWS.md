# robustness 0.1.2

* The Wald statistics (`W`, `p_W`, `Wstar`) now return `NA` with a warning when
  the contrast covariance is rank deficient (collinear specifications), rather
  than inverting it with a generalized inverse and reporting a degenerate Wald.
  Rank deficiency is detected with a relative eigenvalue tolerance, which
  catches structured collinearity (for example one specification a fixed shift
  or linear combination of others) that a default rank check can miss. The
  range statistics (`R`, `p_R`, `Rstar`) do not use the contrast covariance and
  are reported normally throughout.

* Each per-comparison result now carries a logical `wald_ok`, `FALSE` exactly
  when the Wald is undefined through rank deficiency. It is part of the public
  object and is exposed as a column by `as.data.frame()` and `summary()`. The
  printed output notes the rank-deficient case rather than showing a bare `NA`.

* `robustness()` now validates its inputs up front, matching the Stata
  command: `draws` must be numeric, `theta` finite, every `alpha` strictly
  between 0 and 1, and `max_drop` in `[0, 1)`. The printed header labels the
  supplied draw count as `B (supplied)`, since the count entering each
  comparison is shown per comparison after incomplete replications are dropped.
