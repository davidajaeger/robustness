# robustness 0.2.0

* New: `read_robustness_dta()` reads the three canonical `.dta` files
  produced by the Stata package's bootstrap-generation step (meta, draws,
  comps) and returns a list ready to pass to `robustness()`. Convenience for
  users who run the Stata-side pipeline and want to consume its outputs from
  R; `robustness()` itself still takes only in-memory inputs. Uses `haven`
  (in `Suggests`).

* New: each per-comparison result carries `theta_bar` (the mean of the
  comparison's full-sample estimates) and `ratio` (the robustness ratio
  `R*(.95) / |theta_bar|`, `NA` when `theta_bar = 0`). Both are reported in
  the printed output and as columns of `as.data.frame()`. The ratio always
  uses the .95 quantile of the bootstrap range, regardless of what is
  supplied to `alpha`, so it has a fixed meaning and is comparable across
  calls. `R*(.95)` is also exposed as a top-level field on each result.

* New: `robustness()` accepts optional `se`, `labels`, `n_full`, and `n_boot`
  arguments. When `se` is supplied, the result carries a Panel A data frame
  with one row per specification (label, theta, SE, full-sample n,
  bootstrap-average n), matching Panel A in the Stata counterpart. Without
  `se`, Panel A is omitted and only Panel B is produced. The new accessor
  `panel_a()` returns the Panel A data frame, or `NULL` when absent.

* The print method now emits a two-panel layout (Panel A when available,
  Panel B always) followed by per-comparison detail blocks. Panel B has one
  row per comparison with `K`, `theta_bar`, `R(theta)`, `R*(.50)`,
  `R*(.95)`, `p_R`, and the robustness ratio, matching the Stata layout.

* Stricter input validation: `draws` may contain `NA` (incomplete bootstrap
  replications, dropped downstream by `max_drop`) but must not contain `Inf`
  or `-Inf`. Comparison indices must be whole numbers; non-integer values
  (for example `c(1.5, 2.7)`) are rejected rather than silently truncated.

* Breaking: the optional `n` argument is renamed to `n_boot`. The
  `n_boot` matrix is also supported as a length-`K` vector of per-spec
  averages (the old `n` accepted either). Reproducing earlier results is
  unaffected because the simulations and applications shipped with the
  package do not pass anything through this argument.

* Reproducibility: the underlying mathematics is unchanged from 0.1.x. The
  recentring, the type-1 quantile for `R*`, the `(1 + r)/(B + 1)` p-value,
  and the eigenvalue-based `wald_ok` test all match earlier behaviour
  exactly. Simulation results produced with 0.1.x are reproduced bit-for-bit
  by 0.2.0 from the same seed and draws.


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
