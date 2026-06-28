#' Read robustness inputs from Stata .dta files
#'
#' Reads the three files produced by the Stata package's bootstrap-generation
#' step (metadata, draws, comparisons) and returns the in-memory objects that
#' \code{\link{robustness}} expects. A convenience for users who run the
#' Stata-side generation pipeline and want to consume the canonical .dta
#' outputs from R; \code{robustness} itself takes only in-memory inputs, so
#' users with their own R-based bootstrap can skip this function entirely.
#'
#' @details
#' The reader expects the three-file convention documented in the Stata
#' \code{robustness} command:
#'
#' \itemize{
#'   \item \code{meta}: one row per specification, variables \code{k label
#'     theta se}, with optional \code{n}. \code{k} indexes specifications
#'     \code{1..K} and the rows are sorted by it.
#'   \item \code{draws}: one row per bootstrap replication, variables
#'     \code{coef1 coef2 ...} (required), optionally \code{n1 n2 ...}. Any
#'     other variables in the file (\code{rep}, \code{se1 se2 ...}) are
#'     ignored.
#'   \item \code{comps}: one row per comparison, variables \code{comp_name}
#'     (string) and \code{comp_cols} (space-separated 1-indexed column list).
#' }
#'
#' All three files must be passed; there is no default. The function reads
#' them with \code{haven::read_dta} and returns a plain list with no class
#' attached.
#'
#' @param meta  Path to the metadata .dta.
#' @param draws Path to the draws .dta.
#' @param comps Path to the comparisons .dta.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{theta}}{Length-\eqn{K} numeric vector of full-sample estimates.}
#'   \item{\code{se}}{Length-\eqn{K} numeric vector of full-sample SEs.}
#'   \item{\code{labels}}{Length-\eqn{K} character vector of specification labels.}
#'   \item{\code{draws}}{\eqn{B} by \eqn{K} numeric matrix of bootstrap draws.}
#'   \item{\code{comparisons}}{Named list of integer vectors of column indices.}
#'   \item{\code{n_full}}{Length-\eqn{K} numeric vector of full-sample sizes,
#'     or \code{NULL} if absent from \code{meta}.}
#'   \item{\code{n_boot}}{\eqn{B} by \eqn{K} numeric matrix of per-rep sample
#'     sizes, or \code{NULL} if absent from \code{draws}.}
#' }
#' Pass the list directly into \code{\link{robustness}}; see examples.
#'
#' @examples
#' \dontrun{
#' x <- read_robustness_dta(
#'   meta  = "output/draws/lee_draws_meta.dta",
#'   draws = "output/draws/lee_draws.dta",
#'   comps = "output/draws/lee_draws_comps.dta"
#' )
#' r <- robustness(x$theta, x$draws, x$comparisons,
#'                 se = x$se, labels = x$labels,
#'                 n_full = x$n_full, n_boot = x$n_boot)
#' print(r)
#' }
#'
#' @seealso \code{\link{robustness}}
#' @export
read_robustness_dta <- function(meta, draws, comps) {
  if (!requireNamespace("haven", quietly = TRUE))
    stop("The 'haven' package is required to read .dta files; install it with ",
         "install.packages(\"haven\").")

  # --- meta ---
  if (!file.exists(meta))
    stop("meta file not found: ", meta)
  meta_df <- as.data.frame(haven::read_dta(meta))
  req_meta <- c("k", "label", "theta", "se")
  miss_meta <- setdiff(req_meta, names(meta_df))
  if (length(miss_meta))
    stop("meta file missing required variables: ",
         paste(miss_meta, collapse = ", "))

  K <- nrow(meta_df)
  if (K < 2L)
    stop("meta file has ", K, " row(s); need at least 2 specifications.")

  # k must be a permutation of 1..K. Sort by k so position equals spec number.
  k_vals <- as.integer(meta_df$k)
  if (anyNA(k_vals))
    stop("meta variable 'k' must be a nonmissing integer in every row.")
  if (!setequal(k_vals, seq_len(K)))
    stop("meta variable 'k' must list each specification exactly once over ",
         "1..", K, " (found values: ",
         paste(sort(unique(k_vals)), collapse = ", "), ").")
  meta_df <- meta_df[order(k_vals), , drop = FALSE]

  theta  <- as.numeric(meta_df$theta)
  se     <- as.numeric(meta_df$se)
  labels <- as.character(meta_df$label)
  n_full <- if ("n" %in% names(meta_df)) as.numeric(meta_df$n) else NULL

  # --- draws ---
  if (!file.exists(draws))
    stop("draws file not found: ", draws)
  draws_df <- as.data.frame(haven::read_dta(draws))

  coef_vars <- paste0("coef", seq_len(K))
  miss_coef <- setdiff(coef_vars, names(draws_df))
  if (length(miss_coef))
    stop("draws file missing required variables: ",
         paste(miss_coef, collapse = ", "),
         " (expected coef1..coef", K, " matching ", K,
         " specifications in meta).")

  draws_mat <- as.matrix(draws_df[, coef_vars, drop = FALSE])
  storage.mode(draws_mat) <- "double"

  # Optional per-rep n columns: only attached if present for every spec.
  n_vars   <- paste0("n", seq_len(K))
  if (all(n_vars %in% names(draws_df))) {
    n_boot <- as.matrix(draws_df[, n_vars, drop = FALSE])
    storage.mode(n_boot) <- "double"
  } else {
    n_boot <- NULL
  }

  # --- comps ---
  if (!file.exists(comps))
    stop("comps file not found: ", comps)
  comps_df <- as.data.frame(haven::read_dta(comps))
  req_comps <- c("comp_name", "comp_cols")
  miss_comps <- setdiff(req_comps, names(comps_df))
  if (length(miss_comps))
    stop("comps file missing required variables: ",
         paste(miss_comps, collapse = ", "))
  if (nrow(comps_df) < 1L)
    stop("comps file has no rows.")

  cnames  <- as.character(comps_df$comp_name)
  ccols   <- as.character(comps_df$comp_cols)
  comparisons <- lapply(ccols, function(s) {
    tok <- strsplit(trimws(s), "\\s+")[[1]]
    if (!length(tok))
      stop("comps file has a row with empty comp_cols.")
    v <- suppressWarnings(as.integer(tok))
    if (anyNA(v))
      stop("comps file has non-integer entries in comp_cols: '", s, "'.")
    if (any(v < 1L | v > K))
      stop("comps file references a column outside 1..", K,
           " in row: '", s, "'.")
    v
  })
  names(comparisons) <- cnames

  list(
    theta       = theta,
    se          = se,
    labels      = labels,
    draws       = draws_mat,
    comparisons = comparisons,
    n_full      = n_full,
    n_boot      = n_boot
  )
}
