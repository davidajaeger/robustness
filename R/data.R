#' Example bootstrap draws in the shape of the Lee (2008) application
#'
#' A small, synthetic illustration shaped like the Lee (2008) incumbency
#' regression-discontinuity application used in Jaeger (2026): seven
#' specifications of an effect near 0.077, with correlated bootstrap draws.
#' Provided for documentation and examples. This is NOT the actual Lee
#' bootstrap draws, which are large and external and are produced by the Stata
#' replication package. The numbers here are illustrative, not a replication.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{theta}{Named numeric vector of 7 full-sample point estimates.}
#'   \item{draws}{A 9999-by-7 numeric matrix of uncentred bootstrap draws,
#'     columns matching \code{theta}.}
#'   \item{comparisons}{A named list of four comparisons (integer column
#'     vectors): \code{all}, \code{main}, \code{main_resid},
#'     \code{main_firstdiff}.}
#' }
#'
#' @references
#' Jaeger, David A. (2026). Robustness? Range Tests for Equality and
#' Equivalence Across Specifications.
#'
#' @examples
#' robustness(lee_example$theta, lee_example$draws,
#'            comparisons = lee_example$comparisons)
"lee_example"
