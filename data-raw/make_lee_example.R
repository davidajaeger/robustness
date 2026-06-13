# Builds the bundled example dataset `lee_example`.
#
# A small synthetic illustration shaped like the Lee (2008) incumbency RDD
# application: 7 specifications of an effect near 0.077 with correlated
# bootstrap draws. For documentation and tests only; NOT the real Lee draws,
# which are large, external, and produced by the Stata package.
#
# Run from the package root: Rscript data-raw/make_lee_example.R

set.seed(2026)

labels <- c("No controls", "Prev vote", "Office exp", "Elect exp",
            "All controls", "Resid outcome", "First diff")
theta  <- c(0.0766, 0.0777, 0.0766, 0.0765, 0.0775, 0.0808, 0.0788)
K      <- length(theta)
B      <- 9999
se     <- c(0.0114, 0.0108, 0.0114, 0.0114, 0.0109, 0.0136, 0.0134)

common <- rnorm(B, 0, 0.6)
draws  <- sapply(seq_len(K), function(k) {
  theta[k] + se[k] * (sqrt(0.7) * common + sqrt(0.3) * rnorm(B))
})
colnames(draws) <- labels

lee_example <- list(
  theta       = stats::setNames(theta, labels),
  draws       = draws,
  comparisons = list(
    all            = 1:7,
    main           = 1:5,
    main_resid     = 1:6,
    main_firstdiff = c(1, 2, 3, 4, 5, 7)
  )
)

save(lee_example, file = "data/lee_example.rda", compress = "xz")
