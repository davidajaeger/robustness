test_that("R_obs equals the range of the estimates", {
  set.seed(1)
  theta <- c(0.10, 0.12, 0.11)
  draws <- sapply(theta, function(m) rnorm(2000, m, 0.02))
  r <- robustness(theta, draws)
  expect_equal(r$results$all$R, max(theta) - min(theta))
})

test_that("statistics match independent recomputation", {
  set.seed(7)
  theta <- c(0.10, 0.12, 0.11)
  B <- 5000; K <- 3
  draws <- sapply(theta, function(m) rnorm(B, m, 0.02))
  r <- robustness(theta, draws, alpha = c(0.50, 0.05))$results$all

  R_obs <- max(theta) - min(theta)
  R_unc <- apply(draws, 1, function(x) max(x) - min(x))
  Rstar_50 <- quantile(R_unc, 0.50, names = FALSE, type = 1)
  Rstar_95 <- quantile(R_unc, 0.95, names = FALSE, type = 1)
  expect_equal(r$equivalence$Rstar[1], Rstar_50)
  expect_equal(r$equivalence$Rstar[2], Rstar_95)

  thbar <- mean(theta)
  Dc <- sweep(draws, 2, theta, "-") + thbar
  R_rc <- apply(Dc, 1, function(x) max(x) - min(x))
  # The package uses the Davison-Hinkley (1 + #)/(B + 1) convention.
  expect_equal(r$p_R, (1 + sum(R_rc >= R_obs)) / (B + 1))

  Rmat <- matrix(-1 / K, K - 1, K); for (i in 1:(K - 1)) Rmat[i, i] <- Rmat[i, i] + 1
  V <- var(draws); RVRinv <- solve(Rmat %*% V %*% t(Rmat))
  d <- Rmat %*% theta
  expect_equal(r$W, as.numeric(t(d) %*% RVRinv %*% d))
})

test_that("named comparisons are honoured and labelled", {
  set.seed(2)
  theta <- c(0.2, 0.22, 0.19, 0.21)
  draws <- sapply(theta, function(m) rnorm(999, m, 0.03))
  r <- robustness(theta, draws,
                  comparisons = list(all = 1:4, extremes = c(1, 4)))
  expect_named(r$results, c("all", "extremes"))
  expect_equal(r$results$extremes$K, 2L)
  expect_equal(r$results$all$K, 4L)
})

test_that("single integer vector is treated as one comparison", {
  set.seed(3)
  theta <- c(0.2, 0.22, 0.19)
  draws <- sapply(theta, function(m) rnorm(500, m, 0.03))
  r <- robustness(theta, draws, comparisons = c(1, 3))
  expect_length(r$results, 1L)
  expect_equal(r$results[[1]]$K, 2L)
})

test_that("dimension mismatch errors", {
  expect_error(robustness(c(1, 2, 3), matrix(rnorm(20), 10, 2)),
               "one column per element")
})

test_that("fewer than two specs errors", {
  set.seed(4)
  draws <- matrix(rnorm(100), 50, 2)
  expect_error(robustness(c(0.1, 0.2), draws, comparisons = 1),
               "at least 2 specifications")
})

test_that("too many incomplete reps stops", {
  set.seed(5)
  theta <- c(0.1, 0.2)
  draws <- cbind(rnorm(100, 0.1, 0.02), rnorm(100, 0.2, 0.02))
  draws[1:10, 1] <- NA  # 10% incomplete
  expect_error(robustness(theta, draws, max_drop = 0.01),
               "exceeding max_drop")
})

test_that("incomplete reps within tolerance are dropped and counted", {
  set.seed(6)
  theta <- c(0.1, 0.2)
  draws <- cbind(rnorm(2000, 0.1, 0.02), rnorm(2000, 0.2, 0.02))
  draws[1, 1] <- NA  # 0.05% incomplete, under default 1%
  r <- robustness(theta, draws)
  expect_equal(r$results$all$B_dropped, 1L)
  expect_equal(r$results$all$B, 1999L)
})

test_that("p-values lie in [0, 1] and equivalence bounds are non-negative", {
  set.seed(8)
  theta <- c(0.05, 0.07, 0.06, 0.065)
  draws <- sapply(theta, function(m) rnorm(1500, m, 0.015))
  r <- robustness(theta, draws)$results$all
  expect_true(r$p_R >= 0 && r$p_R <= 1)
  expect_true(r$p_W >= 0 && r$p_W <= 1)
  expect_true(all(r$equivalence$Rstar >= 0))
  expect_true(all(r$equivalence$Wstar >= 0))
})

test_that("as.data.frame returns one row per comparison-by-alpha", {
  set.seed(9)
  theta <- c(0.2, 0.22, 0.19, 0.21)
  draws <- sapply(theta, function(m) rnorm(800, m, 0.03))
  r <- robustness(theta, draws,
                  comparisons = list(all = 1:4, extremes = c(1, 4)),
                  alpha = c(0.50, 0.05))
  df <- as.data.frame(r)
  expect_equal(nrow(df), 4L)  # 2 comparisons x 2 alphas
  expect_true(all(c("comparison", "R", "p_R", "alpha", "Rstar") %in% names(df)))
})

test_that("per-spec sample sizes are averaged and reported", {
  set.seed(10)
  theta <- c(0.1, 0.2, 0.15)
  draws <- sapply(theta, function(m) rnorm(1000, m, 0.02))
  nmat <- cbind(rep(5000, 1000), rep(4000, 1000), rep(4500, 1000))
  r <- robustness(theta, draws, n_boot = nmat)$results$all
  expect_equal(r$avg_n, c(5000, 4000, 4500))
})

test_that("non-finite draws error", {
  set.seed(11)
  theta <- c(0.1, 0.2)
  draws <- matrix(rnorm(200), 100, 2)
  draws[1, 1] <- Inf
  expect_error(robustness(theta, draws), "Inf")
})

test_that("non-integer comparison indices error", {
  set.seed(12)
  theta <- c(0.1, 0.2, 0.3)
  draws <- matrix(rnorm(300), 100, 3)
  # Bare vector path: a single comparison given as a numeric vector.
  expect_error(robustness(theta, draws, comparisons = c(1, 2.5)),
               "whole numbers")
  # Named-list path: same check should fire inside the list.
  expect_error(robustness(theta, draws,
                          comparisons = list(mix = c(1, 2.5))),
               "whole numbers")
})
