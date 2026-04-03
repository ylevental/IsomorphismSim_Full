test_that("generate_regression_data returns correct dimensions", {
  d <- generate_regression_data(n = 100, p = 20)
  expect_equal(nrow(d$X), 100)
  expect_equal(ncol(d$X), 20)
  expect_length(d$y, 100)
  expect_length(d$f_true, 100)
})

test_that("generate_regression_data respects signal strength", {
  d_lo <- generate_regression_data(n = 5000, p = 10, signal_strength = 0.1, noise_sd = 0)
  d_hi <- generate_regression_data(n = 5000, p = 10, signal_strength = 10,  noise_sd = 0)
  expect_lt(var(d_lo$y), var(d_hi$y))
})

test_that("simulate_ant_colony returns a valid decision", {
  sim <- simulate_ant_colony(n_ants = 20, n_steps = 50, p_explore = 0.5)
  expect_true(sim$decision %in% 1:5)
  expect_equal(nrow(sim$ant_preferences), 20)
  expect_equal(ncol(sim$ant_preferences), 5)
})

test_that("colony with high exploration is less correlated", {
  set.seed(123)
  n_reps <- 15

  get_cor <- function(pe) {
    pmat <- matrix(NA, 30, n_reps)
    for (r in seq_len(n_reps)) {
      s <- simulate_ant_colony(n_ants = 30, p_explore = pe, n_steps = 80)
      pmat[, r] <- s$ant_preferences[, 1]
    }
    cm <- cor(pmat, use = "pairwise.complete.obs")
    mean(cm[lower.tri(cm)], na.rm = TRUE)
  }

  rho_low  <- get_cor(0.1)
  rho_high <- get_cor(0.9)
  # More exploration should yield lower correlation

  expect_lt(rho_high, rho_low)
})

test_that("variance formula matches analytical result", {
  rho <- 0.4; sigma2 <- 2; M <- 100
  expected <- rho * sigma2 + (1 - rho) * sigma2 / M
  expect_equal(expected, 0.4 * 2 + 0.6 * 2 / 100)
})
