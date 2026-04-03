library(testthat)

# ---- adaboost ----

test_that("adaboost returns correct structure", {
  set.seed(1)
  X <- matrix(rnorm(200), 100, 2)
  y <- sign(X[, 1] + X[, 2])
  y[y == 0] <- 1
  res <- adaboost(X, y, M = 10)

  expect_type(res, "list")
  expect_length(res$models, 10)
  expect_length(res$alpha, 10)
  expect_equal(nrow(res$weight_history), 10)
  expect_equal(ncol(res$weight_history), 100)
})

test_that("adaboost predictions are in {-1, +1}", {
  set.seed(2)
  X <- matrix(rnorm(200), 100, 2)
  y <- sign(X[, 1] + X[, 2])
  y[y == 0] <- 1
  res <- adaboost(X, y, M = 20)
  preds <- predict_adaboost(res, X)

  expect_true(all(preds %in% c(-1, 1)))
})

test_that("adaboost achieves reasonable accuracy on separable data", {
  set.seed(3)
  X <- matrix(rnorm(400), 200, 2)
  y <- sign(X[, 1] + X[, 2])
  y[y == 0] <- 1
  res <- adaboost(X, y, M = 50)
  preds <- predict_adaboost(res, X)
  acc <- mean(preds == y)

  expect_gt(acc, 0.85)
})


# ---- acar ----

test_that("acar returns correct structure", {
  set.seed(10)
  res <- acar(c(10, 5, 3), n_ants = 50, n_waves = 30, early_stop = FALSE)

  expect_type(res, "list")
  expect_true("pheromone_history" %in% names(res))
  expect_true("final_decision" %in% names(res))
  expect_equal(ncol(res$pheromone_history), 3)
})

test_that("acar picks the best site", {
  set.seed(11)
  correct <- 0
  for (i in 1:20) {
    res <- acar(c(10, 3, 1), n_ants = 100, n_waves = 50)
    if (res$final_decision == 1) correct <- correct + 1
  }
  # Should pick site 1 most of the time

  expect_gt(correct, 14)
})


# ---- quorum margin ----

test_that("quorum margin is positive when best site has most pheromone", {
  set.seed(20)
  res <- acar(c(10, 3, 1), n_ants = 100, n_waves = 50)
  qm <- calculate_quorum_margin(res)
  expect_gt(qm, 0)
})


# ---- generate_classification_data ----

test_that("generate_classification_data returns correct dimensions", {
  dat <- generate_classification_data(n = 100, p = 5, noise = 0.1)
  expect_equal(nrow(dat$X), 100)
  expect_equal(ncol(dat$X), 5)
  expect_length(dat$y, 100)
  expect_true(all(dat$y %in% c(-1, 1)))
})
