#' @title AdaBoost with Decision Stumps
#' @description Implements the AdaBoost algorithm using decision stumps as weak
#'   learners, following Freund & Schapire (1997).
#' @param X Numeric matrix of predictors (n x p).
#' @param y Numeric vector of labels in \{-1, +1\}.
#' @param M Integer number of boosting iterations.
#' @return A list with components \code{models} (list of stumps),
#'   \code{alpha} (learner weights), and \code{weight_history} (M x n matrix
#'   of instance weights at each iteration).
#' @export
adaboost <- function(X, y, M = 100) {
  n <- length(y)
  D <- rep(1 / n, n)
  models <- vector("list", M)
  alpha <- numeric(M)
  weight_history <- matrix(NA_real_, nrow = M, ncol = n)

  for (m in seq_len(M)) {
    weight_history[m, ] <- D
    best <- find_best_stump(X, y, D)
    predictions <- predict_stump(best, X)
    error <- sum(D * (predictions != y))
    error <- max(error, 1e-10)
    error <- min(error, 1 - 1e-10)

    alpha[m] <- 0.5 * log((1 - error) / error)
    D <- D * exp(-alpha[m] * y * predictions)
    D <- D / sum(D)

    models[[m]] <- best
  }

  list(models = models, alpha = alpha, weight_history = weight_history)
}


#' @title Find the Best Decision Stump
#' @description Searches over all features and thresholds to find the stump
#'   that minimizes the weighted classification error.
#' @param X Numeric matrix of predictors.
#' @param y Label vector in \{-1, +1\}.
#' @param D Weight vector (sums to 1).
#' @return A list describing the best stump: \code{feature}, \code{threshold},
#'   and \code{direction}.
#' @export
find_best_stump <- function(X, y, D) {
  n <- nrow(X)
  p <- ncol(X)
  best_error <- Inf
  best_stump <- NULL

  for (j in seq_len(p)) {
    thresholds <- sort(unique(X[, j]))
    # Subsample thresholds for speed if many unique values
    if (length(thresholds) > 50) {
      thresholds <- quantile(X[, j], probs = seq(0, 1, length.out = 50))
    }
    for (thresh in thresholds) {
      left_pred  <- ifelse(X[, j] <= thresh,  1, -1)
      error_left <- sum(D * (left_pred != y))

      right_pred  <- ifelse(X[, j] <= thresh, -1,  1)
      error_right <- sum(D * (right_pred != y))

      if (error_left < best_error) {
        best_error <- error_left
        best_stump <- list(feature = j, threshold = thresh, direction = "left")
      }
      if (error_right < best_error) {
        best_error <- error_right
        best_stump <- list(feature = j, threshold = thresh, direction = "right")
      }
    }
  }
  best_stump
}


#' @title Predict with a Decision Stump
#' @param stump A stump list from \code{find_best_stump}.
#' @param X Numeric predictor matrix.
#' @return Numeric prediction vector in \{-1, +1\}.
#' @export
predict_stump <- function(stump, X) {
  if (stump$direction == "left") {
    ifelse(X[, stump$feature] <= stump$threshold, 1, -1)
  } else {
    ifelse(X[, stump$feature] <= stump$threshold, -1, 1)
  }
}


#' @title Predict with an AdaBoost Ensemble
#' @param result Output of \code{adaboost()}.
#' @param X Numeric predictor matrix.
#' @param up_to Optional integer; use only the first \code{up_to} iterations.
#' @return Numeric prediction vector in \{-1, +1\}.
#' @export
predict_adaboost <- function(result, X, up_to = NULL) {
  if (is.null(up_to)) up_to <- length(result$models)
  pred_matrix <- matrix(0, nrow = nrow(X), ncol = up_to)
  for (m in seq_len(up_to)) {
    pred_matrix[, m] <- predict_stump(result$models[[m]], X)
  }
  sign(rowSums(result$alpha[1:up_to] * pred_matrix))
}


#' @title Track AdaBoost Weight Evolution
#' @description Convenience wrapper that returns only the weight history.
#' @param X Predictor matrix.
#' @param y Label vector.
#' @param M Number of iterations.
#' @return An M x n matrix of instance weights.
#' @export
track_weights <- function(X, y, M = 50) {
  res <- adaboost(X, y, M = M)
  res$weight_history
}


#' @title Calculate Boosting Margins
#' @description Computes the normalized margin of each training instance.
#' @param result Output of \code{adaboost()}.
#' @param X Predictor matrix.
#' @param y Label vector.
#' @return Numeric vector of margins.
#' @export
calculate_margins <- function(result, X, y) {
  M <- length(result$models)
  pred_matrix <- matrix(0, nrow = nrow(X), ncol = M)
  for (m in seq_len(M)) {
    pred_matrix[, m] <- predict_stump(result$models[[m]], X)
  }
  y * rowSums(result$alpha * pred_matrix) / sum(abs(result$alpha))
}
