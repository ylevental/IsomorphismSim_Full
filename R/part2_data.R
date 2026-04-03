#' @title Generate Synthetic Classification Data
#' @description Creates a binary classification problem where the true boundary
#'   depends on the first three features: \eqn{y = sign(x_1 + x_2 x_3 + 0.5 x_3^2)},
#'   with optional label noise.
#' @param n Number of observations.
#' @param p Number of features.
#' @param noise Fraction of labels flipped.
#' @return A list with \code{X} (n x p matrix), \code{y} (label vector), and
#'   \code{true_labels}.
#' @export
generate_classification_data <- function(n = 500, p = 10, noise = 0.2) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  true_labels <- sign(X[, 1] + X[, 2] * X[, 3] + 0.5 * X[, 3]^2)
  true_labels[true_labels == 0] <- 1
  y <- true_labels
  if (noise > 0) {
    flip_idx <- sample(seq_len(n), size = floor(n * noise))
    y[flip_idx] <- -y[flip_idx]
  }
  list(X = X, y = y, true_labels = true_labels)
}
