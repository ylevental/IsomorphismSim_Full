#' Generational Ant Colony Learning (GACL)
#'
#' Implements the full GACL algorithm with configurable parameters.
#' Pheromone evolution across generations follows update equations
#' isomorphic to stochastic gradient descent.
#'
#' @param site_qualities Numeric vector of true qualities for K sites.
#' @param n_ants Number of ants per generation (default 100).
#' @param n_generations Number of generations to simulate (default 50).
#' @param n_waves Number of recruitment waves per generation (default 10).
#' @param rho_wave Within-generation evaporation rate (default 0.3).
#' @param rho_gen Between-generation evaporation rate, analogous to the
#'   learning rate \eqn{\eta} in neural networks (default 0.1).
#' @param gamma Pheromone deposition rate (default 0.5).
#' @param alpha Pheromone influence exponent (default 1).
#' @param beta Heuristic influence exponent (default 1).
#' @param noise_sd Standard deviation of observation noise (default 0.2).
#'
#' @return A list with components:
#' \describe{
#'   \item{pheromone_history}{Matrix (generations x K) of pheromone levels.}
#'   \item{fitness_history}{Numeric vector of colony fitness per generation.}
#'   \item{error_signal_history}{Numeric vector of negative fitness (loss analog).}
#'   \item{decision_history}{Integer vector; best site per generation.}
#'   \item{final_decision}{Best site at the last generation.}
#' }
#'
#' @examples
#' result <- gacl(c(10, 7, 5, 4, 3), n_generations = 30)
#' plot(result$fitness_history, type = "l")
#'
#' @export
gacl <- function(site_qualities,
                 n_ants = 100,
                 n_generations = 50,
                 n_waves = 10,
                 rho_wave = 0.3,
                 rho_gen = 0.1,
                 gamma = 0.5,
                 alpha = 1,
                 beta = 1,
                 noise_sd = 0.2) {

  K <- length(site_qualities)
  heuristic <- site_qualities / max(site_qualities)

  pheromone_history    <- matrix(NA, nrow = n_generations, ncol = K)
  fitness_history      <- numeric(n_generations)
  error_signal_history <- numeric(n_generations)
  decision_history     <- numeric(n_generations)

  tau <- rep(1, K)
  pheromone_history[1, ] <- tau

  for (gen in 1:n_generations) {
    tau_wave <- tau
    for (wave in 1:n_waves) {
      prob <- (tau_wave^alpha) * (heuristic^beta)
      prob <- prob / sum(prob)
      visits <- stats::rmultinom(1, n_ants, prob)[, 1]
      observed_qualities <- pmax(site_qualities + stats::rnorm(K, 0, noise_sd), 0)
      tau_wave <- (1 - rho_wave) * tau_wave + gamma * visits * observed_qualities
      tau_wave <- pmax(tau_wave, 0.01)
    }
    tau <- tau_wave
    pheromone_history[gen, ] <- tau

    fitness <- sum(visits * observed_qualities) / n_ants
    fitness_history[gen]      <- fitness
    error_signal_history[gen] <- -fitness
    decision_history[gen]     <- which.max(tau)

    if (gen < n_generations) {
      gradient <- (fitness - mean(fitness_history[max(1, gen - 5):gen])) *
                  (tau / sum(tau))
      tau <- (1 - rho_gen) * tau + gamma * gradient
      tau <- pmax(tau, 0.01)
    }
  }

  list(pheromone_history    = pheromone_history,
       fitness_history      = fitness_history,
       error_signal_history = error_signal_history,
       decision_history     = decision_history,
       final_decision       = decision_history[n_generations])
}


#' Simple Neural Network with Stochastic Gradient Descent
#'
#' Implements a single-hidden-layer perceptron trained with mini-batch SGD
#' and binary cross-entropy loss.  Under the isomorphism, weight updates
#' correspond to pheromone updates in [gacl()].
#'
#' @param X Numeric matrix of input features (n x p).
#' @param y Numeric vector of labels (\{-1, 1\} or \{0, 1\}).
#' @param n_epochs Number of training epochs (default 50).
#' @param batch_size Mini-batch size (default 32).
#' @param learning_rate Learning rate \eqn{\eta}, analogous to the
#'   evaporation rate \eqn{\rho} in GACL (default 0.1).
#' @param n_hidden Number of hidden units (default 10).
#' @param validation_split Fraction of data held out for validation (default 0.2).
#'
#' @return A list with components:
#' \describe{
#'   \item{train_loss}{Numeric vector of training loss per epoch.}
#'   \item{val_loss}{Numeric vector of validation loss per epoch.}
#'   \item{val_acc}{Numeric vector of validation accuracy per epoch.}
#'   \item{gradient_norm}{Numeric vector of gradient L2 norm per epoch.}
#'   \item{final_weights}{List of weight matrices (W1, b1, W2, b2).}
#' }
#'
#' @examples
#' d <- generate_synthetic_data(n = 300, p = 5)
#' nn <- simple_neural_network(d$X[1:240, ], d$y[1:240], n_epochs = 30)
#' plot(nn$val_acc, type = "l")
#'
#' @export
simple_neural_network <- function(X, y,
                                  n_epochs = 50,
                                  batch_size = 32,
                                  learning_rate = 0.1,
                                  n_hidden = 10,
                                  validation_split = 0.2) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (all(unique(y) %in% c(-1, 1))) y <- ifelse(y == -1, 0, 1)

  n <- nrow(X); p <- ncol(X)
  n_val <- floor(n * validation_split)
  idx_val   <- sample(1:n, n_val)
  idx_train <- setdiff(1:n, idx_val)
  X_train <- X[idx_train, ]; y_train <- y[idx_train]
  X_val   <- X[idx_val, ];   y_val   <- y[idx_val]

  # Xavier initialisation

  W1 <- matrix(stats::rnorm(p * n_hidden, 0, sqrt(2 / p)), p, n_hidden)
  b1 <- matrix(0, 1, n_hidden)
  W2 <- matrix(stats::rnorm(n_hidden, 0, sqrt(2 / n_hidden)), n_hidden, 1)
  b2 <- 0

  train_loss_history    <- numeric(n_epochs)
  val_loss_history      <- numeric(n_epochs)
  val_acc_history       <- numeric(n_epochs)
  gradient_norm_history <- numeric(n_epochs)

  sigmoid <- function(x) 1 / (1 + exp(-pmax(pmin(x, 50), -50)))

  for (epoch in 1:n_epochs) {
    shuf <- sample(seq_along(y_train))
    X_train <- X_train[shuf, ]; y_train <- y_train[shuf]
    epoch_loss <- 0

    for (i in seq(1, length(y_train), batch_size)) {
      end_idx <- min(i + batch_size - 1, length(y_train))
      Xb <- X_train[i:end_idx, , drop = FALSE]; yb <- y_train[i:end_idx]
      bn <- nrow(Xb)

      z1 <- Xb %*% W1 + matrix(1, bn, 1) %*% b1
      a1 <- sigmoid(z1)
      z2 <- a1 %*% W2 + b2
      a2 <- pmin(pmax(sigmoid(z2), 1e-6), 1 - 1e-6)

      loss <- -mean(yb * log(a2) + (1 - yb) * log(1 - a2))
      epoch_loss <- epoch_loss + loss * bn

      d2  <- (a2 - yb) / bn
      dW2 <- t(a1) %*% d2; db2 <- sum(d2)
      d1  <- (d2 %*% t(W2)) * (a1 * (1 - a1))
      dW1 <- t(Xb) %*% d1;  db1 <- colSums(d1)

      mx <- 5
      dW1 <- pmin(pmax(dW1, -mx), mx)
      dW2 <- pmin(pmax(dW2, -mx), mx)

      W1 <- W1 - learning_rate * dW1; b1 <- b1 - learning_rate * db1
      W2 <- W2 - learning_rate * dW2; b2 <- b2 - learning_rate * db2
    }

    train_loss_history[epoch]    <- epoch_loss / length(y_train)
    gradient_norm_history[epoch] <- sqrt(sum(dW1^2) + sum(dW2^2))

    z1v <- X_val %*% W1 + matrix(1, n_val, 1) %*% b1
    a1v <- sigmoid(z1v)
    z2v <- a1v %*% W2 + b2
    a2v <- pmin(pmax(sigmoid(z2v), 1e-6), 1 - 1e-6)
    val_loss_history[epoch] <- -mean(y_val * log(a2v) + (1 - y_val) * log(1 - a2v))
    val_acc_history[epoch]  <- mean((a2v > 0.5) == y_val)
  }

  list(train_loss    = train_loss_history,
       val_loss      = val_loss_history,
       val_acc       = val_acc_history,
       gradient_norm = gradient_norm_history,
       final_weights = list(W1 = W1, b1 = b1, W2 = W2, b2 = b2))
}


#' Generate Synthetic Classification Data
#'
#' Creates a binary classification dataset with a non-linear decision
#' boundary whose complexity can be varied.
#'
#' @param n Number of samples (default 1000).
#' @param p Number of features (default 5).
#' @param noise Label-flip noise rate, between 0 and 0.5 (default 0.1).
#' @param complexity Complexity of the boundary: 1 = linear, 2 = quadratic,
#'   3 = complex with interactions (default 2).
#'
#' @return A list with components \code{X} (matrix), \code{y} (labels with
#'   noise), and \code{true_labels}.
#'
#' @examples
#' d <- generate_synthetic_data(n = 500, p = 5, complexity = 2)
#' table(d$y)
#'
#' @export
generate_synthetic_data <- function(n = 1000, p = 5, noise = 0.1, complexity = 2) {
  X <- matrix(stats::rnorm(n * p), n, p)
  if (complexity == 1) {
    true_labels <- sign(X[, 1] + X[, 2])
  } else if (complexity == 2) {
    true_labels <- sign(X[, 1]^2 + X[, 2]^2 - 1)
  } else {
    true_labels <- sign(X[, 1] + X[, 2] * X[, 3] + sin(X[, 4]) + X[, 5]^2)
  }
  y <- true_labels
  flip_idx <- sample(1:n, floor(n * noise))
  y[flip_idx] <- -y[flip_idx]
  list(X = X, y = y, true_labels = true_labels)
}
