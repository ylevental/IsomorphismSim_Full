#' Experiment 1: Random Forest Variance Decomposition
#'
#' Validates that \eqn{\mathrm{Var}[\hat{f}_{\mathrm{rf}}] =
#' \rho\sigma^2 + (1-\rho)\sigma^2/M}{Var[f_rf] = rho*sigma^2 + (1-rho)*sigma^2/M}
#' holds empirically and that \eqn{\rho} varies with \code{m_try/p}.
#'
#' @param n_train Training set size
#' @param n_test Test set size
#' @param p Number of features
#' @param m_try_values Integer vector of \code{mtry} values to sweep
#' @param n_trees Number of trees per forest
#' @param n_replicates Number of Monte Carlo replicates
#' @param signal_strength Signal multiplier for \code{\link{generate_regression_data}}
#' @param noise_sd Noise level for \code{\link{generate_regression_data}}
#' @return A data.frame with columns \code{m_try, rep, rho, sigma2,
#'   ensemble_var, pred_error, theoretical_var}
#' @export
variance_decomposition_experiment <- function(
    n_train      = 500,
    n_test       = 1000,
    p            = 50,
    m_try_values = c(1, 2, 5, 10, 20, 50),
    n_trees      = 500,
    n_replicates = 100,
    signal_strength = 2,
    noise_sd     = 1
) {
  results <- expand.grid(m_try = m_try_values, rep = seq_len(n_replicates))
  results$rho <- results$sigma2 <- results$ensemble_var <- NA_real_
  results$pred_error <- results$theoretical_var <- NA_real_

  for (i in seq_len(nrow(results))) {
    data <- generate_regression_data(n_train + n_test, p, signal_strength, noise_sd)
    X_train     <- data$X[seq_len(n_train), ]
    y_train     <- data$y[seq_len(n_train)]
    X_test      <- data$X[(n_train + 1):(n_train + n_test), ]
    f_true_test <- data$f_true[(n_train + 1):(n_train + n_test)]

    rf <- ranger::ranger(
      y ~ ., data = cbind(y = y_train, X_train),
      num.trees = n_trees, mtry = results$m_try[i],
      keep.inbag = TRUE, write.forest = TRUE
    )

    tree_preds <- predict(rf, data = X_test, predict.all = TRUE)$predictions
    sigma2_est <- mean(apply(tree_preds, 2, stats::var))

    # Subsample for speed
    k <- min(ncol(tree_preds), 100)
    idx <- sample(ncol(tree_preds), k)
    cm <- stats::cor(tree_preds[, idx])
    rho_est <- mean(cm[lower.tri(cm)])

    ens <- rowMeans(tree_preds)
    results$sigma2[i]          <- sigma2_est
    results$rho[i]             <- rho_est
    results$ensemble_var[i]    <- stats::var(ens)
    results$theoretical_var[i] <- rho_est * sigma2_est +
      (1 - rho_est) * sigma2_est / n_trees
    results$pred_error[i]      <- mean((ens - f_true_test)^2)

    if (i %% 50 == 0) message("  RF exp: ", i, "/", nrow(results))
  }
  results
}

#' Experiment 2: Ant Colony Variance Decomposition
#'
#' Sweeps colony size and exploration probability, computing decision
#' accuracy and the empirical variance/correlation of ant assessments.
#'
#' @param n_ants_values Integer vector of colony sizes
#' @param p_explore_values Numeric vector of exploration probabilities
#' @param n_replicates Monte Carlo replicates
#' @param n_sites Number of candidate nest sites
#' @param site_qualities True site quality vector
#' @return A data.frame with columns \code{n_ants, p_explore, accuracy_mean,
#'   accuracy_sd, var_mean, cor_mean}
#' @export
colony_variance_experiment <- function(
    n_ants_values    = c(10, 20, 30, 50, 100),
    p_explore_values = seq(0, 1, by = 0.2),
    n_replicates     = 50,
    n_sites          = 5,
    site_qualities   = c(10, 8, 6, 4, 2)
) {
  best_site <- which.max(site_qualities)
  results <- expand.grid(n_ants = n_ants_values, p_explore = p_explore_values,
                          stringsAsFactors = FALSE)
  results$accuracy_mean <- results$accuracy_sd <- NA_real_
  results$var_mean <- results$cor_mean <- NA_real_

  for (i in seq_len(nrow(results))) {
    decs <- integer(n_replicates)
    rhos <- numeric(n_replicates)
    vars <- numeric(n_replicates)
    for (r in seq_len(n_replicates)) {
      sim <- simulate_ant_colony(
        n_ants = results$n_ants[i], p_explore = results$p_explore[i],
        n_sites = n_sites, site_qualities = site_qualities
      )
      decs[r] <- as.integer(sim$decision == best_site)
      # CORRECTED: within-colony correlation (ants' preference vectors over sites)
      rhos[r] <- within_colony_correlation(sim$ant_preferences)
      vars[r] <- stats::var(sim$ant_preferences[, best_site], na.rm = TRUE)
    }
    results$accuracy_mean[i] <- mean(decs)
    results$accuracy_sd[i]   <- stats::sd(decs)
    results$var_mean[i]      <- mean(vars, na.rm = TRUE)
    results$cor_mean[i]      <- mean(rhos, na.rm = TRUE)
    if (i %% 10 == 0) message("  Colony exp: ", i, "/", nrow(results))
  }
  results
}

#' Experiment 3: Direct Isomorphism Test
#'
#' Runs random forest and ant colony under matched \eqn{\theta} values
#' (\code{m_try/p = p_explore}) and compares the resulting pairwise
#' correlation and ensemble variance.
#'
#' @param n_replicates Monte Carlo replicates per \eqn{\theta}
#' @param p_vals Numeric vector of \eqn{\theta} values
#' @return A data.frame with \code{theta, system, correlation, cor_sd,
#'   ensemble_var, var_sd}
#' @export
isomorphism_test <- function(
    n_replicates = 50,
    p_vals       = c(0.1, 0.3, 0.5, 0.7, 0.9)
) {
  p <- 50
  results <- data.frame()
  for (theta in p_vals) {
    rf_cors <- rf_vars <- numeric(n_replicates)
    for (r in seq_len(n_replicates)) {
      d <- generate_regression_data(n = 500, p = p)
      rf <- ranger::ranger(
        y ~ ., data = cbind(y = d$y[1:400], d$X[1:400, ]),
        num.trees = 200, mtry = max(1, round(theta * p))
      )
      pa <- predict(rf, data = d$X[401:500, ], predict.all = TRUE)$predictions
      cm <- stats::cor(pa)
      rf_cors[r] <- mean(cm[lower.tri(cm)])
      rf_vars[r] <- stats::var(rowMeans(pa))
    }
    results <- rbind(results, data.frame(
      theta = theta, system = "Random Forest",
      correlation = mean(rf_cors), cor_sd = stats::sd(rf_cors),
      ensemble_var = mean(rf_vars), var_sd = stats::sd(rf_vars)
    ))

    apm_rhos <- numeric(n_replicates)
    ant_vars <- numeric(n_replicates)
    sq <- seq(10, 2, length.out = 20)  # 20 sites, quality gap = 8
    for (r in seq_len(n_replicates)) {
      sim <- simulate_ant_colony(n_ants = 50, p_explore = theta,
                                  n_sites = 20, site_qualities = sq,
                                  n_steps = 50)
      # CORRECTED: within-colony correlation
      apm_rhos[r] <- within_colony_correlation(sim$ant_preferences)
      ant_vars[r] <- stats::var(sim$ant_preferences[, which.max(sq)], na.rm = TRUE)
    }
    results <- rbind(results, data.frame(
      theta = theta, system = "Ant Colony",
      correlation = mean(apm_rhos, na.rm = TRUE),
      cor_sd = stats::sd(apm_rhos, na.rm = TRUE),
      ensemble_var = mean(ant_vars, na.rm = TRUE), var_sd = NA
    ))
  }
  results
}

#' Experiment 4: Optimal Decorrelation
#'
#' Sweeps ensemble size and \eqn{\theta} to find the performance-optimal
#' decorrelation parameter. Uses MSE for random forests and 1−accuracy
#' for ant colonies so both metrics point downward.
#'
#' @param ensemble_sizes Integer vector of ensemble sizes
#' @param theta_values Numeric vector of \eqn{\theta} values
#' @param n_replicates Monte Carlo replicates
#' @return A data.frame with \code{M, theta, system, performance, se}
#' @export
optimal_decorrelation_experiment <- function(
    ensemble_sizes = c(10, 30, 100),
    theta_values   = seq(0.1, 0.9, by = 0.1),
    n_replicates   = 30
) {
  p <- 50
  results <- data.frame()
  for (M in ensemble_sizes) {
    for (theta in theta_values) {
      rf_err <- replicate(n_replicates, {
        d <- generate_regression_data(n = 500, p = p)
        rf <- ranger::ranger(y ~ ., data = cbind(y = d$y[1:400], d$X[1:400, ]),
                              num.trees = M, mtry = max(1, round(theta * p)))
        mean((predict(rf, data = d$X[401:500, ])$predictions - d$f_true[401:500])^2)
      })
      ant_acc <- replicate(n_replicates, {
        sim <- simulate_ant_colony(n_ants = M, p_explore = theta)
        as.numeric(sim$decision == which.max(c(10, 8, 6, 4, 2)))
      })
      results <- rbind(results, data.frame(
        M = M, theta = theta, system = "Random Forest",
        performance = mean(rf_err), se = stats::sd(rf_err) / sqrt(n_replicates)
      ), data.frame(
        M = M, theta = theta, system = "Ant Colony",
        performance = 1 - mean(ant_acc), se = stats::sd(ant_acc) / sqrt(n_replicates)
      ))
    }
  }
  results
}

#' Experiment 5: Sensitivity Analysis
#'
#' Tests whether the isomorphism holds under varying signal strengths,
#' noise levels, and \eqn{\theta} values.
#'
#' @param signal_strengths Numeric vector of signal multipliers
#' @param noise_levels Numeric vector of noise standard deviations
#' @param n_replicates Monte Carlo replicates
#' @return A data.frame with \code{signal, noise, theta, system,
#'   correlation, cor_sd}
#' @export
sensitivity_analysis <- function(
    signal_strengths = c(0.5, 1, 2, 5),
    noise_levels     = c(0.5, 1, 2, 4),
    n_replicates     = 20
) {
  p <- 50
  results <- data.frame()
  for (signal in signal_strengths) {
    for (noise in noise_levels) {
      for (theta in c(0.2, 0.5, 0.8)) {
        rf_cors <- replicate(n_replicates, {
          d <- generate_regression_data(500, p, signal, noise)
          rf <- ranger::ranger(y ~ ., data = cbind(y = d$y[1:400], d$X[1:400, ]),
                                num.trees = 100, mtry = max(1, round(theta * p)))
          pa <- predict(rf, d$X[401:500, ], predict.all = TRUE)$predictions
          cm <- stats::cor(pa); mean(cm[lower.tri(cm)])
        })

        sq <- c(10, 8, 6, 4, 2) * signal
        apm <- matrix(NA_real_, 50, n_replicates)
        for (r in seq_len(n_replicates)) {
          sim <- simulate_ant_colony(n_ants = 50, p_explore = theta,
                                      site_qualities = sq, noise_sd = noise)
          apm[, r] <- sim$ant_preferences[, which.max(sq)]
        }
        acm <- stats::cor(apm, use = "pairwise.complete.obs")

        results <- rbind(results, data.frame(
          signal = signal, noise = noise, theta = theta,
          system = "Random Forest",
          correlation = mean(rf_cors, na.rm = TRUE),
          cor_sd = stats::sd(rf_cors, na.rm = TRUE)
        ), data.frame(
          signal = signal, noise = noise, theta = theta,
          system = "Ant Colony",
          correlation = mean(acm[lower.tri(acm)], na.rm = TRUE),
          cor_sd = NA
        ))
      }
    }
  }
  results
}
