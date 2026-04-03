#' @title Weak Learnability Experiment
#' @description Tests whether gamma-weak ant colonies converge to correct
#'   decisions with sufficient recruitment waves, as predicted by Theorem 4.
#'   Parameters are calibrated so that small gamma is genuinely hard (slow
#'   convergence) while large gamma converges quickly.
#' @param gamma_values Numeric vector of weakness parameters.
#' @param wave_counts Integer vector of wave counts to test.
#' @param n_replicates Number of Monte Carlo replicates per combination.
#' @return A data.frame with columns \code{gamma}, \code{waves}, \code{rep},
#'   and \code{accuracy}.
#' @export
weak_learnability_experiment <- function(
    gamma_values = c(0.05, 0.1, 0.2, 0.3),
    wave_counts  = c(3, 5, 8, 12, 18, 25, 40, 60, 90, 130, 200, 300),
    n_replicates = 200) {

  results <- expand.grid(
    gamma = gamma_values,
    waves = wave_counts,
    rep   = seq_len(n_replicates)
  )
  results$accuracy <- NA_real_

  cat("  Figure 1: ", nrow(results), " total simulations\n")
  for (i in seq_len(nrow(results))) {
    if (i %% 1000 == 0) cat(sprintf("\r  Progress: %d / %d", i, nrow(results)))
    gv <- results$gamma[i]
    wc <- results$waves[i]
    # Quality gap = 5 * gamma (ranges from 0.25 to 1.5)
    sq <- c(10, 10 - 5 * gv)
    # n_ants=5: small committee per wave (like a shallow decision stump).
    # noise_sd=3: each ant's observation is noisy but 5 ants per wave
    # give a weak aggregate signal. Per-wave P(majority correct):
    #   gamma=0.05 (gap=0.25): per-ant ~53.3% -> 5-ant majority ~54%
    #   gamma=0.10 (gap=0.50): per-ant ~56.6% -> 5-ant majority ~60%
    #   gamma=0.20 (gap=1.00): per-ant ~63.1% -> 5-ant majority ~71%
    #   gamma=0.30 (gap=1.50): per-ant ~69.1% -> 5-ant majority ~81%
    # Pheromone feedback amplifies these over waves, giving clearly
    # separated convergence curves.
    # rho=0.02: slow evaporation so signal accumulates across waves.
    # beta=0.5: mild heuristic helps establish initial direction.
    sim <- acar(sq, n_ants = 5, n_waves = wc, noise_sd = 3.0,
                rho = 0.02, gamma = 0.5, alpha = 1, beta = 0.5,
                early_stop = FALSE)
    results$accuracy[i] <- as.numeric(sim$final_decision == 1)
  }
  cat("\n")
  results
}


#' @title Convergence Rate Experiment
#' @description Compares convergence of AdaBoost and ACAR across iterations.
#'   For AdaBoost: test-set accuracy using first m iterations.
#'   For ACAR: P(correct final decision) using m waves, averaged across
#'   independent replicates.
#' @param n_replicates Number of Monte Carlo replicates.
#' @param max_iters Maximum iterations / waves.
#' @return A data.frame with columns \code{iteration}, \code{accuracy},
#'   \code{system}, and \code{rep}.
#' @export
convergence_experiment_boost <- function(n_boost_reps = 30, n_acar_reps = 200,
                                   max_iters = 80) {
  results <- data.frame()

  # Measure at these wave/iteration counts
  measure_at <- sort(unique(c(
    seq(1, 10, by = 1),
    seq(12, 30, by = 2),
    seq(35, max_iters, by = 5)
  )))
  measure_at <- measure_at[measure_at <= max_iters]

  # --- AdaBoost: continuous accuracy, fewer reps needed ---
  cat("  AdaBoost convergence (", n_boost_reps, "reps)...\n")
  for (rep in seq_len(n_boost_reps)) {
    cat(sprintf("\r  AdaBoost replicate %d / %d", rep, n_boost_reps))
    set.seed(rep * 137)
    data_tr <- generate_classification_data(n = 200, p = 5, noise = 0.1)
    data_te <- generate_classification_data(n = 500, p = 5, noise = 0.0)

    boost_res <- adaboost(data_tr$X, data_tr$y, M = max_iters)
    for (m in measure_at) {
      pred <- predict_adaboost(boost_res, data_te$X, up_to = m)
      results <- rbind(results,
        data.frame(iteration = m,
                   accuracy  = mean(pred == data_te$y),
                   system    = "AdaBoost",
                   rep       = rep))
    }
  }

  # --- ACAR: binary outcomes, many reps needed for smooth curve ---
  cat("\n  ACAR convergence (", n_acar_reps, "reps)...\n")
  sq <- c(10, 9.5)  # gap = 0.5
  for (rep in seq_len(n_acar_reps)) {
    if (rep %% 20 == 0) cat(sprintf("\r  ACAR replicate %d / %d", rep, n_acar_reps))
    set.seed(rep * 251)
    for (m in measure_at) {
      sim <- acar(sq, n_ants = 5, n_waves = m, noise_sd = 5.0,
                  rho = 0.02, gamma = 0.5, alpha = 1, beta = 0.5,
                  early_stop = FALSE)
      correct <- as.numeric(sim$final_decision == 1)
      results <- rbind(results,
        data.frame(iteration = m,
                   accuracy  = correct,
                   system    = "ACAR",
                   rep       = rep))
    }
  }
  cat("\n")
  results
}


#' @title Noise Robustness Experiment
#' @description Compares how AdaBoost and ACAR degrade under increasing noise.
#'   Both systems are measured as P(correct) across independent replicates.
#'   Parameters are calibrated so both show visible, parallel degradation.
#' @param noise_levels Numeric vector of noise fractions (0 to 0.45).
#' @param n_replicates Number of Monte Carlo replicates per level.
#' @return A data.frame with columns \code{noise}, \code{accuracy}, and
#'   \code{system}.
#' @export
noise_experiment_boost <- function(noise_levels = seq(0, 0.45, by = 0.05),
                             n_replicates = 50) {
  results <- data.frame()

  for (noise in noise_levels) {
    cat(sprintf("\r  Noise level %.2f", noise))
    for (rep in seq_len(n_replicates)) {
      set.seed(rep * 257 + which(noise_levels == noise) * 31)

      # --- AdaBoost ---
      data_tr <- generate_classification_data(n = 200, p = 5, noise = noise)
      data_te <- generate_classification_data(n = 500, p = 5, noise = 0.0)
      boost_res <- adaboost(data_tr$X, data_tr$y, M = 50)
      pred      <- predict_adaboost(boost_res, data_te$X)
      boost_acc <- mean(pred == data_te$y)

      # --- ACAR ---
      # Binary choice, tight gap = 1.  Few ants + few waves so noise
      # actually degrades performance.
      # Calibration:
      #   noise=0.0 -> sd=0.01 -> trivial        -> ~100%
      #   noise=0.2 -> sd=2.4  -> moderate        -> ~80-90%
      #   noise=0.4 -> sd=4.8  -> hard (sd>>gap)  -> ~55-65%
      sq <- c(10, 9)  # gap = 1
      ant_noise_sd <- max(noise * 12, 0.01)
      sim <- acar(sq, n_ants = 8, n_waves = 15,
                  noise_sd = ant_noise_sd, early_stop = FALSE)
      ant_acc <- as.numeric(sim$final_decision == 1)

      results <- rbind(results,
        data.frame(noise = noise, accuracy = boost_acc, system = "AdaBoost"),
        data.frame(noise = noise, accuracy = ant_acc,   system = "ACAR")
      )
    }
  }
  cat("\n")
  results
}
