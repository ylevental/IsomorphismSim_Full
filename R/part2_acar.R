#' @title Ant Colony Adaptive Recruitment (ACAR)
#' @description Simulates the adaptive recruitment dynamics of an ant colony,
#'   implementing Algorithm 3 (ACAR) from the Part II manuscript.
#' @param site_qualities Numeric vector of true site qualities.
#' @param n_ants Number of ants released per wave.
#' @param n_waves Maximum number of recruitment waves.
#' @param rho Evaporation rate in \eqn{(0,1]}.
#' @param gamma Pheromone deposition rate.
#' @param alpha Pheromone influence exponent.
#' @param beta Heuristic influence exponent.
#' @param noise_sd Standard deviation of observation noise.
#' @param early_stop Logical; stop early when consensus is reached.
#' @return A list with:
#'   \describe{
#'     \item{pheromone_history}{Matrix (waves x K) of pheromone concentrations.}
#'     \item{fitness_history}{Vector of colony fitness per wave.}
#'     \item{decision_history}{Integer vector of colony decisions per wave.}
#'     \item{final_decision}{Index of the chosen site.}
#'   }
#' @export
acar <- function(site_qualities,
                 n_ants    = 50,
                 n_waves   = 100,
                 rho       = 0.1,
                 gamma     = 0.5,
                 alpha     = 1,
                 beta      = 1,
                 noise_sd  = 0.5,
                 early_stop = TRUE) {

  K  <- length(site_qualities)
  sq <- as.numeric(site_qualities)
  tau <- rep(1, K)

  pheromone_history <- matrix(NA_real_, nrow = n_waves, ncol = K)
  fitness_history   <- numeric(n_waves)
  decision_history  <- integer(n_waves)

  actual_waves <- n_waves
  for (t in seq_len(n_waves)) {
    # Ant decision probabilities (Equation 8 in manuscript)
    heuristic <- sq / max(sq)
    prob <- (tau^alpha) * (heuristic^beta)
    prob <- prob / sum(prob)

    # Multinomial ant visits
    visits <- as.integer(rmultinom(1, n_ants, prob))

    # Noisy quality observations
    observed <- sq + rnorm(K, mean = 0, sd = noise_sd)
    observed <- pmax(observed, 0)

    # Pheromone update (Equation 7 in manuscript)
    tau <- (1 - rho) * tau + gamma * visits * observed
    tau <- pmax(tau, 0.01)  # prevent extinction

    pheromone_history[t, ] <- tau
    fitness_history[t]     <- sum(visits * observed) / n_ants
    decision_history[t]    <- which.max(tau)

    # Early stopping: consensus
    if (early_stop && max(tau) > 10 * median(tau[tau > 0])) {
      actual_waves <- t
      break
    }
  }

  list(
    pheromone_history = pheromone_history[seq_len(actual_waves), , drop = FALSE],
    fitness_history   = fitness_history[seq_len(actual_waves)],
    decision_history  = decision_history[seq_len(actual_waves)],
    final_decision    = which.max(tau)
  )
}


#' @title Calculate Quorum Margin
#' @description Computes the normalized pheromone margin between the best and
#'   second-best sites, analogous to boosting margins.
#' @param result Output of \code{acar()}.
#' @return Numeric scalar: the quorum margin.
#' @export
calculate_quorum_margin <- function(result) {
  final_ph  <- result$pheromone_history[nrow(result$pheromone_history), ]
  best      <- result$final_decision
  best_val  <- final_ph[best]
  second    <- max(final_ph[-best])
  (best_val - second) / sum(final_ph)
}
