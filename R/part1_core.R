#' Generate Synthetic Regression Data
#'
#' Creates data with a known sparse signal (first 5 features active)
#' and controllable signal-to-noise ratio. Used to benchmark both
#' random forest and ant colony simulation experiments.
#'
#' @param n Number of observations
#' @param p Number of features
#' @param signal_strength Multiplier for the true signal
#' @param noise_sd Standard deviation of additive Gaussian noise
#' @return A list with components:
#'   \describe{
#'     \item{X}{A data.frame of features (\code{n} rows, \code{p} columns)}
#'     \item{y}{Numeric response vector}
#'     \item{f_true}{True function values (without noise)}
#'   }
#' @export
#' @examples
#' d <- generate_regression_data(n = 200, p = 20, signal_strength = 2, noise_sd = 1)
#' plot(d$f_true, d$y, xlab = "True signal", ylab = "Observed response")
generate_regression_data <- function(n = 1000, p = 50, signal_strength = 1, noise_sd = 1) {
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  f_true <- signal_strength * (X[, 1] + X[, 2] * X[, 3] +
                                 sin(X[, 4]) + X[, 5]^2)
  y <- f_true + stats::rnorm(n, 0, noise_sd)
  list(X = as.data.frame(X), y = y, f_true = f_true)
}

#' Simulate an Ant Colony Decision Process
#'
#' Agent-based simulation of \emph{Temnothorax}-style nest-site selection.
#' Each ant probabilistically explores or follows pheromone trails, makes
#' noisy observations of site quality, and contributes to recruitment.
#' The colony reaches a decision when recruitment exceeds a quorum threshold.
#'
#' @param n_ants Number of ants in the colony
#' @param n_sites Number of candidate nest sites
#' @param site_qualities Numeric vector of true site quality values
#' @param p_explore Probability that an ant explores independently
#'   (vs following pheromone). This is the decorrelation parameter,
#'   analogous to \code{m_try/p} in random forests.
#' @param n_steps Maximum number of simulation time steps
#' @param noise_sd Standard deviation of observation noise
#' @param quorum_threshold Recruitment level that triggers a colony decision
#' @param alpha Learning rate for pheromone updates
#' @param beta Pheromone evaporation rate (0 to 1)
#' @param gamma Recruitment strength multiplier
#' @return A list with components:
#'   \describe{
#'     \item{history}{List of per-step snapshots (pheromone, recruitment, preferences)}
#'     \item{final_preferences}{Colony-average quality estimates per site}
#'     \item{final_recruitment}{Final recruitment levels per site}
#'     \item{decision}{Index of chosen site}
#'     \item{ant_states}{Matrix of visit counts (\code{n_ants} x \code{n_sites})}
#'     \item{ant_preferences}{Matrix of estimated qualities (\code{n_ants} x \code{n_sites})}
#'   }
#' @export
#' @examples
#' sim <- simulate_ant_colony(n_ants = 30, p_explore = 0.4)
#' cat("Colony chose site", sim$decision, "\n")
simulate_ant_colony <- function(
    n_ants           = 50,
    n_sites          = 5,
    site_qualities   = c(10, 8, 6, 4, 2),
    p_explore        = 0.3,
    n_steps          = 100,
    noise_sd         = 2,
    quorum_threshold = 20,
    alpha            = 0.1,
    beta             = 0.05,
    gamma            = 0.2
) {
  ant_states      <- matrix(0, nrow = n_ants, ncol = n_sites)
  ant_preferences <- matrix(0, nrow = n_ants, ncol = n_sites)
  pheromone       <- rep(0, n_sites)
  recruitment     <- rep(0, n_sites)
  decision        <- NA
  history         <- list()

  for (t in seq_len(n_steps)) {
    for (i in seq_len(n_ants)) {
      if (stats::runif(1) < p_explore) {
        site <- sample(n_sites, 1)
      } else {
        if (sum(pheromone) > 0) {
          site <- sample(n_sites, 1, prob = pheromone / sum(pheromone))
        } else {
          site <- sample(n_sites, 1)
        }
      }

      observation <- site_qualities[site] + stats::rnorm(1, 0, noise_sd)
      ant_states[i, site] <- ant_states[i, site] + 1
      n_visits <- ant_states[i, site]
      current_mean <- ant_preferences[i, site]
      ant_preferences[i, site] <- (current_mean * (n_visits - 1) + observation) / n_visits

      if (ant_preferences[i, site] > mean(site_qualities)) {
        recruitment[site] <- recruitment[site] +
          gamma * (ant_preferences[i, site] - mean(site_qualities))
      }
    }

    pheromone <- pheromone * (1 - beta) + recruitment * alpha

    if (max(recruitment) > quorum_threshold) {
      decision <- which.max(recruitment)
      break
    }

    history[[t]] <- list(
      time = t, pheromone = pheromone, recruitment = recruitment,
      ant_preferences = ant_preferences, ant_states = ant_states
    )

    recruitment <- recruitment * 0.8
  }

  if (is.na(decision)) decision <- which.max(recruitment)

  list(
    history           = history,
    final_preferences = colMeans(ant_preferences),
    final_recruitment = recruitment,
    decision          = decision,
    ant_states        = ant_states,
    ant_preferences   = ant_preferences
  )
}

#' Compute Within-Colony Ant Correlation
#'
#' Measures the mean pairwise correlation between ants' preference vectors
#' over all candidate sites within a single colony. This is the correct
#' analogue of tree-tree correlation in a random forest: each ant's
#' preference vector over \eqn{K} sites corresponds to a tree's prediction
#' vector over test points.
#'
#' @param ant_preferences An \code{n_ants x n_sites} matrix of quality
#'   estimates, as returned by \code{\link{simulate_ant_colony}}.
#' @return Mean pairwise correlation (scalar), or \code{NA} if fewer than
#'   2 ants have visited at least 2 sites.
#' @export
#' @examples
#' sim <- simulate_ant_colony(n_ants = 30, p_explore = 0.3)
#' within_colony_correlation(sim$ant_preferences)
within_colony_correlation <- function(ant_preferences) {
  visited <- rowSums(ant_preferences != 0) >= 2
  if (sum(visited) < 2) return(NA_real_)
  prefs <- ant_preferences[visited, , drop = FALSE]
  cm <- stats::cor(t(prefs), use = "pairwise.complete.obs")
  mean(cm[lower.tri(cm)], na.rm = TRUE)
}

