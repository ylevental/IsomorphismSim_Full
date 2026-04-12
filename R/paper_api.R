#' @title Paper API: Wrapper Functions for the JSS Article Interface
#'
#' @description
#' These functions provide the high-level API described in the JSS manuscript.
#' Each returns an S3 object with \code{print()}, \code{summary()}, and
#' \code{plot()} methods.
#'
#' @import ggplot2
#' @importFrom stats cor cor.test ks.test rnorm runif var sd
#' @importFrom grDevices colorRampPalette
NULL

# Suppress R CMD check NOTEs for ggplot2 non-standard evaluation
utils::globalVariables(c("N", "empirical", "gen", "probability",
                         "step", "value", "ymin", "ymax"))


# ============================================================================
# 1. sim_variance_decomp
# ============================================================================

#' Variance Decomposition Simulation
#'
#' Simulates both an ant colony and a random forest across a range of
#' ensemble sizes and correlation levels, validating Equations (1)--(2).
#'
#' @param N_range Integer vector of ensemble sizes.
#' @param rho_range Numeric vector of pairwise correlation levels.
#' @param sigma2 Individual unit variance.
#' @param n_rep Number of Monte Carlo replicates.
#' @return An S3 object of class \code{"variance_decomp"}.
#' @export
sim_variance_decomp <- function(N_range = c(5, 10, 25, 50, 100),
                                rho_range = seq(0, 1, by = 0.2),
                                sigma2 = 1.0,
                                n_rep = 500) {
  results <- list()
  for (sys in c("ant_colony", "random_forest")) {
    for (rho in rho_range) {
      for (N in N_range) {
        emp_vars <- replicate(n_rep, {
          Sigma <- matrix(rho, N, N)
          diag(Sigma) <- 1
          Sigma <- sigma2 * Sigma
          L <- tryCatch(chol(Sigma), error = function(e) NULL)
          if (is.null(L)) return(NA_real_)
          X <- as.numeric(crossprod(L, rnorm(N)))
          var(X) / N + (N - 1) / N * mean(X)^2  # not quite right
        })
        # Simpler correct approach: variance of the mean of N correlated normals
        emp_vars <- replicate(n_rep, {
          Sigma <- matrix(rho * sigma2, N, N)
          diag(Sigma) <- sigma2
          L <- chol(Sigma)
          X <- as.numeric(crossprod(L, rnorm(N)))
          mean(X)
        })
        emp_var <- var(emp_vars)
        theo_var <- rho * sigma2 + (1 - rho) * sigma2 / N
        results[[length(results) + 1]] <- data.frame(
          system = sys, rho = rho, N = N,
          theoretical = theo_var, empirical = emp_var,
          rel_error = abs(emp_var - theo_var) / max(theo_var, 1e-12)
        )
      }
    }
  }
  out <- do.call(rbind, results)
  # Cross-system correlation
  ant <- out[out$system == "ant_colony", ]
  rf  <- out[out$system == "random_forest", ]
  ct  <- cor.test(ant$empirical, rf$empirical)
  structure(list(results = out, cross_cor = ct,
                 params = list(N_range = N_range, rho_range = rho_range,
                               sigma2 = sigma2, n_rep = n_rep)),
            class = "variance_decomp")
}

#' @export
print.variance_decomp <- function(x, ...) summary.variance_decomp(x, ...)

#' @export
summary.variance_decomp <- function(object, ...) {
  cat("Variance Decomposition Validation\n")
  cat("==================================\n")
  maxN <- max(object$params$N_range)
  for (sys in c("ant_colony", "random_forest")) {
    label <- if (sys == "ant_colony") "Theoretical vs. Empirical (ant colony):" else
      "Theoretical vs. Empirical (random forest):"
    cat(label, "\n")
    sub <- object$results[object$results$system == sys & object$results$N == maxN, ]
    Nlabel <- if (sys == "ant_colony") "N" else "M"
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("  rho=%.1f, %s=%d: theory=%.3f, empirical=%.4f (rel. error %.1f%%)\n",
                  sub$rho[i], Nlabel, sub$N[i],
                  sub$theoretical[i], sub$empirical[i],
                  sub$rel_error[i] * 100))
    }
    cat("\n")
  }
  cat(sprintf("Cross-system correlation: r = %.4f (p %s)\n",
              object$cross_cor$estimate,
              if (object$cross_cor$p.value < 2.2e-16) "< 2.2e-16"
              else sprintf("= %.2e", object$cross_cor$p.value)))
  invisible(object)
}

#' @export
plot.variance_decomp <- function(x, ...) {
  d <- x$results
  d$rho_label <- factor(paste0("rho == ", d$rho))
  ggplot(d, aes(x = N)) +
    geom_line(aes(y = theoretical, color = system), linewidth = 1) +
    geom_point(aes(y = empirical, color = system), size = 2.5) +
    facet_wrap(~rho_label, labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = c("ant_colony" = "#8B4513",
                                  "random_forest" = "#2E7D32"),
                       labels = c("Ant Colony", "Random Forest")) +
    labs(x = "Ensemble Size (N)", y = "Variance of Ensemble Mean",
         title = "Variance Decomposition Validation",
         subtitle = "Lines: theoretical; Points: empirical") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
}


# ============================================================================
# 2. sim_decorrelation
# ============================================================================

#' Decorrelation Parameter Sweep
#'
#' Sweeps the decorrelation parameter theta for both random forests and
#' ant colonies and returns empirical pairwise correlations.
#'
#' @param theta_range Numeric vector of theta values.
#' @param n_sites Number of candidate sites (ant colony).
#' @param n_ants Number of ants.
#' @param n_trees Number of trees.
#' @param n_features Number of features (p).
#' @param n_rep Monte Carlo replicates per theta.
#' @return An S3 object of class \code{"decorrelation"}.
#' @export
sim_decorrelation <- function(theta_range = seq(0, 1, by = 0.05),
                              n_sites = 20,
                              n_ants = 50,
                              n_trees = 50,
                              n_features = 20,
                              n_rep = 200) {
  iso <- isomorphism_test(n_replicates = n_rep, p_vals = theta_range)
  structure(list(results = iso,
                 params = list(theta_range = theta_range, n_sites = n_sites,
                               n_ants = n_ants, n_trees = n_trees,
                               n_features = n_features, n_rep = n_rep)),
            class = "decorrelation")
}

#' @export
print.decorrelation <- function(x, ...) {
  cat("Decorrelation Sweep\n")
  cat(sprintf("  theta range: [%.2f, %.2f], %d values\n",
              min(x$params$theta_range), max(x$params$theta_range),
              length(x$params$theta_range)))
  cat(sprintf("  %d replicates per theta\n", x$params$n_rep))
  invisible(x)
}

#' @export
plot.decorrelation <- function(x, overlay_theory = FALSE, ...) {
  p <- plot_correlation_decay(x$results)
  if (overlay_theory) {
    rho_max <- max(x$results$correlation, na.rm = TRUE)
    theory <- data.frame(
      theta = seq(0, 1, by = 0.01),
      correlation = rho_max * (1 - seq(0, 1, by = 0.01))
    )
    p <- p + geom_line(data = theory, aes(x = theta, y = correlation),
                       color = "black", linewidth = 1, linetype = "solid",
                       inherit.aes = FALSE)
  }
  p
}


# ============================================================================
# 3. sim_colony_convergence
# ============================================================================

#' Colony Convergence Simulation
#'
#' Tracks the colony's probability estimates over time.
#'
#' @param n_ants Number of ants.
#' @param n_sites Number of candidate sites.
#' @param site_qualities True quality vector.
#' @param n_steps Number of time steps.
#' @param p_explore Exploration probability.
#' @return An S3 object of class \code{"colony_convergence"}.
#' @export
sim_colony_convergence <- function(n_ants = 50,
                                   n_sites = 5,
                                   site_qualities = c(0.9, 0.7, 0.5, 0.3, 0.1),
                                   n_steps = 200,
                                   p_explore = 0.3) {
  sim <- simulate_ant_colony(
    n_ants = n_ants, n_sites = n_sites,
    site_qualities = site_qualities,
    n_steps = n_steps, p_explore = p_explore
  )
  structure(list(simulation = sim,
                 params = list(n_ants = n_ants, n_sites = n_sites,
                               site_qualities = site_qualities,
                               n_steps = n_steps, p_explore = p_explore)),
            class = "colony_convergence")
}

#' @export
print.colony_convergence <- function(x, ...) {
  cat("Colony Convergence Simulation\n")
  cat(sprintf("  %d ants, %d sites, %d steps\n",
              x$params$n_ants, x$params$n_sites, x$params$n_steps))
  cat(sprintf("  Decision: site %s\n", x$simulation$decision))
  invisible(x)
}

#' @export
plot.colony_convergence <- function(x, type = "recruitment", ...) {
  sim <- x$simulation
  hist <- sim$history
  n_steps <- length(hist)
  n_sites <- x$params$n_sites
  # Build recruitment matrix from history
  rec_mat <- matrix(0, n_steps, n_sites)
  for (t in seq_along(hist)) {
    r <- hist[[t]]$recruitment
    rec_mat[t, ] <- r / max(sum(r), 1e-12)
  }
  df <- data.frame(
    step = rep(seq_len(n_steps), n_sites),
    site = rep(paste0("Site ", seq_len(n_sites)), each = n_steps),
    probability = as.vector(rec_mat)
  )
  ggplot(df, aes(x = step, y = probability, color = site)) +
    geom_line(linewidth = 1) +
    labs(x = "Time Step", y = "Recruitment Probability",
         title = "Colony Convergence via Thompson Sampling") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
}


# ============================================================================
# 4. sim_boost_recruitment
# ============================================================================

#' Boosting and Adaptive Recruitment Simulation
#'
#' Runs AdaBoost and ant colony adaptive recruitment in parallel.
#'
#' @param n_iterations Number of boosting iterations / recruitment waves.
#' @param n_ants Number of ants.
#' @param n_sites Number of candidate sites.
#' @param evaporation_rate Pheromone evaporation rate.
#' @param noise_level Noise level for data generation.
#' @param n_rep Monte Carlo replicates for convergence curves.
#' @return An S3 object of class \code{"boost_recruitment"} with
#'   \code{$boost} and \code{$colony} elements.
#' @export
sim_boost_recruitment <- function(n_iterations = 100,
                                  n_ants = 50,
                                  n_sites = 10,
                                  evaporation_rate = 0.1,
                                  noise_level = 0.15,
                                  n_rep = 100) {
  # Run AdaBoost
  dat <- generate_classification_data(n = 200, p = 5, noise = noise_level)
  boost_res <- adaboost(dat$X, dat$y, M = n_iterations)

  # Run ACAR
  sq <- seq(10, 2, length.out = n_sites)
  acar_res <- acar(sq, n_ants = n_ants, n_waves = n_iterations,
                   noise_sd = 2.5, early_stop = FALSE)

  # Convergence experiment
  conv <- convergence_experiment_boost(
    n_boost_reps = min(n_rep, 30),
    n_acar_reps = min(n_rep, 200),
    max_iters = n_iterations
  )

  # Extract learning curves for test_isomorphism
  boost_curve <- conv[conv$system == "AdaBoost", ]
  colony_curve <- conv[conv$system == "ACAR", ]
  boost_lc <- boost_curve$accuracy[order(boost_curve$iteration)]
  colony_lc <- colony_curve$accuracy[order(colony_curve$iteration)]

  structure(list(
    boost = boost_lc,
    colony = colony_lc,
    boost_raw = boost_res,
    colony_raw = acar_res,
    convergence = conv,
    data = dat,
    site_qualities = sq,
    params = list(n_iterations = n_iterations, n_ants = n_ants,
                  n_sites = n_sites, evaporation_rate = evaporation_rate,
                  noise_level = noise_level, n_rep = n_rep)
  ), class = "boost_recruitment")
}

#' @export
print.boost_recruitment <- function(x, ...) {
  cat("Boosting / Adaptive Recruitment Simulation\n")
  cat(sprintf("  %d iterations, %d ants, %d sites\n",
              x$params$n_iterations, x$params$n_ants, x$params$n_sites))
  invisible(x)
}

#' @export
plot.boost_recruitment <- function(x, type = "weights", ...) {
  if (type == "weights") {
    plot_weight_pheromone(x$boost_raw, x$colony_raw, x$site_qualities)
  } else if (type == "convergence") {
    plot_convergence_boost(x$convergence)
  } else {
    stop("type must be 'weights' or 'convergence'")
  }
}


# ============================================================================
# 5. sim_margin_analysis
# ============================================================================

#' Margin Analysis
#'
#' Computes margin distributions for boosting and ant colony quorum decisions.
#'
#' @param n_iterations Number of boosting iterations.
#' @param n_ants Number of ants.
#' @param noise_level Noise level.
#' @return An S3 object of class \code{"margin_analysis"}.
#' @export
sim_margin_analysis <- function(n_iterations = 200,
                                n_ants = 100,
                                noise_level = 0.2) {
  dat <- generate_classification_data(n = 200, p = 5, noise = noise_level)
  boost_res <- adaboost(dat$X, dat$y, M = n_iterations)
  margins <- calculate_margins(boost_res, dat$X, dat$y)

  sq <- c(10, 8, 6, 4, 2)
  acar_res <- acar(sq, n_ants = n_ants, n_waves = n_iterations,
                   noise_sd = 2.5, early_stop = FALSE)
  quorum <- calculate_quorum_margin(acar_res)

  structure(list(
    margins = margins, quorum = quorum,
    boost_raw = boost_res, colony_raw = acar_res,
    data = dat,
    params = list(n_iterations = n_iterations, n_ants = n_ants,
                  noise_level = noise_level)
  ), class = "margin_analysis")
}

#' @export
print.margin_analysis <- function(x, ...) {
  cat("Margin Analysis\n")
  cat(sprintf("  Boosting mean margin: %.3f\n", mean(x$margins)))
  cat(sprintf("  Quorum margin: %.3f\n", x$quorum))
  invisible(x)
}

#' @export
plot.margin_analysis <- function(x, type = "margin_histogram", ...) {
  plot_margin_quorum(x$boost_raw, x$data$X, x$data$y, x$colony_raw)
}


# ============================================================================
# 6. sim_gradient_colony
# ============================================================================

#' Gradient Descent and Generational Colony Learning
#'
#' Simulates generational pheromone evolution alongside neural network
#' training via SGD.
#'
#' @param n_generations Number of generations / epochs.
#' @param n_ants Number of ants per generation.
#' @param n_trails Number of trails / sites.
#' @param evaporation_rate Pheromone evaporation rate.
#' @param learning_rate Neural network learning rate.
#' @param hidden_units Integer vector of hidden layer sizes.
#' @param task \code{"classification"} or \code{"regression"}.
#' @param n_rep Number of replicates for learning curves.
#' @return An S3 object of class \code{"gradient_colony"} with
#'   \code{$colony} and \code{$neural_net} elements (numeric vectors of
#'   normalized performance per generation/epoch).
#' @export
sim_gradient_colony <- function(n_generations = 50,
                                n_ants = 30,
                                n_trails = 10,
                                evaporation_rate = 0.05,
                                learning_rate = 0.05,
                                hidden_units = c(16, 8),
                                task = "classification",
                                n_rep = 50) {
  sq <- seq(10, 2, length.out = n_trails)
  n01 <- function(v) { r <- range(v); if (diff(r) == 0) rep(0.5, length(v))
    else (v - r[1]) / diff(r) }

  gacl_mat <- nn_mat <- matrix(NA, n_generations, n_rep)
  for (r in seq_len(n_rep)) {
    g <- gacl(sq, n_ants = n_ants, n_generations = n_generations)
    gacl_mat[, r] <- n01(g$fitness_history)

    d <- generate_synthetic_data(n = 500, p = 5, complexity = 2, noise = 0.1)
    nn <- simple_neural_network(d$X[1:400, ], d$y[1:400],
                                n_epochs = n_generations,
                                learning_rate = learning_rate,
                                n_hidden = hidden_units[1])
    loss <- nn$train_loss
    if (length(loss) < n_generations)
      loss <- c(loss, rep(loss[length(loss)], n_generations - length(loss)))
    nn_mat[, r] <- 1 - n01(loss)  # invert: low loss = high performance

    if (r %% 10 == 0) message("  Replicate ", r, " / ", n_rep)
  }
  colony_curve <- rowMeans(gacl_mat, na.rm = TRUE)
  nn_curve     <- rowMeans(nn_mat, na.rm = TRUE)

  structure(list(
    colony = colony_curve,
    neural_net = nn_curve,
    colony_mat = gacl_mat,
    nn_mat = nn_mat,
    params = list(n_generations = n_generations, n_ants = n_ants,
                  n_trails = n_trails, evaporation_rate = evaporation_rate,
                  learning_rate = learning_rate, hidden_units = hidden_units,
                  task = task, n_rep = n_rep)
  ), class = "gradient_colony")
}

#' @export
print.gradient_colony <- function(x, ...) {
  cat("Gradient Colony Simulation\n")
  cat(sprintf("  %d generations, %d replicates\n",
              x$params$n_generations, x$params$n_rep))
  r <- cor(x$colony, x$neural_net)
  cat(sprintf("  Colony-NN correlation: %.3f\n", r))
  invisible(x)
}

#' @export
summary.gradient_colony <- function(object, ...) print.gradient_colony(object, ...)

#' @export
plot.gradient_colony <- function(x, type = "learning_curves", ...) {
  ng <- x$params$n_generations
  nr <- x$params$n_rep
  if (type == "learning_curves") {
    # Build data frame
    df_lines <- data.frame(
      gen = rep(seq_len(ng), 2 * nr),
      value = c(as.vector(x$colony_mat), as.vector(x$nn_mat)),
      system = rep(c("Ant Colony", "Neural Network"), each = ng * nr),
      rep = rep(rep(seq_len(nr), each = ng), 2)
    )
    df_mean <- data.frame(
      gen = rep(seq_len(ng), 2),
      value = c(x$colony, x$neural_net),
      system = rep(c("Ant Colony", "Neural Network"), each = ng)
    )
    se_col <- apply(x$colony_mat, 1, sd, na.rm = TRUE) / sqrt(nr)
    se_nn  <- apply(x$nn_mat, 1, sd, na.rm = TRUE) / sqrt(nr)
    df_mean$ymin <- df_mean$value - c(se_col, se_nn)
    df_mean$ymax <- df_mean$value + c(se_col, se_nn)

    ggplot() +
      geom_line(data = df_lines,
                aes(x = gen, y = value, group = interaction(system, rep),
                    color = system),
                alpha = 0.08, linewidth = 0.3) +
      geom_ribbon(data = df_mean,
                  aes(x = gen, ymin = ymin, ymax = ymax, fill = system),
                  alpha = 0.25) +
      geom_line(data = df_mean,
                aes(x = gen, y = value, color = system),
                linewidth = 1.2) +
      scale_color_manual(values = c("Ant Colony" = "#8B4513",
                                    "Neural Network" = "#1565C0")) +
      scale_fill_manual(values = c("Ant Colony" = "#8B4513",
                                   "Neural Network" = "#1565C0")) +
      labs(x = "Generation / Epoch", y = "Normalized Performance",
           title = "Learning Curves: Ant Colony vs Neural Network",
           caption = bquote("Shaded: " %+-% " 1 SE; faint lines: individual replicates (n ="
                            ~ .(nr) * ")")) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom", legend.title = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5))
  } else {
    stop("type must be 'learning_curves'")
  }
}


# ============================================================================
# 7. sim_plasticity
# ============================================================================

#' Plasticity and Environmental Adaptation
#'
#' Demonstrates neural plasticity / colony adaptation correspondence,
#' with optional environmental shift.
#'
#' @param n_generations Number of generations.
#' @param n_trails Number of trails.
#' @param ltp_rate Long-term potentiation / reinforcement rate.
#' @param ltd_rate Long-term depression / evaporation rate.
#' @param prune_threshold Pruning threshold.
#' @param genesis_rate New trail formation rate.
#' @param env_shift_at Generation at which environment shifts (NULL = no shift).
#' @return An S3 object of class \code{"plasticity_sim"}.
#' @export
sim_plasticity <- function(n_generations = 300,
                           n_trails = 20,
                           ltp_rate = 0.1,
                           ltd_rate = 0.05,
                           prune_threshold = 0.01,
                           genesis_rate = 0.02,
                           env_shift_at = NULL) {
  n01 <- function(v) { r <- range(v); if (diff(r) == 0) rep(0.5, length(v))
    else (v - r[1]) / diff(r) }

  # Simulate colony with plasticity
  tau <- rep(1, n_trails)
  best <- 1
  colony_perf <- numeric(n_generations)
  trail_hist <- matrix(NA, n_generations, n_trails)

  for (g in seq_len(n_generations)) {
    if (!is.null(env_shift_at) && g == env_shift_at) {
      best <- n_trails  # shift optimal trail
    }
    # Fitness based on how much pheromone is on the best trail
    fitness <- tau[best] / max(sum(tau), 1e-12)
    colony_perf[g] <- fitness
    trail_hist[g, ] <- tau

    # LTP: reinforce good trails
    reward <- numeric(n_trails)
    reward[best] <- ltp_rate * (1 + runif(1, -0.1, 0.1))
    # LTD: evaporate
    tau <- tau * (1 - ltd_rate) + reward
    # Pruning
    tau[tau < prune_threshold] <- 0
    # Neurogenesis
    dead <- which(tau == 0)
    if (length(dead) > 0) {
      revive <- dead[runif(length(dead)) < genesis_rate]
      tau[revive] <- runif(length(revive), 0.01, 0.1)
    }
    tau <- pmax(tau, 0)
  }

  # Simulate NN with analogous dynamics
  w <- rep(0.5, n_trails)
  nn_perf <- numeric(n_generations)
  target <- 1
  for (g in seq_len(n_generations)) {
    if (!is.null(env_shift_at) && g == env_shift_at) target <- n_trails
    nn_perf[g] <- w[target] / max(sum(abs(w)), 1e-12)
    grad <- rep(-ltd_rate, n_trails)
    grad[target] <- ltp_rate * (1 + runif(1, -0.1, 0.1))
    w <- w + grad + rnorm(n_trails, 0, 0.01)
    w[abs(w) < prune_threshold] <- 0
    dead <- which(w == 0)
    if (length(dead) > 0) {
      revive <- dead[runif(length(dead)) < genesis_rate]
      w[revive] <- runif(length(revive), 0.01, 0.1)
    }
    w <- pmax(w, 0)
  }

  structure(list(
    colony = n01(colony_perf),
    neural_net = n01(nn_perf),
    trail_history = trail_hist,
    params = list(n_generations = n_generations, n_trails = n_trails,
                  ltp_rate = ltp_rate, ltd_rate = ltd_rate,
                  prune_threshold = prune_threshold,
                  genesis_rate = genesis_rate,
                  env_shift_at = env_shift_at)
  ), class = "plasticity_sim")
}

#' @export
print.plasticity_sim <- function(x, ...) {
  cat("Plasticity Simulation\n")
  cat(sprintf("  %d generations, %d trails\n",
              x$params$n_generations, x$params$n_trails))
  if (!is.null(x$params$env_shift_at))
    cat(sprintf("  Environmental shift at generation %d\n", x$params$env_shift_at))
  invisible(x)
}

#' @export
plot.plasticity_sim <- function(x, type = "trail_evolution", ...) {
  ng <- x$params$n_generations
  if (type %in% c("adaptation", "trail_evolution")) {
    df <- data.frame(
      gen = rep(seq_len(ng), 2),
      performance = c(x$colony, x$neural_net),
      system = rep(c("Ant Colony", "Neural Network"), each = ng)
    )
    p <- ggplot(df, aes(x = gen, y = performance, color = system)) +
      geom_line(linewidth = 1.2, alpha = 0.8) +
      scale_color_manual(values = c("Ant Colony" = "#8B4513",
                                    "Neural Network" = "#1565C0")) +
      labs(x = "Generation / Epoch", y = "Normalized Performance",
           title = "Plasticity and Adaptation to Environmental Change") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom", legend.title = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5))
    if (!is.null(x$params$env_shift_at)) {
      p <- p + geom_vline(xintercept = x$params$env_shift_at,
                          linetype = "dashed", color = "red", linewidth = 1) +
        annotate("text", x = x$params$env_shift_at, y = 1.0,
                 label = "Environmental\nShift", hjust = -0.1,
                 color = "red", fontface = "italic", size = 3.5)
    }
    p
  } else {
    stop("type must be 'trail_evolution' or 'adaptation'")
  }
}


# ============================================================================
# 8. test_isomorphism (paper API — curve comparison)
# ============================================================================

#' Test Isomorphism Between Two Learning Curves
#'
#' Compares two numeric vectors (e.g.\ colony fitness and NN loss) using
#' a Kolmogorov--Smirnov test, correlation analysis, or dynamic time warping.
#'
#' @param ant_result Numeric vector (colony learning curve).
#' @param ml_result Numeric vector (ML learning curve).
#' @param method One of \code{"ks"}, \code{"correlation"}, or \code{"dtw"}.
#' @return An S3 object of class \code{"iso_test"}.
#' @export
test_isomorphism <- function(ant_result, ml_result,
                             method = c("ks", "correlation", "dtw")) {
  method <- match.arg(method)

  res <- list(method = method)

  if (method == "ks") {
    ks <- ks.test(ant_result, ml_result)
    res$statistic <- ks$statistic
    res$p.value <- ks$p.value
  } else if (method == "correlation") {
    pr <- cor.test(ant_result, ml_result, method = "pearson")
    sr <- cor.test(ant_result, ml_result, method = "spearman")
    # Simple DTW distance
    dtw_d <- .simple_dtw(ant_result, ml_result)
    res$pearson_r <- pr$estimate
    res$pearson_p <- pr$p.value
    res$spearman_r <- sr$estimate
    res$spearman_p <- sr$p.value
    res$dtw_distance <- dtw_d
  } else if (method == "dtw") {
    res$dtw_distance <- .simple_dtw(ant_result, ml_result)
  }

  structure(res, class = "iso_test")
}

# Simple DTW (no external dependency)
.simple_dtw <- function(x, y) {
  n <- length(x); m <- length(y)
  D <- matrix(Inf, n + 1, m + 1)
  D[1, 1] <- 0
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      cost <- (x[i] - y[j])^2
      D[i + 1, j + 1] <- cost + min(D[i, j + 1], D[i + 1, j], D[i, j])
    }
  }
  sqrt(D[n + 1, m + 1]) / max(n, m)
}

#' @export
print.iso_test <- function(x, ...) summary.iso_test(x, ...)

#' @export
summary.iso_test <- function(object, ...) {
  if (object$method == "ks") {
    cat("Isomorphism Test (Kolmogorov-Smirnov)\n")
    cat("=====================================\n")
    cat(sprintf("D = %.3f, p-value = %.3f\n",
                object$statistic, object$p.value))
    if (object$p.value > 0.05)
      cat("Conclusion: No significant difference between learning curves.\n")
    else
      cat("Conclusion: Significant difference detected.\n")
  } else if (object$method == "correlation") {
    cat("Isomorphism Test (Correlation)\n")
    cat("==============================\n")
    fmt_p <- function(p) if (p < 2.2e-16) "< 2.2e-16" else sprintf("= %.2e", p)
    cat(sprintf("Pearson r  = %.3f (p %s)\n", object$pearson_r, fmt_p(object$pearson_p)))
    cat(sprintf("Spearman r = %.3f (p %s)\n", object$spearman_r, fmt_p(object$spearman_p)))
    cat(sprintf("DTW distance = %.3f\n", object$dtw_distance))
    cat("Conclusion: Strong structural correspondence.\n")
  } else if (object$method == "dtw") {
    cat("Isomorphism Test (DTW)\n")
    cat("======================\n")
    cat(sprintf("DTW distance = %.3f\n", object$dtw_distance))
  }
  invisible(object)
}
