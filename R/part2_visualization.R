#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_histogram geom_col
#'   geom_vline geom_point geom_hline geom_smooth stat_summary scale_x_log10
#'   scale_color_manual scale_fill_manual labs theme_minimal theme element_text
#'   element_blank annotate ggsave facet_wrap
#' @import patchwork
#' @importFrom dplyr group_by summarise mutate filter n %>%
#' @importFrom tidyr pivot_longer
NULL

# Colour palette ---------------------------------------------------------
.BLUE   <- "#4472C4"
.RED    <- "#C44E52"
.GREEN  <- "#55A868"
.BROWN  <- "#8C564B"
.ORANGE <- "#DD8452"
.PURPLE <- "#937DC2"

# Common theme -----------------------------------------------------------
.theme_pub <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}


# ========================================================================
# Figure 1: Weak Learnability
# ========================================================================
#' @title Plot Figure 1: Weak Learnability Theorem
#' @param results Output of \code{weak_learnability_experiment()}.
#' @return A \code{ggplot} object.
#' @export
plot_weak_learnability <- function(results) {
  summ <- results %>%
    group_by(gamma, waves) %>%
    summarise(
      accuracy = mean(accuracy),
      se       = sd(accuracy) / sqrt(n()),
      .groups  = "drop"
    )

  # Theoretical curves
  theory <- expand.grid(
    gamma = unique(summ$gamma),
    waves = seq(1, max(summ$waves), length.out = 200)
  )
  theory$theoretical <- 1 - exp(-theory$gamma^2 * theory$waves / 2)

  ggplot(summ, aes(x = waves, y = accuracy,
                   color = factor(gamma), fill = factor(gamma))) +
    geom_ribbon(aes(ymin = accuracy - se, ymax = accuracy + se),
                alpha = 0.15, color = NA) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_line(data = theory,
              aes(x = waves, y = theoretical, color = factor(gamma)),
              linetype = "dashed", alpha = 0.5, size = 0.8,
              inherit.aes = FALSE) +
    scale_x_log10() +
    scale_color_manual(
      values = c(.BLUE, .RED, .GREEN, .ORANGE),
      labels = function(x) parse(text = paste0("gamma == ", x))
    ) +
    scale_fill_manual(
      values = c(.BLUE, .RED, .GREEN, .ORANGE),
      labels = function(x) parse(text = paste0("gamma == ", x))
    ) +
    labs(
      title = "The Strength of Weak Learnability in Ant Colonies",
      x     = "Number of Recruitment Waves (log scale)",
      y     = "Probability of Correct Decision",
      color = NULL, fill = NULL
    ) +
    annotate("text", x = 3, y = 0.95, label = "Dashed: theoretical\n1 - exp(-gamma^2 T / 2)",
             fontface = "italic", size = 3.5, hjust = 0) +
    coord_cartesian(ylim = c(0.4, 1.05)) +
    .theme_pub()
}


# ========================================================================
# Figure 2: Weight vs Pheromone Evolution
# ========================================================================
#' @title Plot Figure 2: Weight vs Pheromone Evolution
#' @param boost_result Output of \code{adaboost()}.
#' @param acar_result Output of \code{acar()}.
#' @param site_qualities Numeric vector used in the ACAR call.
#' @return A combined \code{patchwork} plot.
#' @export
plot_weight_pheromone <- function(boost_result, acar_result, site_qualities) {
  M  <- nrow(boost_result$weight_history)
  wh <- boost_result$weight_history

  # Panel (a): weight trajectories
  # Sample ~80 instances for visibility
  n_inst  <- ncol(wh)
  samp_ix <- sort(sample(seq_len(n_inst), min(80, n_inst)))
  weight_df <- data.frame(
    iteration = rep(seq_len(M), times = length(samp_ix)),
    weight    = as.vector(wh[, samp_ix]),
    instance  = rep(samp_ix, each = M)
  )
  med_df <- data.frame(
    iteration = seq_len(M),
    median_w  = apply(wh, 1, median)
  )

  p1 <- ggplot() +
    geom_line(data = weight_df,
              aes(x = iteration, y = weight, group = instance),
              alpha = 0.12, color = .BLUE, size = 0.4) +
    geom_line(data = med_df,
              aes(x = iteration, y = median_w),
              color = "darkblue", size = 1.5) +
    labs(title = "(a) AdaBoost: Instance Weight Evolution",
         x = "Iteration", y = expression(D[t](i))) +
    .theme_pub()

  # Panel (b): pheromone trajectories
  ph <- acar_result$pheromone_history
  K  <- ncol(ph)
  pheromone_df <- data.frame(
    wave      = rep(seq_len(nrow(ph)), times = K),
    pheromone = as.vector(ph),
    site      = rep(paste0("Site ", seq_len(K), " (Q=", site_qualities, ")"),
                    each = nrow(ph))
  )

  site_cols <- c(.RED, .GREEN, .BLUE, .ORANGE, .PURPLE)[seq_len(K)]
  names(site_cols) <- unique(pheromone_df$site)

  p2 <- ggplot(pheromone_df, aes(x = wave, y = pheromone, color = site)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = site_cols) +
    labs(title = expression("(b) ACAR: Pheromone " * tau[j](t) * " Evolution"),
         x = "Recruitment Wave",
         y = expression(tau[j](t)),
         color = NULL) +
    .theme_pub()

  (p1 | p2) +
    plot_annotation(
      title = expression(bold("Isomorphic Evolution: Instance Weights") ~
                           symbol("\253") ~ bold("Pheromone Concentrations")),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 15))
    )
}


# ========================================================================
# Figure 3: Margin vs Quorum
# ========================================================================
#' @title Plot Figure 3: Margin vs Quorum
#' @param boost_result Output of \code{adaboost()}.
#' @param X Predictor matrix used in boosting.
#' @param y Label vector used in boosting.
#' @param acar_result Output of \code{acar()}.
#' @return A combined \code{patchwork} plot.
#' @export
plot_margin_quorum <- function(boost_result, X, y, acar_result) {
  margins <- calculate_margins(boost_result, X, y)
  qm      <- calculate_quorum_margin(acar_result)

  p1 <- ggplot(data.frame(margin = margins), aes(x = margin)) +
    geom_histogram(bins = 30, fill = .BLUE, alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = .RED, size = 1.2) +
    geom_vline(xintercept = mean(margins), color = "darkblue", size = 1.2) +
    annotate("text", x = mean(margins) + 0.05, y = Inf, vjust = 2,
             label = sprintf("mean = %.3f", mean(margins)), size = 3.8) +
    labs(title = "(a) Boosting: Margin Distribution",
         x = expression("Margin " * rho[i]), y = "Count") +
    .theme_pub()

  # Quorum bar plot
  final_ph <- acar_result$pheromone_history[nrow(acar_result$pheromone_history), ]
  best     <- acar_result$final_decision
  bar_df   <- data.frame(
    site      = seq_along(final_ph),
    pheromone = final_ph,
    is_best   = seq_along(final_ph) == best
  )

  p2 <- ggplot(bar_df, aes(x = factor(site), y = pheromone, fill = is_best)) +
    geom_col(alpha = 0.8, color = "white") +
    scale_fill_manual(values = c("TRUE" = .RED, "FALSE" = .BLUE), guide = "none") +
    annotate("text", x = 3, y = max(final_ph) * 0.6,
             label = sprintf("Quorum margin\nmu = %.3f", qm),
             fontface = "bold", size = 4.5) +
    labs(title = sprintf("(b) Ant Colony: Quorum Margin = %.3f", qm),
         x = "Site", y = expression(tau[j](T))) +
    .theme_pub()

  (p1 | p2) +
    plot_annotation(
      title = expression(bold("Isomorphic Margin Concepts: Boosting Margin") ~
                           symbol("\253") ~ bold("Quorum Margin")),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 15))
    )
}


# ========================================================================
# Figure 4: Convergence Rates
# ========================================================================
#' @title Plot Figure 4: Convergence Rates
#' @param results Output of \code{convergence_experiment_boost()}.
#' @return A \code{ggplot} object.
#' @export
plot_convergence_boost <- function(results) {
  ggplot(results, aes(x = iteration, y = accuracy, color = system,
                      fill = system)) +
    stat_summary(fun = mean, geom = "line", size = 1.5) +
    stat_summary(fun.data = function(y) {
      data.frame(ymin = mean(y) - sd(y) / sqrt(length(y)),
                 ymax = mean(y) + sd(y) / sqrt(length(y)))
    }, geom = "ribbon", alpha = 0.2, color = NA) +
    scale_color_manual(values = c("AdaBoost" = .BLUE, "ACAR" = .BROWN)) +
    scale_fill_manual(values  = c("AdaBoost" = .BLUE, "ACAR" = .BROWN)) +
    labs(
      title = "Convergence Rates: AdaBoost vs ACAR",
      x     = "Iteration / Wave",
      y     = "Accuracy / Probability of Correct Decision",
      color = NULL, fill = NULL
    ) +
    coord_cartesian(ylim = c(0.4, 1.05)) +
    .theme_pub()
}


# ========================================================================
# Figure 5: Noise Robustness
# ========================================================================
#' @title Plot Figure 5: Noise Robustness
#' @param results Output of \code{noise_experiment_boost()}.
#' @return A \code{ggplot} object.
#' @export
plot_noise_robustness_boost <- function(results) {
  ggplot(results, aes(x = noise, y = accuracy, color = system,
                      fill = system)) +
    stat_summary(fun = mean, geom = "line", size = 1.5) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = function(y) {
      data.frame(ymin = mean(y) - sd(y) / sqrt(length(y)),
                 ymax = mean(y) + sd(y) / sqrt(length(y)))
    }, geom = "ribbon", alpha = 0.2, color = NA) +
    scale_color_manual(values = c("AdaBoost" = .BLUE, "ACAR" = .BROWN)) +
    scale_fill_manual(values  = c("AdaBoost" = .BLUE, "ACAR" = .BROWN)) +
    labs(
      title = "Noise Robustness: AdaBoost vs ACAR",
      x     = "Noise Level",
      y     = "Accuracy / Probability of Correct Decision",
      color = NULL, fill = NULL
    ) +
    coord_cartesian(ylim = c(0.3, 1.05)) +
    .theme_pub()
}


# ========================================================================
# Master function: generate all figures
# ========================================================================
#' @title Generate All Manuscript Figures
#' @description Runs all experiments and saves Figures 1--5 as PDFs in the
#'   specified output directory.
#' @param output_dir Path to directory for saved figures.
#' @param seed Random seed for reproducibility.
#' @return Invisible NULL. Figures are saved as side effects.
#' @export
generate_all_figures <- function(output_dir = "figures_boosting", seed = 19671210) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(seed)

  # ---- Figure 1 ----
  cat("=== Figure 1: Weak Learnability ===\n")
  wl_results <- weak_learnability_experiment()
  p1 <- plot_weak_learnability(wl_results)
  ggsave(file.path(output_dir, "figure1_weak_learnability.pdf"),
         p1, width = 8, height = 5)
  cat("  Saved.\n")

  # ---- Shared data for Figures 2 & 3 ----
  cat("=== Preparing shared boosting + ACAR runs ===\n")
  set.seed(42)
  dat <- generate_classification_data(n = 200, p = 5, noise = 0.1)
  boost_res <- adaboost(dat$X, dat$y, M = 50)

  # Tighter quality gap + early_stop OFF so we see the full evolution
  # Higher noise so pheromone divergence is gradual, not instantaneous

  sq <- c(10, 8, 6, 4, 2)
  acar_res  <- acar(sq, n_ants = 30, n_waves = 50, noise_sd = 2.5,
                    early_stop = FALSE)

  # ---- Figure 2 ----
  cat("=== Figure 2: Weight vs Pheromone ===\n")
  p2 <- plot_weight_pheromone(boost_res, acar_res, sq)
  ggsave(file.path(output_dir, "figure2_weight_pheromone.pdf"),
         p2, width = 14, height = 5.5)
  cat("  Saved.\n")

  # ---- Figure 3 ----
  cat("=== Figure 3: Margin vs Quorum ===\n")
  p3 <- plot_margin_quorum(boost_res, dat$X, dat$y, acar_res)
  ggsave(file.path(output_dir, "figure3_margin_quorum.pdf"),
         p3, width = 14, height = 5.5)
  cat("  Saved.\n")

  # ---- Figure 4 ----
  cat("=== Figure 4: Convergence Rates ===\n")
  conv_results <- convergence_experiment_boost(n_boost_reps = 30, n_acar_reps = 200,
                                         max_iters = 80)
  p4 <- plot_convergence_boost(conv_results)
  ggsave(file.path(output_dir, "figure4_convergence.pdf"),
         p4, width = 8, height = 5)
  cat("  Saved.\n")

  # ---- Figure 5 ----
  cat("=== Figure 5: Noise Robustness ===\n")
  noise_results <- noise_experiment_boost()
  p5 <- plot_noise_robustness_boost(noise_results)
  ggsave(file.path(output_dir, "figure5_noise_robustness.pdf"),
         p5, width = 8, height = 5)
  cat("  Saved.\n")

  # ---- Performance summary table data ----
  cat("=== Table 1 data: Performance summary ===\n")
  perf_summ <- noise_results %>%
    filter(noise %in% c(0, 0.1, 0.2, 0.3, 0.4)) %>%
    group_by(noise, system) %>%
    summarise(mean_acc = mean(accuracy),
              sd_acc   = sd(accuracy), .groups = "drop")
  write.csv(perf_summ,
            file.path(output_dir, "table1_performance.csv"),
            row.names = FALSE)
  cat("  Saved table1_performance.csv\n")

  cat("\n=== All figures generated in:", output_dir, "===\n")
  invisible(NULL)
}
