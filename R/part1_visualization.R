#' @import ggplot2
#' @importFrom dplyr group_by summarise first
#' @importFrom viridis scale_fill_viridis
NULL

# Internal theme
theme_iso <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle    = element_text(hjust = 0.5, color = "grey40"),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}

#' Figure 1: Isomorphism Schematic
#'
#' Plots theoretical \eqn{\rho = \rho_{max}(1-\theta)} alongside empirical
#' values from both random forest and ant colony experiments.
#'
#' @param rf_results Output of \code{\link{variance_decomposition_experiment}}
#' @param colony_results Output of \code{\link{colony_variance_experiment}}
#' @return A ggplot object
#' @export
create_isomorphism_schematic <- function(rf_results, colony_results) {
  rf_sum <- rf_results |>
    dplyr::group_by(m_try) |>
    dplyr::summarise(rho_mean = mean(rho, na.rm = TRUE),
                     m_try_ratio = dplyr::first(m_try) / 50, .groups = "drop")

  col_cor <- colony_results |>
    dplyr::group_by(p_explore) |>
    dplyr::summarise(cor_avg = mean(cor_mean, na.rm = TRUE), .groups = "drop")

  rho_max <- max(c(rf_sum$rho_mean, col_cor$cor_avg), na.rm = TRUE)
  th <- seq(0, 1, length.out = 100)

  d <- rbind(
    data.frame(theta = th, correlation = rho_max * (1 - th), system = "Theoretical"),
    data.frame(theta = 1 - rf_sum$m_try_ratio, correlation = rf_sum$rho_mean,
               system = "Random Forest"),
    data.frame(theta = col_cor$p_explore, correlation = col_cor$cor_avg,
               system = "Ant Colony")
  )

  ggplot(d, aes(theta, correlation, color = system, linetype = system,
                shape = system)) +
    geom_line(data = subset(d, system == "Theoretical"), linewidth = 1.2) +
    geom_point(data = subset(d, system != "Theoretical"), size = 3.5) +
    geom_line(data = subset(d, system != "Theoretical"), linewidth = 0.8) +
    scale_color_manual(values = c("Theoretical" = "black",
                                   "Random Forest" = "#2E7D32",
                                   "Ant Colony" = "#BF360C")) +
    scale_linetype_manual(values = c("Theoretical" = "solid",
                                      "Random Forest" = "dashed",
                                      "Ant Colony" = "dotted")) +
    scale_shape_manual(values = c("Theoretical" = NA,
                                   "Random Forest" = 16,
                                   "Ant Colony" = 17)) +
    guides(color = guide_legend(override.aes = list(
             shape = c(17, 16, NA),
             linetype = c("dotted", "dashed", "solid")
           ))) +
    labs(title = "The Isomorphism: Correlation Decay with Decorrelation Parameter",
         subtitle = expression(rho == rho[max](1 - theta)),
         x = expression(theta ~ "(1 - m_try/p  or  p_explore)"),
         y = expression("Pairwise Correlation " ~ rho)) +
    theme_iso()
}

#' Figure 2: Variance Decomposition
#'
#' Side-by-side panels showing (A) correlation vs \code{m_try/p} and
#' (B) theoretical vs empirical ensemble variance.
#'
#' @param rf_results Output of \code{\link{variance_decomposition_experiment}}
#' @return A combined ggplot (via \code{ggpubr::ggarrange})
#' @export
plot_variance_decomposition <- function(rf_results) {
  s <- rf_results |>
    dplyr::group_by(m_try) |>
    dplyr::summarise(rho_mean = mean(rho), rho_sd = sd(rho),
                     ensemble_var_mean = mean(ensemble_var),
                     theoretical_var_mean = mean(theoretical_var),
                     m_try_ratio = dplyr::first(m_try) / 50, .groups = "drop")

  p1 <- ggplot(s, aes(m_try_ratio, rho_mean)) +
    geom_point(size = 3, color = "#2E7D32") +
    geom_line(color = "#2E7D32") +
    geom_errorbar(aes(ymin = rho_mean - 2 * rho_sd,
                      ymax = rho_mean + 2 * rho_sd), width = .02, color = "#2E7D32") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                linetype = "dashed", color = "red") +
    labs(title = "A. Correlation vs Feature Subsampling",
         x = expression(m[try] / p), y = expression(rho)) +
    theme_iso()

  p2 <- ggplot(s, aes(ensemble_var_mean, theoretical_var_mean)) +
    geom_point(size = 3, color = "#1565C0") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = "B. Theoretical vs Empirical Variance",
         x = "Empirical Ensemble Variance",
         y = expression(rho * sigma^2 + (1 - rho) * sigma^2 / M)) +
    theme_iso()

  ggpubr::ggarrange(p1, p2, ncol = 2)
}

#' Figure 3: Correlation Decay Comparison
#'
#' @param iso_results Output of \code{\link{isomorphism_test}}
#' @return A ggplot object
#' @export
plot_correlation_decay <- function(iso_results) {
  plot_data <- iso_results
  plot_data$theta[plot_data$system == "Random Forest"] <-
    1 - plot_data$theta[plot_data$system == "Random Forest"]
  plot_data$cor_sd[is.na(plot_data$cor_sd)] <- 0

  ggplot(plot_data, aes(theta, correlation, color = system,
                         shape = system, fill = system)) +
    geom_point(size = 3.5) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(ymin = correlation - 1.96 * cor_sd,
                    ymax = correlation + 1.96 * cor_sd),
                alpha = .15, color = NA) +
    scale_color_manual(values = c("Random Forest" = "#2E7D32",
                                   "Ant Colony" = "#BF360C")) +
    scale_fill_manual(values = c("Random Forest" = "#2E7D32",
                                  "Ant Colony" = "#BF360C")) +
    scale_shape_manual(values = c("Random Forest" = 16, "Ant Colony" = 17)) +
    labs(title = "Isomorphism Validation: Correlation Decay",
         x = expression(theta ~ "(1 - m_try/p  or  p_explore)"),
         y = expression(rho)) +
    theme_iso()
}

#' Figure 4: Optimal Decorrelation
#'
#' @param optimal_results Output of \code{\link{optimal_decorrelation_experiment}}
#' @return A ggplot object
#' @export
plot_optimal_decorrelation <- function(optimal_results) {
  ggplot(optimal_results, aes(theta, performance, color = factor(M))) +
    geom_point(size = 2) + geom_line(linewidth = .7) +
    geom_errorbar(aes(ymin = performance - 1.96 * se,
                      ymax = performance + 1.96 * se), width = .02, alpha = .5) +
    facet_wrap(~ system, scales = "free_y") +
    viridis::scale_color_viridis(discrete = TRUE, option = "D", end = .85) +
    labs(title = "Optimal Decorrelation: Performance vs Theta",
         x = expression(theta), y = "Error (lower = better)",
         color = "Ensemble Size") +
    theme_iso() + theme(legend.position = "right")
}

#' Figure 5: Sensitivity Heat-map
#'
#' @param sensitivity_results Output of \code{\link{sensitivity_analysis}}
#' @return A ggplot object
#' @export
plot_sensitivity_heatmap <- function(sensitivity_results) {
  d <- subset(sensitivity_results, system == "Random Forest")
  ggplot(d, aes(factor(signal), factor(noise), fill = correlation)) +
    geom_tile(color = "white", linewidth = .5) +
    geom_text(aes(label = sprintf("%.2f", correlation)), color = "white", size = 3.5) +
    facet_wrap(~ paste("theta =", theta), ncol = 3) +
    viridis::scale_fill_viridis(option = "C", direction = -1,
                                 name = expression(rho)) +
    labs(title = "Sensitivity: Tree Correlation Across Conditions",
         x = "Signal Strength", y = "Noise Level") +
    theme_iso() + theme(legend.position = "right")
}

#' Supplementary: Colony Accuracy vs Size
#'
#' @param colony_results Output of \code{\link{colony_variance_experiment}}
#' @return A ggplot object
#' @export
plot_colony_accuracy <- function(colony_results) {
  ggplot(colony_results, aes(n_ants, accuracy_mean,
                              color = factor(p_explore))) +
    geom_point(size = 3) + geom_line(linewidth = .7) +
    viridis::scale_color_viridis(discrete = TRUE, option = "D", end = .9) +
    labs(title = "Colony Decision Accuracy vs Colony Size",
         x = "Number of Ants (N)", y = "P(Correct Decision)",
         color = expression(p[explore])) +
    theme_iso() + theme(legend.position = "right")
}
