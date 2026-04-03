library(shiny)
library(ggplot2)
library(viridis)

# ── Helper: ant colony simulation (self-contained for the app) ──────────
sim_colony <- function(n_ants, n_sites, site_qualities, p_explore,
                       n_steps, noise_sd, quorum_threshold,
                       alpha = 0.1, beta = 0.05, gamma = 0.2) {
  ant_states      <- matrix(0, n_ants, n_sites)
  ant_preferences <- matrix(0, n_ants, n_sites)
  pheromone       <- rep(0, n_sites)
  recruitment     <- rep(0, n_sites)
  decision        <- NA
  phero_history   <- matrix(NA, n_steps, n_sites)
  recruit_history <- matrix(NA, n_steps, n_sites)

  for (t in seq_len(n_steps)) {
    for (i in seq_len(n_ants)) {
      if (runif(1) < p_explore) {
        site <- sample(n_sites, 1)
      } else {
        if (sum(pheromone) > 0) {
          site <- sample(n_sites, 1, prob = pheromone / sum(pheromone))
        } else {
          site <- sample(n_sites, 1)
        }
      }
      obs <- site_qualities[site] + rnorm(1, 0, noise_sd)
      ant_states[i, site] <- ant_states[i, site] + 1
      nv <- ant_states[i, site]
      ant_preferences[i, site] <- (ant_preferences[i, site] * (nv - 1) + obs) / nv
      if (ant_preferences[i, site] > mean(site_qualities)) {
        recruitment[site] <- recruitment[site] +
          gamma * (ant_preferences[i, site] - mean(site_qualities))
      }
    }
    pheromone <- pheromone * (1 - beta) + recruitment * alpha
    phero_history[t, ]   <- pheromone
    recruit_history[t, ] <- recruitment
    if (max(recruitment) > quorum_threshold) { decision <- which.max(recruitment); break }
    recruitment <- recruitment * 0.8
  }
  if (is.na(decision)) decision <- which.max(recruitment)
  list(decision = decision, ant_preferences = ant_preferences,
       phero_history = phero_history[seq_len(t), , drop = FALSE],
       recruit_history = recruit_history[seq_len(t), , drop = FALSE],
       steps_used = t)
}

# ── Theme ───────────────────────────────────────────────────────────────
theme_app <- function(bs = 13) {
  theme_minimal(base_size = bs) +
    theme(plot.title = element_text(face = "bold", hjust = .5),
          legend.position = "bottom", panel.grid.minor = element_blank())
}

# ═══════════════════════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════════════════════
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
    .navbar { background: linear-gradient(135deg, #1b5e20 0%, #4e342e 100%); }
    .navbar-brand { color: #fff !important; font-weight: bold; }
    .nav-link { color: #e0e0e0 !important; }
    .well { background: #fafafa; border: 1px solid #ddd; }
    h4 { color: #2e7d32; }
  "))),

  navbarPage(
    title = "Isomorphism Explorer: Ants \u2194 Ensemble Learning",
    id = "main_nav",

    # ── Tab 1: Single Colony ─────────────────────────────────────────────
    tabPanel("Colony Simulation",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Colony Parameters"),
          sliderInput("n_ants", "Number of ants (N)", 10, 200, 50, step = 10),
          sliderInput("p_explore", "Exploration probability", 0, 1, 0.3, step = 0.05),
          sliderInput("noise_sd", "Observation noise (\u03c3)", 0.5, 5, 2, step = 0.5),
          sliderInput("n_steps", "Max time steps", 50, 500, 100, step = 50),
          numericInput("quorum", "Quorum threshold", 20, min = 5, max = 100),
          hr(),
          h4("Site Qualities"),
          sliderInput("q1", "Site 1 (best)", 1, 20, 10),
          sliderInput("q2", "Site 2", 1, 20, 8),
          sliderInput("q3", "Site 3", 1, 20, 6),
          sliderInput("q4", "Site 4", 1, 20, 4),
          sliderInput("q5", "Site 5 (worst)", 1, 20, 2),
          actionButton("run_colony", "Run Simulation",
                       class = "btn-success btn-block",
                       style = "width:100%; margin-top:10px;")
        ),
        mainPanel(width = 9,
          fluidRow(
            column(6, plotOutput("phero_plot", height = "350px")),
            column(6, plotOutput("recruit_plot", height = "350px"))
          ),
          fluidRow(
            column(6, plotOutput("pref_plot", height = "350px")),
            column(6, verbatimTextOutput("colony_summary"))
          )
        )
      )
    ),

    # ── Tab 2: Isomorphism Sweep ─────────────────────────────────────────
    tabPanel("Isomorphism Sweep",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Sweep Parameters"),
          sliderInput("sweep_reps", "Replicates per \u03b8", 5, 50, 15, step = 5),
          sliderInput("sweep_ants", "Ants per colony", 10, 100, 50, step = 10),
          checkboxGroupInput("sweep_theta", "Theta values",
                             choices = as.character(seq(0.1, 0.9, 0.1)),
                             selected = as.character(c(0.1, 0.3, 0.5, 0.7, 0.9))),
          actionButton("run_sweep", "Run Sweep",
                       class = "btn-success btn-block",
                       style = "width:100%; margin-top:10px;"),
          hr(),
          helpText("This runs ant colony simulations at each \u03b8 and",
                   "plots the resulting correlation decay — the core",
                   "of the isomorphism with random forests.")
        ),
        mainPanel(width = 9,
          plotOutput("sweep_plot", height = "500px"),
          verbatimTextOutput("sweep_summary")
        )
      )
    ),

    # ── Tab 3: Variance Decomposition ────────────────────────────────────
    tabPanel("Variance Decomposition",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Demonstration"),
          helpText("Shows how the theoretical variance formula",
                   "Var = \u03c1\u03c3\u00b2 + (1\u2212\u03c1)\u03c3\u00b2/M",
                   "predicts the actual ensemble variance."),
          sliderInput("vd_rho", "Correlation (\u03c1)", 0, 1, 0.5, step = 0.05),
          sliderInput("vd_sigma2", "Unit variance (\u03c3\u00b2)", 0.1, 5, 1, step = 0.1),
          sliderInput("vd_M_max", "Max ensemble size", 10, 500, 200, step = 10)
        ),
        mainPanel(width = 9,
          plotOutput("vd_plot", height = "500px"),
          helpText("Red dashed line = irreducible variance floor at \u03c1\u03c3\u00b2.",
                   "As M \u2192 \u221e, the ensemble variance converges to this floor.",
                   "Lower \u03c1 (more decorrelation) makes the floor lower.")
        )
      )
    ),

    # ── Tab 4: About ─────────────────────────────────────────────────────
    tabPanel("About",
      fluidRow(column(8, offset = 2,
        h3("Decorrelation, Diversity, and Emergent Intelligence"),
        p("This app accompanies the paper:"),
        tags$blockquote(
          tags$em("The Isomorphism Between Social Insect Colonies",
                  "and Ensemble Machine Learning"),
          tags$br(),
          "Ernest Fokou\u00e9, Gregory Babbitt, Yuval Levental",
          tags$br(), "Rochester Institute of Technology, 2026"
        ),
        h4("Key Insight"),
        p("Ant colonies and random forests are instances of the same",
          "abstract computational system. Both reduce variance by",
          "decorrelating identical units:"),
        tags$ul(
          tags$li(tags$b("Random Forest:"),
                  " bootstrap + random feature selection decorrelate trees"),
          tags$li(tags$b("Ant Colony:"),
                  " stochastic exploration decorrelates ant assessments")
        ),
        p("The variance of the ensemble is governed by the same formula",
          "in both cases:"),
        tags$p(style = "text-align:center; font-size:1.2em;",
               "Var = \u03c1\u03c3\u00b2 + (1 \u2212 \u03c1)\u03c3\u00b2 / M"),
        h4("R Package"),
        p("Install the companion package with:"),
        tags$code("devtools::install('IsomorphismSim')")
      ))
    )
  )
)

# ═══════════════════════════════════════════════════════════════════════════
# SERVER
# ═══════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Tab 1: Colony simulation ───────────────────────────────────────────
  colony_result <- eventReactive(input$run_colony, {
    sq <- c(input$q1, input$q2, input$q3, input$q4, input$q5)
    sim_colony(n_ants = input$n_ants, n_sites = 5, site_qualities = sq,
               p_explore = input$p_explore, n_steps = input$n_steps,
               noise_sd = input$noise_sd, quorum_threshold = input$quorum)
  })

  output$phero_plot <- renderPlot({
    res <- colony_result()
    ph <- as.data.frame(res$phero_history)
    colnames(ph) <- paste("Site", 1:5)
    ph$step <- seq_len(nrow(ph))
    long <- tidyr::pivot_longer(ph, -step, names_to = "site", values_to = "pheromone")
    ggplot(long, aes(step, pheromone, color = site)) +
      geom_line(linewidth = 1) + scale_color_viridis_d() +
      labs(title = "Pheromone Levels Over Time", x = "Step", y = "Pheromone") +
      theme_app()
  })

  output$recruit_plot <- renderPlot({
    res <- colony_result()
    rh <- as.data.frame(res$recruit_history)
    colnames(rh) <- paste("Site", 1:5)
    rh$step <- seq_len(nrow(rh))
    long <- tidyr::pivot_longer(rh, -step, names_to = "site", values_to = "recruitment")
    ggplot(long, aes(step, recruitment, color = site)) +
      geom_line(linewidth = 1) + scale_color_viridis_d() +
      labs(title = "Recruitment Over Time", x = "Step", y = "Recruitment") +
      theme_app()
  })

  output$pref_plot <- renderPlot({
    res <- colony_result()
    prefs <- as.data.frame(res$ant_preferences)
    colnames(prefs) <- paste("Site", 1:5)
    prefs$ant <- seq_len(nrow(prefs))
    long <- tidyr::pivot_longer(prefs, -ant, names_to = "site", values_to = "quality_est")
    ggplot(long, aes(site, quality_est, fill = site)) +
      geom_boxplot(alpha = 0.7) + scale_fill_viridis_d() +
      labs(title = "Ant Quality Estimates by Site", x = NULL, y = "Estimated Quality") +
      theme_app() + theme(legend.position = "none")
  })

  output$colony_summary <- renderText({
    res <- colony_result()
    sq <- c(input$q1, input$q2, input$q3, input$q4, input$q5)
    paste0(
      "=== Colony Decision ===\n",
      "Chosen site: ", res$decision,
      " (true best: ", which.max(sq), ")\n",
      "Correct: ", ifelse(res$decision == which.max(sq), "YES", "NO"), "\n",
      "Steps used: ", res$steps_used, " / ", input$n_steps, "\n\n",
      "Colony avg quality estimates:\n",
      paste(sprintf("  Site %d: %.2f  (true: %d)", 1:5,
                    colMeans(res$ant_preferences), sq), collapse = "\n")
    )
  })

  # ── Tab 2: Isomorphism sweep ───────────────────────────────────────────
  sweep_result <- eventReactive(input$run_sweep, {
    thetas <- sort(as.numeric(input$sweep_theta))
    n_reps <- input$sweep_reps
    n_a    <- input$sweep_ants
    results <- data.frame()

    withProgress(message = "Running sweep...", value = 0, {
      for (j in seq_along(thetas)) {
        th <- thetas[j]
        pmat <- matrix(NA, n_a, n_reps)
        for (r in seq_len(n_reps)) {
          sim <- sim_colony(n_ants = n_a, n_sites = 5,
                            site_qualities = c(10, 8, 6, 4, 2),
                            p_explore = th, n_steps = 100,
                            noise_sd = 2, quorum_threshold = 20)
          pmat[, r] <- sim$ant_preferences[, 1]
        }
        cm <- cor(pmat, use = "pairwise.complete.obs")
        rho <- mean(cm[lower.tri(cm)], na.rm = TRUE)
        results <- rbind(results, data.frame(theta = th, rho = rho))
        incProgress(1 / length(thetas))
      }
    })
    results
  })

  output$sweep_plot <- renderPlot({
    res <- sweep_result()
    rho_max <- max(res$rho, na.rm = TRUE)
    th <- seq(0, 1, length.out = 100)
    theory <- data.frame(theta = th, rho = rho_max * (1 - th), type = "Theoretical")
    emp    <- data.frame(theta = res$theta, rho = res$rho, type = "Ant Colony (empirical)")
    d <- rbind(theory, emp)

    ggplot(d, aes(theta, rho, color = type, linetype = type)) +
      geom_line(data = subset(d, type == "Theoretical"), linewidth = 1.2) +
      geom_point(data = subset(d, type != "Theoretical"), size = 4) +
      geom_line(data = subset(d, type != "Theoretical"), linewidth = 0.8) +
      scale_color_manual(values = c("Theoretical" = "black",
                                     "Ant Colony (empirical)" = "#BF360C")) +
      scale_linetype_manual(values = c("Theoretical" = "solid",
                                        "Ant Colony (empirical)" = "dashed")) +
      labs(title = "Isomorphism Sweep: Correlation Decay with Exploration Probability",
           subtitle = expression("Theoretical:" ~ rho == rho[max](1 - theta)),
           x = expression(theta ~ "(p_explore)"),
           y = expression("Pairwise Correlation " ~ rho)) +
      theme_app(14)
  })

  output$sweep_summary <- renderText({
    res <- sweep_result()
    fit <- lm(rho ~ theta, data = res)
    paste0(
      "Linear fit:  rho = ", sprintf("%.3f", coef(fit)[1]),
      " + (", sprintf("%.3f", coef(fit)[2]), ") * theta\n",
      "R-squared: ", sprintf("%.4f", summary(fit)$r.squared), "\n",
      "Expected slope (negative) confirms decorrelation mechanism."
    )
  })

  # ── Tab 3: Variance decomposition demo ─────────────────────────────────
  output$vd_plot <- renderPlot({
    rho <- input$vd_rho
    s2  <- input$vd_sigma2
    M   <- seq(1, input$vd_M_max, by = 1)
    v   <- rho * s2 + (1 - rho) * s2 / M
    d   <- data.frame(M = M, Variance = v)

    ggplot(d, aes(M, Variance)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = rho * s2, linetype = "dashed",
                 color = "red", linewidth = 0.8) +
      annotate("text", x = max(M) * 0.75, y = rho * s2 + s2 * 0.05,
               label = sprintf("Floor = \u03c1\u03c3\u00b2 = %.3f", rho * s2),
               color = "red", size = 5) +
      labs(title = "Ensemble Variance as a Function of Ensemble Size",
           subtitle = sprintf("\u03c1 = %.2f,  \u03c3\u00b2 = %.1f", rho, s2),
           x = "Ensemble Size (M)", y = "Var[ensemble]") +
      ylim(0, NA) +
      theme_app(14)
  })
}

shinyApp(ui, server)
