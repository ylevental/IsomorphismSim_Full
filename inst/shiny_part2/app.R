# ===========================================================================
# AntBoost Interactive Explorer
# inst/shiny/app.R
#
# Run with:
#   shiny::runApp("inst/shiny")
# ===========================================================================

library(shiny)
library(ggplot2)

# Source package functions if not installed
if (!requireNamespace("AntBoost", quietly = TRUE)) {
  pkg_dir <- if (file.exists("../../R/adaboost.R")) "../../R" else "R"
  for (f in list.files(pkg_dir, pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }
}

ui <- fluidPage(
  titlePanel("AntBoost: Exploring the Boosting \u2194 Ant Colony Isomorphism"),

  sidebarLayout(
    sidebarPanel(
      h4("Ant Colony Parameters"),
      sliderInput("n_ants", "Number of Ants per Wave:", 10, 200, 50, step = 10),
      sliderInput("n_waves", "Number of Waves:", 10, 200, 100, step = 10),
      sliderInput("rho", "Evaporation Rate (\u03C1):", 0.01, 0.5, 0.1, step = 0.01),
      sliderInput("gamma_dep", "Deposition Rate (\u03B3):", 0.1, 2.0, 0.5, step = 0.1),
      sliderInput("noise_sd", "Observation Noise \u03C3:", 0, 3, 0.5, step = 0.1),
      hr(),
      h4("Site Qualities"),
      textInput("qualities", "Comma-separated:", "10, 5, 3, 2, 1"),
      hr(),
      actionButton("run", "Run Simulation", class = "btn-primary btn-lg"),
      width = 3
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Pheromone Evolution",
                 plotOutput("pheromone_plot", height = "500px")),
        tabPanel("Fitness History",
                 plotOutput("fitness_plot", height = "500px")),
        tabPanel("Colony Decision",
                 plotOutput("decision_plot", height = "500px")),
        tabPanel("Quorum Margin",
                 plotOutput("quorum_plot", height = "500px"),
                 verbatimTextOutput("quorum_text"))
      ),
      width = 9
    )
  )
)

server <- function(input, output) {
  sim_result <- eventReactive(input$run, {
    sq <- as.numeric(strsplit(input$qualities, ",\\s*")[[1]])
    acar(
      site_qualities = sq,
      n_ants    = input$n_ants,
      n_waves   = input$n_waves,
      rho       = input$rho,
      gamma     = input$gamma_dep,
      noise_sd  = input$noise_sd,
      early_stop = FALSE
    )
  })

  output$pheromone_plot <- renderPlot({
    res <- sim_result()
    ph  <- res$pheromone_history
    K   <- ncol(ph)
    df  <- data.frame(
      wave      = rep(seq_len(nrow(ph)), K),
      pheromone = as.vector(ph),
      site      = factor(rep(seq_len(K), each = nrow(ph)))
    )
    ggplot(df, aes(x = wave, y = pheromone, color = site)) +
      geom_line(linewidth = 1.2) +
      labs(title = "Pheromone Evolution Across Recruitment Waves",
           x = "Wave", y = expression(tau[j](t)), color = "Site") +
      theme_minimal(base_size = 14)
  })

  output$fitness_plot <- renderPlot({
    res <- sim_result()
    df  <- data.frame(wave = seq_along(res$fitness_history),
                      fitness = res$fitness_history)
    ggplot(df, aes(x = wave, y = fitness)) +
      geom_line(color = "#8C564B", linewidth = 1.2) +
      geom_smooth(method = "loess", color = "red", se = FALSE, linewidth = 0.8) +
      labs(title = "Colony Fitness Over Time", x = "Wave", y = "Fitness") +
      theme_minimal(base_size = 14)
  })

  output$decision_plot <- renderPlot({
    res <- sim_result()
    df  <- data.frame(wave = seq_along(res$decision_history),
                      decision = factor(res$decision_history))
    ggplot(df, aes(x = wave, y = decision)) +
      geom_point(alpha = 0.5, size = 2) +
      labs(title = "Colony Decision per Wave", x = "Wave", y = "Chosen Site") +
      theme_minimal(base_size = 14)
  })

  output$quorum_plot <- renderPlot({
    res <- sim_result()
    final_ph <- res$pheromone_history[nrow(res$pheromone_history), ]
    best <- res$final_decision
    df <- data.frame(
      site = factor(seq_along(final_ph)),
      pheromone = final_ph,
      is_best = seq_along(final_ph) == best
    )
    ggplot(df, aes(x = site, y = pheromone, fill = is_best)) +
      geom_col(alpha = 0.8, color = "white") +
      scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4472C4"),
                        guide = "none") +
      labs(title = "Final Pheromone Distribution", x = "Site",
           y = expression(tau[j](T))) +
      theme_minimal(base_size = 14)
  })

  output$quorum_text <- renderPrint({
    res <- sim_result()
    qm  <- calculate_quorum_margin(res)
    cat(sprintf("Final decision: Site %d\n", res$final_decision))
    cat(sprintf("Quorum margin:  %.4f\n", qm))
    cat(sprintf("Waves used:     %d\n", nrow(res$pheromone_history)))
  })
}

shinyApp(ui = ui, server = server)
