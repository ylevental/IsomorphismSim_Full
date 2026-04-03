#' Launch an Interactive Shiny App
#'
#' Opens an interactive Shiny application for exploring the ant colony /
#' machine learning isomorphisms.
#'
#' @param part Which part of the trilogy to explore: \code{"part1"} for the
#'   random forest isomorphism (decorrelation, variance decomposition) or
#'   \code{"part2"} for the boosting isomorphism (adaptive recruitment,
#'   weak learnability).
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#' @export
#' @examples
#' if (interactive()) launch_app("part1")
#' if (interactive()) launch_app("part2")
launch_app <- function(part = c("part1", "part2"), ...) {
  part <- match.arg(part)
  app_dir <- system.file(paste0("shiny_", part), package = "AntsNet")
  if (app_dir == "") {
    stop("Shiny app for ", part, " not found. ",
         "Try re-installing the AntsNet package.", call. = FALSE)
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install it with install.packages('shiny').",
         call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
