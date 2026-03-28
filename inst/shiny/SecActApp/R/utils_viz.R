# Shared helpers for the SecAct Shiny app
# These are app-internal utilities, not exported from the SecAct package.

# --- Constants ---
DEFAULT_PLOT_WIDTH <- 600
DEFAULT_PLOT_HEIGHT <- 400
DIMENSION_ROUNDING <- 10
CELLTYPE_COLUMN <- "cellType"

MARKER_COLORS <- c("#00FF00", "#FF0000", "#888888", "#0000FF")
MARKER_LABELS <- c("PanCK (Epithelial)", "CD45 (Immune)", "CD3 (T-cells)", "DAPI (Nuclei)")

# --- Helpers ---

#' Placeholder plot for missing data or error states
empty_state_plot <- function(label) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label) +
    ggplot2::theme_void()
}

#' Get responsive plot dimensions from Shiny client data, rounded to reduce redraws
get_responsive_dims <- function(outputId, session) {
  width <- session$clientData$output[[outputId]]$width
  height <- session$clientData$output[[outputId]]$height
  if (is.null(width) || is.null(height)) {
    width <- DEFAULT_PLOT_WIDTH
    height <- DEFAULT_PLOT_HEIGHT
  }
  width <- round(width / DIMENSION_ROUNDING) * DIMENSION_ROUNDING
  height <- round(height / DIMENSION_ROUNDING) * DIMENSION_ROUNDING
  list(width = width, height = height)
}

#' CosMx labels tumor subtypes as "tumor 1", "tumor 2", etc. — collapse to single label
normalize_tumor_labels <- function(cell_types) {
  gsub("^tumor\\s+\\d+$", "tumor", cell_types)
}

#' Radius "0" means target cell itself; all others are neighbor distances in micrometers
format_radius_label <- function(radius) {
  if (radius == "0") "Target" else paste0(radius, " \u03bcm")
}
