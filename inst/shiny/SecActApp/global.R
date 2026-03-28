# SecAct Unified Visualization App — Global Configuration

# Load app-internal helpers
source(file.path("R", "utils_viz.R"), local = TRUE)

# Shiny dependencies — checked by runSecActApp() launcher before we get here
library(shiny)
library(shinyWidgets)
library(ggplot2)
library(DT)
library(bslib)

# Optional packages — loaded if available
if (requireNamespace("shinyjs", quietly = TRUE)) library(shinyjs)
if (requireNamespace("shinyFeedback", quietly = TRUE)) library(shinyFeedback)

# SpaCET is required for spatial visualization
if (!requireNamespace("SpaCET", quietly = TRUE)) {
  stop("SpaCET package is required. Install from: https://github.com/data2intelligence/spacet")
}

# Maximum file upload size (3 GB)
options(shiny.maxRequestSize = 3000 * 1024^2)

# App-wide color palette
UI_COLORS <- list(
  primary = "#3498db",
  success = "#2ecc71",
  warning = "#f39c12",
  danger  = "#e74c3c",
  muted   = "#95a5a6"
)
