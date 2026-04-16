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

#' Swap a SpaCET object's counts slot with an activity matrix for visualization.
#' SpaCET.visualize.spatialFeature() renders from @input$counts, so we temporarily
#' replace it with the activity matrix to reuse the existing rendering pipeline.
swap_activity_matrix <- function(spacet_obj, activity_matrix) {
  common_spots <- intersect(colnames(activity_matrix), colnames(spacet_obj@input$counts))
  temp_mat <- Matrix::Matrix(0, nrow = nrow(activity_matrix), ncol = ncol(spacet_obj@input$counts),
                             sparse = TRUE,
                             dimnames = list(rownames(activity_matrix), colnames(spacet_obj@input$counts)))
  temp_mat[, common_spots] <- activity_matrix[, common_spots]
  spacet_obj@input$counts <- temp_mat
  spacet_obj
}

#' Extract a platform zip upload and find the output directory.
#' Handles nested zip structures by searching for a landmark file pattern.
#' @param upload_input The Shiny fileInput value (must have $datapath)
#' @param platform_name Label for notifications (e.g., "Space Ranger", "CosMx")
#' @param landmark_pattern Regex to find the platform's landmark file/dir
#' @param landmark_type "dir" to search for directories, "file" for files
#' @return Path to the platform output directory, or NULL on failure
extract_platform_zip <- function(upload_input, platform_name, landmark_pattern, landmark_type = "file") {
  zip_path <- upload_input$datapath
  extract_dir <- file.path(tempdir(), paste0(tolower(gsub(" ", "_", platform_name)), "_upload"))
  if (dir.exists(extract_dir)) unlink(extract_dir, recursive = TRUE)
  on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)

  utils::unzip(zip_path, exdir = extract_dir)

  if (landmark_type == "dir") {
    matches <- list.files(extract_dir, pattern = landmark_pattern,
                          recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
  } else {
    matches <- list.files(extract_dir, pattern = landmark_pattern,
                          recursive = TRUE, full.names = TRUE)
  }

  if (length(matches) == 0) {
    shiny::showNotification(
      paste0("No '", landmark_pattern, "' found in zip. Is this ", platform_name, " output?"),
      type = "error"
    )
    return(NULL)
  }

  dirname(matches[1])
}
