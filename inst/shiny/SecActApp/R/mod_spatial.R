# Spatial Visualization Module
# Displays SecAct activity, cell types, gene expression on spatial coordinates.
# Delegates rendering to SpaCET.visualize.spatialFeature() for multi-platform support.

source(file.path("R", "utils_viz.R"), local = TRUE)

# --- UI ---
spatialUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,

        # Dataset loading
        shiny::wellPanel(
          shiny::h4("Dataset", style = paste0("color: ", UI_COLORS$primary)),
          shiny::fileInput(ns("dataUpload"), "Upload SpaCET/Seurat Object (.RDS, .h5seurat)",
                           accept = c(".rds", ".RDS", ".h5seurat")),
          shiny::actionButton(ns("loadBtn"), "Load Dataset",
                              class = "btn-primary btn-block")
        ),

        # Visualization controls (visible after data loaded)
        shiny::conditionalPanel(
          condition = paste0("output['", ns("dataLoaded"), "']"),
          shiny::wellPanel(
            shiny::h4("Visualization", style = paste0("color: ", UI_COLORS$primary)),

            shiny::selectInput(ns("spatialType"), "Display Type:",
                               choices = c("SecAct Activity" = "SecActActivity",
                                           "Gene Expression" = "GeneExpression",
                                           "Cell Fraction" = "CellFraction",
                                           "Cell Type Composition" = "CellTypeComposition")),

            shinyWidgets::pickerInput(
              ns("feature"), "Select Feature:",
              choices = NULL,
              options = list(`live-search` = TRUE, size = 10),
              multiple = FALSE
            ),

            shiny::sliderInput(ns("pointSize"), "Point Size:",
                               min = 0.1, max = 3, value = 1, step = 0.1),

            shiny::checkboxInput(ns("imageBg"), "Show Tissue Image", value = TRUE)
          )
        )
      ),

      shiny::mainPanel(
        width = 9,

        # Welcome screen (before data loaded)
        shiny::conditionalPanel(
          condition = paste0("!output['", ns("dataLoaded"), "']"),
          shiny::div(
            style = "height: 500px; display: flex; justify-content: center; align-items: center; background-color: #f8f9fa; border-radius: 5px;",
            shiny::div(
              style = "text-align: center; max-width: 600px;",
              shiny::h2("Spatial Visualization"),
              shiny::p("Upload a SpaCET or Seurat object to explore spatial secreted protein activity."),
              shiny::p("Supports: Visium, VisiumHD, CosMx, Xenium, Slide-Seq",
                       style = "color: #666;")
            )
          )
        ),

        # Visualization panel (after data loaded)
        shiny::conditionalPanel(
          condition = paste0("output['", ns("dataLoaded"), "']"),
          shiny::fluidRow(
            shiny::column(12,
              shiny::div(
                style = "border: 1px solid #ddd; border-radius: 5px; padding: 10px;",
                shiny::plotOutput(ns("spatialPlot"), height = "600px")
              )
            )
          ),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(12,
              DT::dataTableOutput(ns("dataTable"))
            )
          )
        )
      )
    )
  )
}

# --- Server ---
spatialServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Reactive values
    rv <- shiny::reactiveValues(
      spacet_obj = NULL,
      dataLoaded = FALSE
    )

    # Flag for conditionalPanel
    output$dataLoaded <- shiny::reactive({ rv$dataLoaded })
    shiny::outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)

    # Load dataset
    shiny::observeEvent(input$loadBtn, {
      shiny::req(input$dataUpload)

      file_path <- input$dataUpload$datapath
      file_name <- input$dataUpload$name
      ext <- tolower(tools::file_ext(file_name))

      tryCatch({
        shiny::withProgress(message = "Loading dataset...", {
          if (ext == "h5seurat") {
            # Load h5seurat via SeuratDisk, convert to SpaCET
            if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
              shiny::showNotification("SeuratDisk package required for .h5seurat files", type = "error")
              return()
            }
            seu <- SeuratDisk::LoadH5Seurat(file_path)
            rv$spacet_obj <- SpaCET::convert.Seurat(seu)
          } else {
            # Load RDS — could be SpaCET or Seurat object
            obj <- readRDS(file_path)
            if (inherits(obj, "SpaCET")) {
              rv$spacet_obj <- obj
            } else if (inherits(obj, "Seurat")) {
              rv$spacet_obj <- SpaCET::convert.Seurat(obj)
            } else {
              shiny::showNotification("File must contain a SpaCET or Seurat object", type = "error")
              return()
            }
          }

          rv$dataLoaded <- TRUE

          # Update feature choices based on spatial type
          update_features()

          shiny::showNotification("Dataset loaded successfully", type = "message")
        })
      }, error = function(e) {
        shiny::showNotification(paste("Error loading file:", e$message), type = "error")
      })
    })

    # Update feature picker when spatial type changes
    update_features <- function() {
      shiny::req(rv$spacet_obj)
      obj <- rv$spacet_obj

      choices <- switch(input$spatialType,
        "SecActActivity" = {
          if (!is.null(obj@results$SecAct_output$SecActTarget)) {
            rownames(obj@results$SecAct_output$SecActTarget)
          } else if (!is.null(obj@results$SecAct_output)) {
            # Try alternative slot names
            nms <- names(obj@results$SecAct_output)
            mat_name <- nms[grepl("Target|activity", nms, ignore.case = TRUE)][1]
            if (!is.na(mat_name)) rownames(obj@results$SecAct_output[[mat_name]]) else character(0)
          } else {
            character(0)
          }
        },
        "GeneExpression" = rownames(obj@input$counts),
        "CellFraction" = {
          if (!is.null(obj@results$deconvolution$propMat)) {
            rownames(obj@results$deconvolution$propMat)
          } else character(0)
        },
        "CellTypeComposition" = "All",
        character(0)
      )

      shinyWidgets::updatePickerInput(session, "feature",
                                       choices = choices,
                                       selected = if (length(choices) > 0) choices[1] else NULL)
    }

    shiny::observeEvent(input$spatialType, {
      if (rv$dataLoaded) update_features()
    })

    # Render spatial plot using SpaCET
    output$spatialPlot <- shiny::renderPlot({
      shiny::req(rv$spacet_obj, input$feature)

      tryCatch({
        # Map app spatial types to SpaCET spatial types
        spacet_type <- switch(input$spatialType,
          "SecActActivity" = "GeneExpression",  # SecAct stores activity like gene expression
          "GeneExpression" = "GeneExpression",
          "CellFraction" = "CellFraction",
          "CellTypeComposition" = "CellTypeComposition",
          "GeneExpression"
        )

        # For SecAct activity, temporarily swap the counts matrix
        if (input$spatialType == "SecActActivity") {
          # Get the activity matrix from SecAct results
          act_mat <- rv$spacet_obj@results$SecAct_output$SecActTarget
          if (is.null(act_mat)) {
            return(empty_state_plot("No SecAct activity data found in this object"))
          }

          # Create a temporary SpaCET object with activity as the "counts"
          temp_obj <- rv$spacet_obj
          # Store original counts, replace with activity matrix for visualization
          common_spots <- intersect(colnames(act_mat), colnames(temp_obj@input$counts))
          temp_mat <- Matrix::Matrix(0, nrow = nrow(act_mat), ncol = ncol(temp_obj@input$counts),
                                     sparse = TRUE,
                                     dimnames = list(rownames(act_mat), colnames(temp_obj@input$counts)))
          temp_mat[, common_spots] <- act_mat[, common_spots]
          temp_obj@input$counts <- temp_mat

          SpaCET::SpaCET.visualize.spatialFeature(
            temp_obj,
            spatialType = "GeneExpression",
            spatialFeatures = input$feature,
            scaleTypeForGeneExpression = "RawCounts",
            pointSize = input$pointSize,
            imageBg = input$imageBg
          )
        } else {
          SpaCET::SpaCET.visualize.spatialFeature(
            rv$spacet_obj,
            spatialType = spacet_type,
            spatialFeatures = input$feature,
            pointSize = input$pointSize,
            imageBg = input$imageBg
          )
        }
      }, error = function(e) {
        empty_state_plot(paste("Visualization error:", e$message))
      })
    })

    # Data table
    output$dataTable <- DT::renderDataTable({
      shiny::req(rv$spacet_obj, input$feature, input$spatialType)

      tryCatch({
        obj <- rv$spacet_obj

        if (input$spatialType == "SecActActivity") {
          act_mat <- obj@results$SecAct_output$SecActTarget
          shiny::req(act_mat)
          if (!input$feature %in% rownames(act_mat)) return(NULL)
          values <- act_mat[input$feature, ]
        } else if (input$spatialType == "GeneExpression") {
          counts <- obj@input$counts
          if (!input$feature %in% rownames(counts)) return(NULL)
          values <- counts[input$feature, ]
        } else {
          return(NULL)
        }

        coords <- obj@input$spotCoordinates
        df <- data.frame(
          SpotID = names(values),
          Value = as.numeric(values),
          stringsAsFactors = FALSE
        )

        if (!is.null(coords)) {
          df$x <- coords[match(df$SpotID, rownames(coords)), "x"]
          df$y <- coords[match(df$SpotID, rownames(coords)), "y"]
        }

        if (!is.null(obj@input$metaData) && CELLTYPE_COLUMN %in% colnames(obj@input$metaData)) {
          df$CellType <- obj@input$metaData[match(df$SpotID, rownames(obj@input$metaData)), CELLTYPE_COLUMN]
          df$CellType <- normalize_tumor_labels(df$CellType)
        }

        DT::datatable(df, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
      }, error = function(e) {
        DT::datatable(data.frame(Error = e$message))
      })
    })
  })
}
