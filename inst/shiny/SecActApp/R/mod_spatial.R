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

          shiny::tabsetPanel(
            id = ns("loadTabs"),

            # Tab: Upload pre-built object
            shiny::tabPanel("Upload Object",
              shiny::div(style = "margin-top: 10px;",
                shiny::fileInput(ns("dataUpload"), "SpaCET/Seurat Object (.RDS, .h5seurat)",
                                 accept = c(".rds", ".RDS", ".h5seurat")),
                shiny::actionButton(ns("loadBtn"), "Load Dataset",
                                    class = "btn-primary btn-block")
              )
            ),

            # Tab: Upload Space Ranger output
            shiny::tabPanel("Space Ranger",
              shiny::div(style = "margin-top: 10px;",
                shiny::fileInput(ns("spacerUpload"), "Space Ranger Output (.zip)",
                                 accept = c(".zip")),
                shiny::actionButton(ns("loadSpaceRangerBtn"), "Load Visium Data",
                                    class = "btn-primary btn-block")
              )
            ),

            # Tab: Upload CosMx output
            shiny::tabPanel("CosMx",
              shiny::div(style = "margin-top: 10px;",
                shiny::fileInput(ns("cosmxUpload"), "CosMx Output (.zip)",
                                 accept = c(".zip")),
                shiny::actionButton(ns("loadCosmxBtn"), "Load CosMx Data",
                                    class = "btn-primary btn-block")
              )
            ),

            # Tab: Upload Xenium output
            shiny::tabPanel("Xenium",
              shiny::div(style = "margin-top: 10px;",
                shiny::fileInput(ns("xeniumUpload"), "Xenium Output (.zip)",
                                 accept = c(".zip")),
                shiny::actionButton(ns("loadXeniumBtn"), "Load Xenium Data",
                                    class = "btn-primary btn-block")
              )
            ),

            # Tab: Demo dataset
            shiny::tabPanel("Demo",
              shiny::div(style = "margin-top: 10px;",
                shiny::p("Load the bundled Visium HCC example dataset.",
                         style = "color: #666; font-size: 0.9em;"),
                shiny::actionButton(ns("loadDemoBtn"), "Load Visium HCC Demo",
                                    class = "btn-info btn-block",
                                    icon = shiny::icon("flask"))
              )
            )
          )
        ),

        # SecAct inference controls (visible after data loaded, before inference run)
        shiny::conditionalPanel(
          condition = paste0("output['", ns("dataLoaded"), "'] && !output['", ns("hasSecActResults"), "']"),
          shiny::wellPanel(
            shiny::h4("Run SecAct Inference", style = paste0("color: ", UI_COLORS$warning)),
            shiny::p("No SecAct activity results found. Run inference on this dataset.",
                     style = "color: #666; font-size: 0.9em;"),
            shiny::actionButton(ns("runInferenceBtn"), "Run SecAct Inference",
                                class = "btn-warning btn-block",
                                icon = shiny::icon("play"))
          )
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
              shiny::p("Upload a SpaCET or Seurat object, a Space Ranger output zip, or try the demo dataset."),
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

    # Flags for conditionalPanel
    output$dataLoaded <- shiny::reactive({ rv$dataLoaded })
    shiny::outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)

    output$hasSecActResults <- shiny::reactive({
      !is.null(rv$spacet_obj) && !is.null(rv$spacet_obj@results$SecAct_output)
    })
    shiny::outputOptions(output, "hasSecActResults", suspendWhenHidden = FALSE)

    # Helper: finalize after any load path succeeds
    finish_load <- function(obj, source_label) {
      rv$spacet_obj <- obj
      rv$dataLoaded <- TRUE
      update_features()
      shiny::showNotification(paste(source_label, "loaded successfully"), type = "message")
    }

    # --- Load path 1: Pre-built object (.RDS / .h5seurat) ---
    shiny::observeEvent(input$loadBtn, {
      shiny::req(input$dataUpload)

      file_path <- input$dataUpload$datapath
      ext <- tolower(tools::file_ext(input$dataUpload$name))

      tryCatch({
        shiny::withProgress(message = "Loading dataset...", {
          if (ext == "h5seurat") {
            if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
              shiny::showNotification("SeuratDisk package required for .h5seurat files", type = "error")
              return()
            }
            seu <- SeuratDisk::LoadH5Seurat(file_path)
            finish_load(SpaCET::convert.Seurat(seu), "Seurat object")
          } else {
            obj <- readRDS(file_path)
            if (inherits(obj, "SpaCET")) {
              finish_load(obj, "SpaCET object")
            } else if (inherits(obj, "Seurat")) {
              finish_load(SpaCET::convert.Seurat(obj), "Seurat object")
            } else {
              shiny::showNotification("File must contain a SpaCET or Seurat object", type = "error")
              return()
            }
          }
        })
      }, error = function(e) {
        shiny::showNotification(paste("Error loading file:", e$message), type = "error")
      })
    })

    # --- Load path 2: Space Ranger zip ---
    shiny::observeEvent(input$loadSpaceRangerBtn, {
      shiny::req(input$spacerUpload)

      tryCatch({
        shiny::withProgress(message = "Loading Space Ranger output...", {
          zip_path <- input$spacerUpload$datapath
          extract_dir <- file.path(tempdir(), "spaceranger_upload")
          if (dir.exists(extract_dir)) unlink(extract_dir, recursive = TRUE)
          on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)

          shiny::incProgress(0.2, detail = "Extracting zip...")
          utils::unzip(zip_path, exdir = extract_dir)

          # Zip structure varies; search for spatial/ subdirectory to locate the data
          spatial_dirs <- list.files(extract_dir, pattern = "^spatial$",
                                     recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
          if (length(spatial_dirs) == 0) {
            shiny::showNotification("No 'spatial/' directory found in zip. Is this a Space Ranger output?", type = "error")
            return()
          }
          visium_path <- dirname(spatial_dirs[1])

          shiny::incProgress(0.4, detail = "Building SpaCET object...")
          obj <- SpaCET::create.SpaCET.object.10X(visiumPath = visium_path)

          shiny::incProgress(0.4, detail = "Done")
          finish_load(obj, "Visium dataset")
        })
      }, error = function(e) {
        shiny::showNotification(paste("Error loading Space Ranger data:", e$message), type = "error")
      })
    })

    # --- Load path: CosMx zip ---
    shiny::observeEvent(input$loadCosmxBtn, {
      shiny::req(input$cosmxUpload)

      tryCatch({
        shiny::withProgress(message = "Loading CosMx data...", {
          zip_path <- input$cosmxUpload$datapath
          extract_dir <- file.path(tempdir(), "cosmx_upload")
          if (dir.exists(extract_dir)) unlink(extract_dir, recursive = TRUE)
          on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)

          shiny::incProgress(0.2, detail = "Extracting zip...")
          utils::unzip(zip_path, exdir = extract_dir)

          # Find the CosMx output directory (look for metadata or expression files)
          meta_files <- list.files(extract_dir, pattern = "metadata", recursive = TRUE, full.names = TRUE)
          if (length(meta_files) == 0) {
            shiny::showNotification("No metadata file found in zip. Is this CosMx output?", type = "error")
            return()
          }
          cosmx_path <- dirname(meta_files[1])

          shiny::incProgress(0.4, detail = "Building SpaCET object...")
          obj <- SpaCET::create.SpaCET.object.CosMx(cosmxPath = cosmx_path)

          shiny::incProgress(0.4, detail = "Done")
          finish_load(obj, "CosMx dataset")
        })
      }, error = function(e) {
        shiny::showNotification(paste("Error loading CosMx data:", e$message), type = "error")
      })
    })

    # --- Load path: Xenium zip ---
    shiny::observeEvent(input$loadXeniumBtn, {
      shiny::req(input$xeniumUpload)

      tryCatch({
        shiny::withProgress(message = "Loading Xenium data...", {
          zip_path <- input$xeniumUpload$datapath
          extract_dir <- file.path(tempdir(), "xenium_upload")
          if (dir.exists(extract_dir)) unlink(extract_dir, recursive = TRUE)
          on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)

          shiny::incProgress(0.2, detail = "Extracting zip...")
          utils::unzip(zip_path, exdir = extract_dir)

          # Find the Xenium output directory (look for cells file)
          cells_files <- list.files(extract_dir, pattern = "cells\\.(csv\\.gz|parquet)$", recursive = TRUE, full.names = TRUE)
          if (length(cells_files) == 0) {
            shiny::showNotification("No cells file found in zip. Is this Xenium output?", type = "error")
            return()
          }
          xenium_path <- dirname(cells_files[1])

          shiny::incProgress(0.4, detail = "Building SpaCET object...")
          obj <- SpaCET::create.SpaCET.object.Xenium(xeniumPath = xenium_path)

          shiny::incProgress(0.4, detail = "Done")
          finish_load(obj, "Xenium dataset")
        })
      }, error = function(e) {
        shiny::showNotification(paste("Error loading Xenium data:", e$message), type = "error")
      })
    })

    # --- Load path 3: Demo dataset ---
    shiny::observeEvent(input$loadDemoBtn, {
      tryCatch({
        shiny::withProgress(message = "Loading Visium HCC demo...", {
          demo_path <- system.file("extdata", "Visium_HCC", package = "SecAct")
          if (demo_path == "") {
            shiny::showNotification("Demo data not found. Is SecAct installed?", type = "error")
            return()
          }

          shiny::incProgress(0.3, detail = "Building SpaCET object...")
          obj <- SpaCET::create.SpaCET.object.10X(visiumPath = demo_path)

          shiny::incProgress(0.6, detail = "Done")
          finish_load(obj, "Visium HCC demo")
        })
      }, error = function(e) {
        shiny::showNotification(paste("Error loading demo:", e$message), type = "error")
      })
    })

    # --- Run SecAct inference ---
    shiny::observeEvent(input$runInferenceBtn, {
      shiny::req(rv$spacet_obj)

      tryCatch({
        shiny::withProgress(message = "Running SecAct inference...", detail = "This may take a few minutes", {
          shiny::incProgress(0.1)
          rv$spacet_obj <- SecAct::SecAct.activity.inference.ST(rv$spacet_obj)
          shiny::incProgress(0.9)
          update_features()

          shiny::showNotification("SecAct inference complete!", type = "message")
        })
      }, error = function(e) {
        shiny::showNotification(paste("Inference error:", e$message), type = "error")
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
        spacet_type <- switch(input$spatialType,
          "SecActActivity" = "GeneExpression",
          "GeneExpression" = "GeneExpression",
          "CellFraction" = "CellFraction",
          "CellTypeComposition" = "CellTypeComposition",
          "GeneExpression"
        )

        if (input$spatialType == "SecActActivity") {
          act_mat <- rv$spacet_obj@results$SecAct_output$SecActTarget
          if (is.null(act_mat)) {
            return(empty_state_plot("No SecAct activity data found. Run inference first."))
          }

          temp_obj <- swap_activity_matrix(rv$spacet_obj, act_mat)

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
