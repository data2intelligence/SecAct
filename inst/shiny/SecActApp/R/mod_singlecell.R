# Single-Cell Analysis Module
# Upload scRNA-seq data, run SecAct inference, visualize activity by cell type.

source(file.path("R", "utils_viz.R"), local = TRUE)

singlecellUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shiny::column(3,
        shiny::wellPanel(
          shiny::h4("Upload scRNA-seq Data", style = paste0("color: ", UI_COLORS$primary)),
          shiny::fileInput(ns("scFile"), "Seurat object or expression matrix (.RDS)",
                           accept = c(".rds", ".RDS")),
          shiny::textInput(ns("cellTypeCol"), "Cell type column:", value = "cellType"),
          shiny::actionButton(ns("runScBtn"), "Run SC Inference",
                              class = "btn-primary btn-block",
                              icon = shiny::icon("play")),
          shiny::hr(),
          shiny::textOutput(ns("statusText")),

          shiny::conditionalPanel(
            condition = paste0("output['", ns("hasResults"), "']"),
            shiny::hr(),
            shiny::h4("Visualization", style = paste0("color: ", UI_COLORS$primary)),
            shinyWidgets::pickerInput(
              ns("protein"), "Select Protein:",
              choices = NULL,
              options = list(`live-search` = TRUE, size = 10),
              multiple = FALSE
            ),
            shiny::selectInput(ns("plotType"), "Plot Type:",
                               choices = c("Activity Bar" = "bar",
                                           "Activity Heatmap" = "heatmap",
                                           "Activity Lollipop" = "lollipop"))
          )
        )
      ),
      shiny::column(9,
        shiny::conditionalPanel(
          condition = paste0("!output['", ns("hasResults"), "']"),
          shiny::div(
            style = "text-align: center; padding: 80px 0; color: #999;",
            shiny::h3("Single-Cell SecAct Analysis"),
            shiny::p("Upload an scRNA-seq dataset and run cell-type level inference."),
            shiny::p("Provide a Seurat object (.RDS) with cell type annotations.")
          )
        ),
        shiny::conditionalPanel(
          condition = paste0("output['", ns("hasResults"), "']"),
          shiny::fluidRow(
            shiny::column(12,
              shiny::div(
                style = "border: 1px solid #ddd; border-radius: 5px; padding: 10px;",
                shiny::plotOutput(ns("activityPlot"), height = "500px")
              )
            )
          ),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(12,
              DT::dataTableOutput(ns("resultsTable"))
            )
          )
        )
      )
    )
  )
}

singlecellServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(
      sc_result = NULL,
      zscore_mat = NULL
    )

    output$hasResults <- shiny::reactive({ !is.null(rv$zscore_mat) })
    shiny::outputOptions(output, "hasResults", suspendWhenHidden = FALSE)

    output$statusText <- shiny::renderText({ "" })

    shiny::observeEvent(input$runScBtn, {
      shiny::req(input$scFile)
      output$statusText <- shiny::renderText({ "Running SC inference..." })

      tryCatch({
        shiny::withProgress(message = "Running single-cell SecAct inference...",
                            detail = "This may take several minutes", {
          shiny::incProgress(0.1, detail = "Loading data...")
          obj <- readRDS(input$scFile$datapath)

          shiny::incProgress(0.3, detail = "Running inference...")

          # Determine if Seurat object or raw matrix
          if (inherits(obj, "Seurat")) {
            result <- SecAct::SecAct.activity.inference.scRNAseq(obj,
                        cellTypeColumn = input$cellTypeCol)
          } else if (is.matrix(obj) || inherits(obj, "dgCMatrix")) {
            result <- SecAct::SecAct.activity.inference(
                        inputProfile = obj)
          } else {
            shiny::showNotification("Unsupported object type", type = "error")
            return()
          }

          shiny::incProgress(0.5, detail = "Processing results...")

          rv$sc_result <- result
          if (!is.null(result$zscore)) {
            rv$zscore_mat <- as.data.frame(result$zscore)
            proteins <- rownames(rv$zscore_mat)
            shinyWidgets::updatePickerInput(session, "protein",
                                            choices = proteins,
                                            selected = proteins[1])
          }

          shiny::incProgress(0.1)
          output$statusText <- shiny::renderText({ "Inference complete!" })
        })
      }, error = function(e) {
        output$statusText <- shiny::renderText({ paste("Error:", e$message) })
      })
    })

    output$activityPlot <- shiny::renderPlot({
      shiny::req(rv$zscore_mat, input$protein, input$plotType)

      tryCatch({
        zscore_vec <- as.numeric(rv$zscore_mat[input$protein, ])
        names(zscore_vec) <- colnames(rv$zscore_mat)

        if (input$plotType == "bar") {
          SecAct::SecAct.bar.plot(zscore_vec, title = input$protein)
        } else if (input$plotType == "heatmap") {
          SecAct::SecAct.heatmap.plot(rv$zscore_mat, title = "SecAct Activity")
        } else if (input$plotType == "lollipop") {
          SecAct::SecAct.lollipop.plot(zscore_vec, title = input$protein)
        }
      }, error = function(e) {
        empty_state_plot(paste("Plot error:", e$message))
      })
    })

    output$resultsTable <- DT::renderDataTable({
      shiny::req(rv$zscore_mat)
      DT::datatable(rv$zscore_mat,
                    options = list(pageLength = 15, scrollX = TRUE))
    })
  })
}
