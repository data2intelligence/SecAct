# Run Inference Module
# Users upload expression data (CSV/TSV/RDS) and run SecAct inference.

source(file.path("R", "utils_viz.R"), local = TRUE)

inferenceUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Upload Data", style = paste0("color: ", UI_COLORS$primary)),
          shiny::fileInput(ns("userFile"), "Expression data (CSV, TSV, or RDS):",
                           accept = c(".csv", ".tsv", ".txt", ".rds", ".RDS")),
          shiny::radioButtons(ns("inputType"), "Input type:",
                              choices = c("Differential expression (logFC)" = "logFC",
                                          "Raw expression matrix" = "expression"),
                              selected = "logFC"),
          shiny::actionButton(ns("submitBtn"), "Run SecAct Inference",
                              class = "btn-primary btn-block",
                              icon = shiny::icon("play")),
          shiny::hr(),
          shiny::textOutput(ns("statusText"))
        )
      ),
      shiny::column(8,
        shiny::conditionalPanel(
          condition = paste0("output['", ns("hasResults"), "']"),
          shiny::wellPanel(
            shiny::h4("Results", style = paste0("color: ", UI_COLORS$primary)),
            DT::dataTableOutput(ns("resultsTable")),
            shiny::br(),
            shiny::downloadButton(ns("downloadResults"), "Download Results (CSV)",
                                  class = "btn-success")
          )
        ),
        shiny::conditionalPanel(
          condition = paste0("!output['", ns("hasResults"), "']"),
          shiny::div(
            style = "text-align: center; padding: 80px 0; color: #999;",
            shiny::h3("Upload expression data and run inference"),
            shiny::p("Rows = genes, Columns = samples/spots. First column = gene names.")
          )
        )
      )
    )
  )
}

inferenceServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(results = NULL)

    output$hasResults <- shiny::reactive({ !is.null(rv$results) })
    shiny::outputOptions(output, "hasResults", suspendWhenHidden = FALSE)

    output$statusText <- shiny::renderText({ "" })

    shiny::observeEvent(input$submitBtn, {
      shiny::req(input$userFile)

      output$statusText <- shiny::renderText({ "Running inference..." })

      tryCatch({
        shiny::withProgress(message = "Running SecAct inference...", detail = "This may take a few minutes", {
          file_path <- input$userFile$datapath
          ext <- tolower(tools::file_ext(input$userFile$name))

          shiny::incProgress(0.1, detail = "Loading data...")

          # Load expression data
          if (ext %in% c("rds")) {
            expr_data <- readRDS(file_path)
          } else if (ext %in% c("csv")) {
            expr_data <- read.csv(file_path, row.names = 1, check.names = FALSE)
          } else {
            expr_data <- read.delim(file_path, row.names = 1, check.names = FALSE)
          }

          shiny::incProgress(0.3, detail = "Running inference...")

          # Run SecAct inference
          is_diff <- input$inputType == "logFC"
          result <- SecAct::SecAct.activity.inference(
            inputProfile = as.matrix(expr_data),
            is.differential = is_diff
          )

          shiny::incProgress(0.5, detail = "Processing results...")

          # Extract z-scores
          if (!is.null(result$zscore)) {
            rv$results <- as.data.frame(result$zscore)
          } else if (!is.null(result)) {
            rv$results <- as.data.frame(result)
          }

          shiny::incProgress(0.1)
          output$statusText <- shiny::renderText({ "Inference complete!" })
        })
      }, error = function(e) {
        output$statusText <- shiny::renderText({ paste("Error:", e$message) })
      })
    })

    output$resultsTable <- DT::renderDataTable({
      shiny::req(rv$results)
      DT::datatable(rv$results, options = list(pageLength = 15, scrollX = TRUE))
    })

    output$downloadResults <- shiny::downloadHandler(
      filename = function() {
        paste0("secact_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        shiny::req(rv$results)
        write.csv(rv$results, file)
      }
    )
  })
}
