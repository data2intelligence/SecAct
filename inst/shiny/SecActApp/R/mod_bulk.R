# Bulk Analysis Module
# Two workflows: Activity Change (treatment vs control) and Cohort Survival.

source(file.path("R", "utils_viz.R"), local = TRUE)

bulkUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::tabsetPanel(
      id = ns("bulkSubTabs"),

      # Sub-tab 1: Activity Change
      shiny::tabPanel("Activity Change",
        shiny::fluidRow(
          shiny::column(4,
            shiny::wellPanel(
              shiny::h4("Activity Change", style = paste0("color: ", UI_COLORS$primary)),
              shiny::p("Compare secreted protein activity between treatment and control.",
                       style = "color: #666; font-size: 0.9em;"),
              shiny::fileInput(ns("treatmentFile"), "Treatment expression (CSV/TSV/RDS):"),
              shiny::fileInput(ns("controlFile"), "Control expression (CSV/TSV/RDS):"),
              shiny::checkboxInput(ns("singleSample"), "Single-sample level analysis", value = FALSE),
              shiny::actionButton(ns("runChangeBtn"), "Run Activity Change",
                                  class = "btn-primary btn-block", icon = shiny::icon("play")),
              shiny::hr(),
              shiny::textOutput(ns("changeStatus"))
            )
          ),
          shiny::column(8,
            shiny::conditionalPanel(
              condition = paste0("output['", ns("hasChangeResults"), "']"),
              shiny::wellPanel(
                shiny::h4("Activity Change Results", style = paste0("color: ", UI_COLORS$primary)),
                shiny::plotOutput(ns("barPlot"), height = "400px"),
                shiny::br(),
                DT::dataTableOutput(ns("changeTable")),
                shiny::br(),
                shiny::downloadButton(ns("downloadChange"), "Download Results (CSV)", class = "btn-success")
              )
            ),
            shiny::conditionalPanel(
              condition = paste0("!output['", ns("hasChangeResults"), "']"),
              shiny::div(
                style = "text-align: center; padding: 80px 0; color: #999;",
                shiny::h3("Upload treatment and control expression data"),
                shiny::p("Rows = genes, Columns = samples. Values should be log2(x+1) transformed.")
              )
            )
          )
        )
      ),

      # Sub-tab 2: Cohort Survival
      shiny::tabPanel("Cohort Survival",
        shiny::fluidRow(
          shiny::column(4,
            shiny::wellPanel(
              shiny::h4("Cohort Survival", style = paste0("color: ", UI_COLORS$primary)),
              shiny::p("Link secreted protein activity to clinical outcomes.",
                       style = "color: #666; font-size: 0.9em;"),
              shiny::fileInput(ns("cohortExprFile"), "Expression matrix (CSV/TSV/RDS):"),
              shiny::fileInput(ns("clinicalFile"), "Clinical data (CSV/TSV):"),
              shiny::p("Clinical file must have 'Time' and 'Event' columns.",
                       style = "color: #999; font-size: 0.85em;"),
              shiny::actionButton(ns("runSurvivalBtn"), "Run Survival Analysis",
                                  class = "btn-primary btn-block", icon = shiny::icon("play")),
              shiny::hr(),
              shiny::textOutput(ns("survivalStatus")),

              shiny::conditionalPanel(
                condition = paste0("output['", ns("hasSurvivalResults"), "']"),
                shiny::hr(),
                shinyWidgets::pickerInput(
                  ns("survivalProtein"), "Select Protein for Survival Curve:",
                  choices = NULL,
                  options = list(`live-search` = TRUE, size = 10),
                  multiple = FALSE
                )
              )
            )
          ),
          shiny::column(8,
            shiny::conditionalPanel(
              condition = paste0("output['", ns("hasSurvivalResults"), "']"),
              shiny::wellPanel(
                shiny::h4("Risk Scores", style = paste0("color: ", UI_COLORS$primary)),
                shiny::plotOutput(ns("lollipopPlot"), height = "400px"),
                shiny::hr(),
                shiny::h4("Survival Curve", style = paste0("color: ", UI_COLORS$primary)),
                shiny::plotOutput(ns("survivalPlot"), height = "400px"),
                shiny::br(),
                DT::dataTableOutput(ns("survivalTable")),
                shiny::br(),
                shiny::downloadButton(ns("downloadSurvival"), "Download Risk Scores (CSV)", class = "btn-success")
              )
            ),
            shiny::conditionalPanel(
              condition = paste0("!output['", ns("hasSurvivalResults"), "']"),
              shiny::div(
                style = "text-align: center; padding: 80px 0; color: #999;",
                shiny::h3("Upload expression and clinical data"),
                shiny::p("Expression: rows = genes, columns = patients."),
                shiny::p("Clinical: must have 'Time' (survival time) and 'Event' (0/1) columns.")
              )
            )
          )
        )
      )
    )
  )
}

bulkServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    rv <- shiny::reactiveValues(
      change_results = NULL,
      activity_mat = NULL,
      clinical_data = NULL,
      risk_mat = NULL
    )

    output$hasChangeResults <- shiny::reactive({ !is.null(rv$change_results) })
    shiny::outputOptions(output, "hasChangeResults", suspendWhenHidden = FALSE)

    output$hasSurvivalResults <- shiny::reactive({ !is.null(rv$risk_mat) })
    shiny::outputOptions(output, "hasSurvivalResults", suspendWhenHidden = FALSE)

    output$changeStatus <- shiny::renderText({ "" })
    output$survivalStatus <- shiny::renderText({ "" })

    # Helper to load expression data from uploaded file
    load_expr <- function(file_input) {
      ext <- tolower(tools::file_ext(file_input$name))
      if (ext == "rds") return(readRDS(file_input$datapath))
      if (ext == "csv") return(as.matrix(read.csv(file_input$datapath, row.names = 1, check.names = FALSE)))
      as.matrix(read.delim(file_input$datapath, row.names = 1, check.names = FALSE))
    }

    # --- Activity Change workflow ---
    shiny::observeEvent(input$runChangeBtn, {
      shiny::req(input$treatmentFile, input$controlFile)
      output$changeStatus <- shiny::renderText({ "Running..." })

      tryCatch({
        shiny::withProgress(message = "Computing activity change...", {
          shiny::incProgress(0.1, detail = "Loading data...")
          treatment <- load_expr(input$treatmentFile)
          control <- load_expr(input$controlFile)

          shiny::incProgress(0.3, detail = "Running inference...")
          res <- SecAct::SecAct.activity.inference(
            inputProfile = treatment,
            inputProfile_control = control,
            is.singleSampleLevel = input$singleSample
          )

          shiny::incProgress(0.5, detail = "Done")
          rv$change_results <- as.data.frame(res$zscore)
          output$changeStatus <- shiny::renderText({ "Complete!" })
        })
      }, error = function(e) {
        output$changeStatus <- shiny::renderText({ paste("Error:", e$message) })
      })
    })

    output$barPlot <- shiny::renderPlot({
      shiny::req(rv$change_results)
      tryCatch({
        if (ncol(rv$change_results) == 1) {
          fg_vec <- rv$change_results[, 1]
          names(fg_vec) <- rownames(rv$change_results)
          SecAct::SecAct.bar.plot(fg_vec, title = "Activity Change")
        } else {
          SecAct::SecAct.heatmap.plot(as.matrix(rv$change_results), title = "Activity Change")
        }
      }, error = function(e) { empty_state_plot(paste("Plot error:", e$message)) })
    })

    output$changeTable <- DT::renderDataTable({
      shiny::req(rv$change_results)
      DT::datatable(rv$change_results, options = list(pageLength = 15, scrollX = TRUE))
    })

    output$downloadChange <- shiny::downloadHandler(
      filename = function() paste0("activity_change_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) { shiny::req(rv$change_results); write.csv(rv$change_results, file) }
    )

    # --- Cohort Survival workflow ---
    shiny::observeEvent(input$runSurvivalBtn, {
      shiny::req(input$cohortExprFile, input$clinicalFile)
      output$survivalStatus <- shiny::renderText({ "Running..." })

      tryCatch({
        shiny::withProgress(message = "Running survival analysis...", {
          shiny::incProgress(0.1, detail = "Loading expression data...")
          expr <- load_expr(input$cohortExprFile)

          shiny::incProgress(0.1, detail = "Loading clinical data...")
          ext <- tolower(tools::file_ext(input$clinicalFile$name))
          if (ext == "csv") {
            clinical <- read.csv(input$clinicalFile$datapath, row.names = 1, check.names = FALSE)
          } else {
            clinical <- read.delim(input$clinicalFile$datapath, row.names = 1, check.names = FALSE)
          }
          rv$clinical_data <- clinical

          shiny::incProgress(0.3, detail = "Running SecAct inference...")
          res <- SecAct::SecAct.activity.inference(inputProfile = expr, inputProfile_control = NULL)
          rv$activity_mat <- res$zscore

          shiny::incProgress(0.3, detail = "Computing risk scores...")
          rv$risk_mat <- SecAct::SecAct.coxph.regression(rv$activity_mat, clinical)

          proteins <- rownames(rv$risk_mat)
          shinyWidgets::updatePickerInput(session, "survivalProtein",
                                          choices = proteins, selected = proteins[1])

          shiny::incProgress(0.2, detail = "Done")
          output$survivalStatus <- shiny::renderText({ "Complete!" })
        })
      }, error = function(e) {
        output$survivalStatus <- shiny::renderText({ paste("Error:", e$message) })
      })
    })

    output$lollipopPlot <- shiny::renderPlot({
      shiny::req(rv$risk_mat)
      tryCatch({
        risk_vec <- rv$risk_mat[, "zscore"]
        names(risk_vec) <- rownames(rv$risk_mat)
        SecAct::SecAct.lollipop.plot(risk_vec, title = "Risk Score")
      }, error = function(e) { empty_state_plot(paste("Plot error:", e$message)) })
    })

    output$survivalPlot <- shiny::renderPlot({
      shiny::req(rv$activity_mat, rv$clinical_data, input$survivalProtein)
      tryCatch({
        SecAct::SecAct.survival.plot(rv$activity_mat, rv$clinical_data,
                                     input$survivalProtein, x.title = "Time")
      }, error = function(e) { empty_state_plot(paste("Plot error:", e$message)) })
    })

    output$survivalTable <- DT::renderDataTable({
      shiny::req(rv$risk_mat)
      DT::datatable(as.data.frame(rv$risk_mat), options = list(pageLength = 15, scrollX = TRUE))
    })

    output$downloadSurvival <- shiny::downloadHandler(
      filename = function() paste0("risk_scores_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) { shiny::req(rv$risk_mat); write.csv(rv$risk_mat, file) }
    )
  })
}
