# SecAct Unified Visualization App
# Entry point — loads modules and launches the tabbed interface.
# Run via SecAct::runSecActApp() or shiny::runApp("inst/shiny/SecActApp")

# Load global config and dependencies
source("global.R", local = TRUE)

# Load modules
source(file.path("R", "mod_spatial.R"), local = TRUE)

# --- UI ---
ui <- shiny::fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),
  if (requireNamespace("shinyjs", quietly = TRUE)) shinyjs::useShinyjs(),

  # Header
  shiny::div(
    style = paste0("background-color: ", UI_COLORS$primary, "; color: white; padding: 15px; margin-bottom: 20px;"),
    shiny::h2("SecAct", style = "margin-top: 0; display: inline;"),
    shiny::span("Secreted Protein Activity Analysis", style = "margin-left: 15px; opacity: 0.8;")
  ),

  # Tabbed interface — Phase 1 has spatial only; more tabs added in later phases
  shiny::tabsetPanel(
    id = "mainTabs",
    shiny::tabPanel("Spatial", spatialUI("spatial"))
  )
)

# --- Server ---
server <- function(input, output, session) {
  spatialServer("spatial")
}

shiny::shinyApp(ui = ui, server = server)
