#' @title Launch SeuratVisPro Shiny app
#' @description Launch the bundled Shiny application showcasing key functions and examples.
#' @param port Port to run the app on.
#' @param host Host interface.
#' @return No return value; launches a local Shiny app.
#' @examples
#' \dontrun{
#' launchSeuratVisPro()
#' }
#' @export
launchSeuratVisPro <- function(port = 8888, host = "127.0.0.1") {
  if (requireNamespace("bs4Dash", quietly = TRUE)) {
    app_dir <- system.file("shiny", package = "SeuratVisPro")
    if (app_dir == "" || !dir.exists(app_dir)) {
      app_dir <- normalizePath(file.path("inst", "shiny"), mustWork = FALSE)
    }
    if (!dir.exists(app_dir)) stop("Shiny app directory not found")
    return(shiny::runApp(appDir = app_dir, port = port, host = host, launch.browser = FALSE))
  }
  obj <- SeuratVisProExample(n_cells = 300, n_genes = 800, n_clusters = 4, spatial = TRUE)
  obj$batch <- sample(c('A','B'), ncol(obj), replace = TRUE)
  ui <- shiny::fluidPage(
    shiny::titlePanel("SeuratVisPro (fallback UI)"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("group", label = "Group by", choices = colnames(obj@meta.data), selected = "seurat_clusters"),
        shiny::textInput("features", "Features (comma)", value = "G1,G2,G3"),
        shiny::sliderInput("bins", "Hex bins", min = 10, max = 60, value = 30, step = 5)
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(id = "tabs",
          shiny::tabPanel("QC", shiny::plotOutput("qc", height = "500px")),
          shiny::tabPanel("Markers", shiny::plotOutput("markers", height = "500px")),
          shiny::tabPanel("Hex Entropy", shiny::plotOutput("hex", height = "500px")),
          shiny::tabPanel("Spatial", shiny::plotOutput("sp", height = "500px"))
        )
      )
    )
  )
  server <- function(input, output, session) {
    output$qc <- shiny::renderPlot({
      print(VisQCPanel(obj, group.by = input$group))
    })
    output$markers <- shiny::renderPlot({
      print(VisMarkerAtlas(obj)$plot)
    })
    output$hex <- shiny::renderPlot({
      print(VisHexEntropy(obj, group.by = input$group, bins = input$bins))
    })
    output$sp <- shiny::renderPlot({
      feats <- trimws(unlist(strsplit(input$features, ",")))
      print(VisSpatialOverlay(obj, features = feats))
    })
  }
  shiny::runApp(shiny::shinyApp(ui, server), port = port, host = host, launch.browser = FALSE)
}
