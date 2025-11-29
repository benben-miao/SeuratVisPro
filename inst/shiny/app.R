library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)
library(DT)

# Load SeuratVisPro
if (requireNamespace("SeuratVisPro", quietly = TRUE)) {
  library(SeuratVisPro)
} else {
  rdir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  if (dir.exists(rdir)) {
    rfiles <- list.files(rdir, pattern = "\\.R$", full.names = TRUE)
    for (f in rfiles) source(f)
  }
}

# Create example object
obj <- SeuratVisProExample(n_cells = 400, n_genes = 800, n_clusters = 4, spatial = TRUE)
obj$batch <- sample(c('A','B'), ncol(obj), replace = TRUE)

# Check bs4Dash
if (!requireNamespace("bs4Dash", quietly = TRUE)) {
  stop("Please install 'bs4Dash': install.packages('bs4Dash')")
}
library(bs4Dash)

# -----------------------------
# UI: Optimized layout
# -----------------------------
ui <- bs4DashPage(
  title = "SeuratVisPro",
  # ✅ FIX: navbar -> header
  header = bs4DashNavbar(skin = "light"),
  sidebar = bs4DashSidebar(
    skin = "light",
    title = "SeuratVisPro",
    bs4SidebarMenu(
      bs4SidebarMenuItem("QC", tabName = "qc", icon = icon("chart-area")),
      bs4SidebarMenuItem("Cluster Stability", tabName = "stab", icon = icon("sliders-h")),
      bs4SidebarMenuItem("Markers", tabName = "markers", icon = icon("th")),
      bs4SidebarMenuItem("Batch Mixing", tabName = "batch", icon = icon("object-group")),
      bs4SidebarMenuItem("Gene Trend", tabName = "trend", icon = icon("project-diagram")),
      bs4SidebarMenuItem("Lig-Rec", tabName = "lr", icon = icon("link")),
      bs4SidebarMenuItem("Module Score", tabName = "ms", icon = icon("layer-group")),
      bs4SidebarMenuItem("Cell Cycle", tabName = "cc", icon = icon("spinner")),
      bs4SidebarMenuItem("Embedding Contour", tabName = "ec", icon = icon("draw-polygon")),
      bs4SidebarMenuItem("Ranked Dotmap", tabName = "rd", icon = icon("dot-circle")),
      bs4SidebarMenuItem("Coexp Hive", tabName = "ch", icon = icon("braille")),
      bs4SidebarMenuItem("Spatial Overlay", tabName = "sp", icon = icon("map"))
      ,bs4SidebarMenuItem("Hex Entropy", tabName = "hex", icon = icon("th"))
      ,bs4SidebarMenuItem("Local Moran", tabName = "lisa", icon = icon("dot-circle"))
      ,bs4SidebarMenuItem("Flow Graph", tabName = "flow", icon = icon("project-diagram"))
    )
  ),
  body = bs4DashBody(
    bs4TabItems(
      # ---- QC Tab ----
      bs4TabItem(tabName = "qc",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  textInput("mt", label = "MT Regex", value = "^MT-"),
                                  textInput("ribo", label = "Ribo Regex", value = "^RPL|^RPS"),
                                  selectInput("group", label = "Group by", choices = colnames(obj@meta.data), selected = "seurat_clusters")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "QC Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("qc_plot", height = "600px")
                          )
                   )
                 )
      ),

      # ---- Cluster Stability Tab ----
      bs4TabItem(tabName = "stab",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  sliderInput("resMin", "Resolution Min", min = 0.1, max = 2, value = 0.2, step = 0.1),
                                  sliderInput("resMax", "Resolution Max", min = 0.1, max = 2, value = 1.0, step = 0.1),
                                  sliderInput("resStep", "Step", min = 0.1, max = 0.5, value = 0.2, step = 0.1),
                                  numericInput("reps", "Repetitions", value = 3, min = 1)
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Stability Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("stab_plot", height = "400px")
                          ),
                          bs4Card(title = "Summary Table", status = "warning", solidHeader = TRUE,
                                  DTOutput("stab_table", height = "200px")
                          )
                   )
                 )
      ),

      # ---- Markers Tab ----
      bs4TabItem(tabName = "markers",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  numericInput("topn", "Top N Markers", value = 5, min = 1)
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Marker Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("marker_plot", height = "400px")
                          ),
                          bs4Card(title = "Marker Table", status = "warning", solidHeader = TRUE,
                                  DTOutput("marker_table", height = "200px")
                          )
                   )
                 )
      ),

      # ---- Batch Mixing Tab ----
      bs4TabItem(tabName = "batch",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  selectInput("batch", "Batch Column", choices = colnames(obj@meta.data), selected = "batch")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Batch Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("batch_plot", height = "400px")
                          ),
                          bs4Card(title = "Batch Summary", status = "warning", solidHeader = TRUE,
                                  DTOutput("batch_table", height = "200px")
                          )
                   )
                 )
      ),

      # ---- Gene Trend Tab ----
      bs4TabItem(tabName = "trend",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  textInput("genes", "Genes (comma)", value = "G10,G20,G30"),
                                  selectInput("trendBy", "Trend By", choices = c("pseudotime", colnames(obj@meta.data)), selected = "pseudotime")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Trend Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("trend_plot", height = "600px")
                          )
                   )
                 )
      ),

      # ---- Lig-Rec Tab ----
      bs4TabItem(tabName = "lr",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  actionButton("lr_demo", "Load Example LR Table"),
                                  selectInput("lr_group", "Group by", choices = colnames(obj@meta.data), selected = "seurat_clusters")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Lig-Rec Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("lr_plot", height = "400px")
                          ),
                          bs4Card(title = "LR Scores", status = "warning", solidHeader = TRUE,
                                  DTOutput("lr_table", height = "200px")
                          )
                   )
                 )
      ),

      # ---- Module Score Tab ----
      bs4TabItem(tabName = "ms",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  textInput("setA", "Set A Genes", value = paste0("G", 1:10)),
                                  textInput("setB", "Set B Genes", value = paste0("G", 11:20)),
                                  selectInput("ms_group", "Group by", choices = colnames(obj@meta.data), selected = "seurat_clusters")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Module Score Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("ms_plot", height = "600px")
                          )
                   )
                 )
      ),

      # ---- Cell Cycle Tab ----
      bs4TabItem(tabName = "cc",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  textInput("sGenes", "S-phase Genes", value = paste0("G", 1:10)),
                                  textInput("g2mGenes", "G2M-phase Genes", value = paste0("G", 11:20))
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Cell Cycle Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("cc_plot", height = "600px")
                          )
                   )
                 )
      ),

      # ---- Embedding Contour Tab ----
      bs4TabItem(tabName = "ec",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  selectInput("ec_group", "Group by", choices = colnames(obj@meta.data), selected = "seurat_clusters")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Contour Plot", status = "info", solidHeader = TRUE,
                                  plotOutput("ec_plot", height = "600px")
                          )
                   )
                 )
      ),

      # ---- Ranked Dotmap Tab ----
      bs4TabItem(tabName = "rd",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  numericInput("rd_topn", "Top N per Group", value = 6, min = 1),
                                  selectInput("rd_group", "Group by", choices = colnames(obj@meta.data), selected = "seurat_clusters")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Ranked Dotmap", status = "info", solidHeader = TRUE,
                                  plotOutput("rd_plot", height = "600px")
                          )
                   )
                 )
      ),

      # ---- Coexp Hive Tab ----
      bs4TabItem(tabName = "ch",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  textInput("ch_genes", "Genes (comma)", value = paste0("G", 1:12)),
                                  numericInput("ch_thr", "Correlation Threshold", value = 0.2, min = 0, max = 1, step = 0.05)
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Co-expression Hive", status = "info", solidHeader = TRUE,
                                  plotOutput("ch_plot", height = "600px")
                          )
                   )
                 )
      ),

      # ---- Spatial Overlay Tab ----
      bs4TabItem(tabName = "sp",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  textInput("sp_feats", "Features (comma)", value = "G1,G2,G3")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Spatial Overlay", status = "info", solidHeader = TRUE,
                                  plotOutput("sp_plot", height = "600px")
                          )
                   )
                 )
      )
      ,bs4TabItem(tabName = "hex",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  selectInput("hex_group", "Group by", choices = colnames(obj@meta.data), selected = "seurat_clusters"),
                                  sliderInput("hex_bins", "Bins", min = 10, max = 60, value = 30, step = 5)
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Hex Entropy", status = "info", solidHeader = TRUE,
                                  plotOutput("hex_plot", height = "600px")
                          )
                   )
                 )
      )
      ,bs4TabItem(tabName = "lisa",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  textInput("lisa_gene", "Gene", value = "G10"),
                                  sliderInput("lisa_k", "Neighbors (k)", min = 5, max = 50, value = 15, step = 5)
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Local Moran's I", status = "info", solidHeader = TRUE,
                                  plotOutput("lisa_plot", height = "600px")
                          )
                   )
                 )
      )
      ,bs4TabItem(tabName = "flow",
                 fluidRow(
                   column(width = 4,
                          bs4Card(title = "Parameters", status = "primary", solidHeader = TRUE,
                                  selectInput("flow_group", "Group by", choices = colnames(obj@meta.data), selected = "seurat_clusters")
                          )
                   ),
                   column(width = 8,
                          bs4Card(title = "Cluster Flow Graph", status = "info", solidHeader = TRUE,
                                  plotOutput("flow_plot", height = "600px")
                          )
                   )
                 )
      )
    )
  )
)

# -----------------------------
# Server logic (unchanged)
# -----------------------------
server <- function(input, output, session) {
  output$qc_plot <- renderPlot({
    p <- VisQCPanel(obj, genes_mt = input$mt, genes_ribo = input$ribo, group.by = input$group, interactive = FALSE)
    print(p)
  })

  output$stab_plot <- renderPlot({
    res <- seq(input$resMin, input$resMax, by = input$resStep)
    r <- VisClusterStability(obj, resolution_range = res, reps = input$reps)
    print(r$plot)
  })
  output$stab_table <- renderDT({
    res <- seq(input$resMin, input$resMax, by = input$resStep)
    r <- VisClusterStability(obj, resolution_range = res, reps = input$reps)
    datatable(r$summary, options = list(scrollY = "150px", pageLength = 5))
  })

  output$marker_plot <- renderPlot({
    r <- VisMarkerAtlas(obj, top_n = input$topn)
    print(r$plot)
  })
  output$marker_table <- renderDT({
    r <- VisMarkerAtlas(obj, top_n = input$topn)
    datatable(r$markers, options = list(scrollY = "150px", pageLength = 5))
  })

  output$batch_plot <- renderPlot({
    r <- VisBatchAlign(obj, batch = input$batch)
    print(r$plot)
  })
  output$batch_table <- renderDT({
    r <- VisBatchAlign(obj, batch = input$batch)
    datatable(r$summary, options = list(scrollY = "150px", pageLength = 5))
  })

  output$trend_plot <- renderPlot({
    feats <- trimws(unlist(strsplit(input$genes, ",")))
    p <- VisGeneTrend(obj, features = feats, by = input$trendBy)
    print(p)
  })

  lr_demo <- reactiveVal(NULL)
  observeEvent(input$lr_demo, {
    lr_demo(data.frame(ligand = paste0('G', 1:5), receptor = paste0('G', 6:10)))
  })
  output$lr_plot <- renderPlot({
    lr <- lr_demo()
    req(lr)
    r <- VisLigRec(obj, lr_table = lr, group.by = input$lr_group)
    print(r$plot)
  })
  output$lr_table <- renderDT({
    lr <- lr_demo()
    req(lr)
    r <- VisLigRec(obj, lr_table = lr, group.by = input$lr_group)
    datatable(r$scores, options = list(scrollY = "150px", pageLength = 5))
  })

  output$ms_plot <- renderPlot({
    setA <- trimws(unlist(strsplit(input$setA, ",")))
    setB <- trimws(unlist(strsplit(input$setB, ",")))
    r <- VisMetaFeature(obj, feature_sets = list(SetA = setA, SetB = setB), group.by = input$ms_group)
    print(r$plot)
  })

  output$cc_plot <- renderPlot({
    s.genes <- trimws(unlist(strsplit(input$sGenes, ",")))
    g2m.genes <- trimws(unlist(strsplit(input$g2mGenes, ",")))
    r <- VisCellCycle(obj, s.genes = s.genes, g2m.genes = g2m.genes)
    print(r$plot)
  })

  output$ec_plot <- renderPlot({
    print(VisEmbeddingContour(obj, group.by = input$ec_group))
  })

  output$rd_plot <- renderPlot({
    print(VisRankedDotmap(obj, group.by = input$rd_group, top_n = input$rd_topn))
  })

  output$ch_plot <- renderPlot({
    genes <- trimws(unlist(strsplit(input$ch_genes, ",")))
    print(VisGeneCoexpHive(obj, genes = genes, threshold = input$ch_thr))
  })

  output$sp_plot <- renderPlot({
    feats <- trimws(unlist(strsplit(input$sp_feats, ",")))
    print(VisSpatialOverlay(obj, features = feats))
  })
  output$hex_plot <- renderPlot({
    print(VisHexEntropy(obj, group.by = input$hex_group, bins = input$hex_bins))
  })
  output$lisa_plot <- renderPlot({
    print(VisLocalMoran(obj, gene = input$lisa_gene, k = input$lisa_k))
  })
  output$flow_plot <- renderPlot({
    print(VisClusterFlowGraph(obj, group.by = input$flow_group))
  })
}

# Run app
shinyApp(ui, server)
