library(shiny)
library(bs4Dash)
library(SeuratVisPro)
library(Seurat)
library(ggplot2)
library(patchwork)
library(DT)

# Create example object (aligned with README examples)
obj <- SeuratVisProExample(
  n_cells = 300,
  n_genes = 1000,
  n_clusters = 10,
  seed = 123,
  genes_mt = "^MT-",
  neighbor_dims = 10,
  cluster_res = 0.5,
  umap_dims = 10,
  spatial = TRUE
)
obj$batch <- sample(c('A', 'B'), ncol(obj), replace = TRUE)

ui <- bs4DashPage(
  title = "SeuratVisPro",
  skin = NULL,
  freshTheme = NULL,
  preloader = NULL,
  options = NULL,
  fullscreen = TRUE,
  help = TRUE,
  dark = FALSE,
  scrollToTop = TRUE,
  header = bs4DashNavbar(skin = "light"),
  sidebar = bs4DashSidebar(
    disable = FALSE,
    width = NULL,
    skin = "dark",
    status = "warning",
    elevation = 3,
    collapsed = FALSE,
    minified = TRUE,
    expandOnHover = TRUE,
    fixed = TRUE,
    id = NULL,
    customArea = NULL,
    bs4SidebarMenu(
      id = NULL,
      .list = NULL,
      flat = FALSE,
      compact = FALSE,
      childIndent = FALSE,
      legacy = FALSE,
      bs4SidebarMenuItem(
        text = "QC Panel",
        tabName = "qc",
        icon = icon("gear"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Cluster Stability",
        tabName = "stab",
        icon = icon("sliders-h"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Marker Atlas",
        tabName = "markers",
        icon = icon("th"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Batch Mixing",
        tabName = "batch",
        icon = icon("object-group"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Gene Trend",
        tabName = "trend",
        icon = icon("project-diagram"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Lig-Rec",
        tabName = "lr",
        icon = icon("link"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Meta Feature",
        tabName = "ms",
        icon = icon("layer-group"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Cell Cycle",
        tabName = "cc",
        icon = icon("spinner"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Embedding Contour",
        tabName = "ec",
        icon = icon("draw-polygon"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Gene Coexp Hive",
        tabName = "ch",
        icon = icon("braille"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Spatial Overlay",
        tabName = "sp",
        icon = icon("map"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Hex Entropy",
        tabName = "hex",
        icon = icon("th"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Local Moran",
        tabName = "lisa",
        icon = icon("dot-circle"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      ),
      bs4SidebarMenuItem(
        text = "Flow Graph",
        tabName = "flow",
        icon = icon("project-diagram"),
        badgeLabel = NULL,
        badgeColor = NULL,
        href = NULL,
        newTab = TRUE,
        selected = NULL,
        expandedName = NULL,
        startExpanded = FALSE,
        condition = NULL
      )
    )
  ),
  body = bs4DashBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    bs4TabItems(
      bs4TabItem(tabName = "qc", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            textInput("mt", label = "MT Regex", value = "^MT-"),
            textInput("ribo", label = "Ribo Regex", value = "^RPL|^RPS"),
            selectInput(
              "group",
              label = "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "QC Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("qc_plot", height = "600px")
          )
        )
      )),

      # ---- Cluster Stability Tab ----
      bs4TabItem(tabName = "stab", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            sliderInput(
              "resMin",
              "Resolution Min",
              min = 0.1,
              max = 2,
              value = 0.2,
              step = 0.1
            ),
            sliderInput(
              "resMax",
              "Resolution Max",
              min = 0.1,
              max = 2,
              value = 1.0,
              step = 0.1
            ),
            sliderInput(
              "resStep",
              "Step",
              min = 0.1,
              max = 0.5,
              value = 0.2,
              step = 0.1
            ),
            numericInput("reps", "Repetitions", value = 3, min = 1)
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Stability Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("stab_plot", height = "400px")
          ),
          bs4Card(
            title = "Summary Table",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            DTOutput("stab_table", height = "200px")
          )
        )
      )),

      # ---- Markers Tab ----
      bs4TabItem(tabName = "markers", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            numericInput("topn", "Top N Markers", value = 5, min = 1)
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Marker Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("marker_plot", height = "400px")
          ),
          bs4Card(
            title = "Marker Table",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            DTOutput("marker_table", height = "200px")
          )
        )
      )),

      # ---- Batch Mixing Tab ----
      bs4TabItem(tabName = "batch", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            selectInput(
              "batch",
              "Batch Column",
              choices = colnames(obj@meta.data),
              selected = "batch"
            )
          )
        ), column(
          width = 8,
          bs4Card(
            title = "Batch Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("batch_plot", height = "400px")
          ),
          bs4Card(
            title = "Batch Summary",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            DTOutput("batch_table", height = "200px")
          )
        )
      )),

      # ---- Gene Trend Tab ----
      bs4TabItem(tabName = "trend", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            textInput("genes", "Genes (comma)", value = "G10,G20,G30"),
            selectInput(
              "trendBy",
              "Trend By",
              choices = c("pseudotime", colnames(obj@meta.data)),
              selected = "pseudotime"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Trend Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("trend_plot", height = "600px")
          )
        )
      )),

      # ---- Lig-Rec Tab ----
      bs4TabItem(tabName = "lr", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            actionButton("lr_demo", "Load Example LR Table"),
            selectInput(
              "lr_group",
              "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Lig-Rec Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("lr_plot", height = "400px")
          ),
          bs4Card(
            title = "LR Scores",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            DTOutput("lr_table", height = "200px")
          )
        )
      )),

      # ---- Module Score Tab ----
      bs4TabItem(tabName = "ms", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            textInput("setA", "Set A Genes", value = paste0("G", 1:10)),
            textInput("setB", "Set B Genes", value = paste0("G", 11:20)),
            selectInput(
              "ms_group",
              "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Module Score Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("ms_plot", height = "600px")
          )
        )
      )),

      # ---- Cell Cycle Tab ----
      bs4TabItem(tabName = "cc", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            textInput("sGenes", "S-phase Genes", value = paste0("G", 1:10)),
            textInput("g2mGenes", "G2M-phase Genes", value = paste0("G", 11:20))
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Cell Cycle Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("cc_plot", height = "600px")
          )
        )
      )),

      # ---- Embedding Contour Tab ----
      bs4TabItem(tabName = "ec", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            selectInput(
              "ec_group",
              "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Contour Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("ec_plot", height = "600px")
          )
        )
      )),



      # ---- Coexp Hive Tab ----
      bs4TabItem(tabName = "ch", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            textInput("ch_genes", "Genes (comma)", value = paste0("G", 1:12)),
            numericInput(
              "ch_thr",
              "Correlation Threshold",
              value = 0.2,
              min = 0,
              max = 1,
              step = 0.05
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Co-expression Hive",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("ch_plot", height = "600px")
          )
        )
      )),

      # ---- Spatial Overlay Tab ----
      bs4TabItem(tabName = "sp", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            textInput("sp_feats", "Features (comma)", value = "G1,G2,G3")
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Spatial Overlay",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("sp_plot", height = "600px")
          )
        )
      ))
      ,
      bs4TabItem(tabName = "hex", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            selectInput(
              "hex_group",
              "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            ),
            sliderInput(
              "hex_bins",
              "Bins",
              min = 10,
              max = 60,
              value = 30,
              step = 5
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Hex Entropy",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("hex_plot", height = "600px")
          )
        )
      ))
      ,
      bs4TabItem(tabName = "lisa", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            textInput("lisa_gene", "Gene", value = "G10"),
            sliderInput(
              "lisa_k",
              "Neighbors (k)",
              min = 5,
              max = 50,
              value = 15,
              step = 5
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Local Moran's I",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("lisa_plot", height = "600px")
          )
        )
      ))
      ,
      bs4TabItem(tabName = "flow", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            selectInput(
              "flow_group",
              "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Cluster Flow Graph",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            height = NULL,
            plotOutput("flow_plot", height = "600px")
          )
        )
      ))
    )
  )
)

# -----------------------------
# Server logic (unchanged)
# -----------------------------
server <- function(input, output, session) {
  output$qc_plot <- renderPlot({
    p <- VisQCPanel(
      obj,
      genes_mt = input$mt,
      genes_ribo = input$ribo,
      group.by = input$group,
      interactive = FALSE
    )
    print(p)
  })

  output$stab_plot <- renderPlot({
    res <- seq(input$resMin, input$resMax, by = input$resStep)
    r <- VisClusterStability(
      obj,
      resolution_range = res,
      dims = 1:10,
      reps = input$reps,
      prop = 0.8,
      palette = "C"
    )
    print(r$plot)
  })
  output$stab_table <- renderDT({
    res <- seq(input$resMin, input$resMax, by = input$resStep)
    r <- VisClusterStability(
      obj,
      resolution_range = res,
      dims = 1:10,
      reps = input$reps,
      prop = 0.8,
      palette = "C"
    )
    datatable(r$summary, options = list(scrollY = "150px", pageLength = 10))
  })

  output$marker_plot <- renderPlot({
    r <- VisMarkerAtlas(obj, markers_top = input$topn)
    print(r$plot)
  })
  output$marker_table <- renderDT({
    r <- VisMarkerAtlas(obj, markers_top = input$topn)
    datatable(r$markers, options = list(scrollY = "150px", pageLength = 10))
  })

  output$batch_plot <- renderPlot({
    r <- VisBatchAlign(
      obj,
      batch = input$batch,
      reduction = 'pca',
      dims = 1:10,
      k = 20,
      palette = "C"
    )
    print(r$plot)
  })
  output$batch_table <- renderDT({
    r <- VisBatchAlign(
      obj,
      batch = input$batch,
      reduction = 'pca',
      dims = 1:10,
      k = 20,
      palette = "C"
    )
    datatable(r$summary, options = list(scrollY = "150px", pageLength = 10))
  })

  output$trend_plot <- renderPlot({
    feats <- trimws(unlist(strsplit(input$genes, ",")))
    p <- VisGeneTrend(
      obj,
      features = feats,
      by = input$trendBy,
      reduction = 'umap',
      dims = 1:2,
      smooth.method = 'loess',
      palette = 'C'
    )
    print(p)
  })

  lr_demo <- reactiveVal(NULL)
  observeEvent(input$lr_demo, {
    lr_demo(data.frame(
      ligand = paste0('G', 1:5),
      receptor = paste0('G', 6:10)
    ))
  })
  output$lr_plot <- renderPlot({
    lr <- lr_demo()
    req(lr)
    r <- VisLigRec(
      obj,
      lr_table = lr,
      group.by = input$lr_group,
      palette = 'C',
      tile_alpha = 0.8
    )
    print(r$plot)
  })
  output$lr_table <- renderDT({
    lr <- lr_demo()
    req(lr)
    r <- VisLigRec(obj, lr_table = lr, group.by = input$lr_group)
    datatable(r$scores, options = list(scrollY = "150px", pageLength = 10))
  })

  output$ms_plot <- renderPlot({
    setA <- trimws(unlist(strsplit(input$setA, ",")))
    setB <- trimws(unlist(strsplit(input$setB, ",")))
    r <- VisMetaFeature(
      obj,
      feature_sets = list(SetA = setA, SetB = setB),
      group.by = input$ms_group
    )
    print(r$plot)
  })

  output$cc_plot <- renderPlot({
    s.genes <- trimws(unlist(strsplit(input$sGenes, ",")))
    g2m.genes <- trimws(unlist(strsplit(input$g2mGenes, ",")))
    r <- VisCellCycle(obj, genes_s = s.genes, genes_g2m = g2m.genes)
    print(r$plot)
  })

  output$ec_plot <- renderPlot({
    print(
      VisEmbeddingContour(
        obj,
        group.by = input$ec_group,
        reduction = 'umap',
        levels = 5,
        palette = 'C'
      )
    )
  })



  output$ch_plot <- renderPlot({
    genes <- trimws(unlist(strsplit(input$ch_genes, ",")))
    print(
      VisGeneCoexpHive(
        obj,
        genes = genes,
        reduction = 'pca',
        threshold = input$ch_thr,
        palette = 'C'
      )
    )
  })

  output$sp_plot <- renderPlot({
    feats <- trimws(unlist(strsplit(input$sp_feats, ",")))
    print(
      VisSpatialOverlay(
        obj,
        features = feats,
        image = NULL,
        coords_cols = c('x', 'y'),
        palette = 'C',
        point_size = 2,
        alpha = 0.5
      )
    )
  })
  output$hex_plot <- renderPlot({
    print(
      VisHexEntropy(
        obj,
        group.by = input$hex_group,
        reduction = 'umap',
        bins = input$hex_bins,
        palette = 'C'
      )
    )
  })
  output$lisa_plot <- renderPlot({
    print(
      VisLocalMoran(
        obj,
        gene = input$lisa_gene,
        reduction = 'umap',
        k = input$lisa_k,
        palette = 'C',
        point_size = 2,
        point_alpha = 0.8
      )
    )
  })
  output$flow_plot <- renderPlot({
    print(
      VisClusterFlowGraph(
        obj,
        group.by = input$flow_group,
        reduction = 'umap',
        palette = 'C',
        point_size = 7,
        point_alpha = 0.9,
        label_size = 5
      )
    )
  })
}

# Run app
shinyApp(ui, server)
