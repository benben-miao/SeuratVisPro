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
              "qc_assay",
              label = "Assay",
              choices = names(obj@assays),
              selected = Seurat::DefaultAssay(obj)
            ),
            selectInput(
              "group",
              label = "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            ),
            selectInput(
              "qc_palette",
              label = "Palette",
              choices = c("A","B","C","D","E","F","G","H"),
              selected = "C"
            ),
            numericInput("qc_violin_width", "Violin width", value = 0.8, min = 0, step = 0.1),
            sliderInput("qc_violin_alpha", "Violin alpha", min = 0, max = 1, value = 0.3, step = 0.1),
            numericInput("qc_box_width", "Box width", value = 0.3, min = 0, step = 0.1),
            sliderInput("qc_box_alpha", "Box alpha", min = 0, max = 1, value = 0.5, step = 0.1)
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
            numericInput("topn", "Top N Markers", value = 5, min = 1),
            numericInput("logfc", "logFC threshold", value = 0.25, step = 0.05),
            numericInput("minpct", "Min percent", value = 0.1, step = 0.05),
            selectInput("test", "Test method", choices = c("wilcox","t","LR"), selected = "wilcox"),
            selectInput("markers_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C")
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
            ),
            selectInput("batch_reduction", "Reduction", choices = c("pca","umap"), selected = "pca"),
            numericInput("batch_k", "Neighbors k", value = 20, min = 1),
            selectInput("batch_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C"),
            numericInput("batch_violin_width", "Violin width", value = 0.8, step = 0.1),
            sliderInput("batch_violin_alpha", "Violin alpha", min = 0, max = 1, value = 0.3, step = 0.1),
            numericInput("batch_box_width", "Box width", value = 0.3, step = 0.1),
            sliderInput("batch_box_alpha", "Box alpha", min = 0, max = 1, value = 0.5, step = 0.1)
          )
        ), column(
          width = 9,
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
            ),
            selectInput("trend_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C"),
            numericInput("trend_point_size", "Point size", value = 2, min = 0.1, step = 0.1),
            sliderInput("trend_point_alpha", "Point alpha", min = 0, max = 1, value = 0.3, step = 0.1),
            sliderInput("trend_smooth_alpha", "Smooth alpha", min = 0, max = 1, value = 0.3, step = 0.1),
            numericInput("trend_smooth_linewidth", "Smooth linewidth", value = 1.5, min = 0.1, step = 0.1)
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
              "lr_assay",
              "Assay",
              choices = names(obj@assays),
              selected = Seurat::DefaultAssay(obj)
            ),
            selectInput(
              "lr_group",
              "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            ),
            selectInput("lr_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C"),
            sliderInput("lr_tile_alpha", "Tile alpha", min = 0, max = 1, value = 0.8, step = 0.1)
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
            textInput("setA", "Set A Genes", value = paste0("G", 1:10, collapse = ",")),
            textInput("setB", "Set B Genes", value = paste0("G", 11:20, collapse = ",")),
            selectInput(
              "ms_group",
              "Group by",
              choices = colnames(obj@meta.data),
              selected = "seurat_clusters"
            ),
            numericInput("ms_nbin", "nbin", value = 24, min = 1),
            numericInput("ms_min_size", "min.size", value = 3, min = 1),
            selectInput("ms_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C"),
            numericInput("ms_violin_width", "Violin width", value = 0.8, step = 0.1),
            sliderInput("ms_violin_alpha", "Violin alpha", min = 0, max = 1, value = 0.3, step = 0.1),
            numericInput("ms_box_width", "Box width", value = 0.3, step = 0.1),
            sliderInput("ms_box_alpha", "Box alpha", min = 0, max = 1, value = 0.5, step = 0.1)
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
            textInput("sGenes", "S-phase Genes", value = paste0("G", 1:10, collapse = ",")),
            textInput("g2mGenes", "G2M-phase Genes", value = paste0("G", 11:20, collapse = ",")),
            selectInput("cc_reduction", "Reduction", choices = c("umap","pca"), selected = "umap"),
            numericInput("cc_dims", "Dims", value = 10, min = 2),
            selectInput("cc_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C"),
            sliderInput("cc_alpha", "Alpha", min = 0, max = 1, value = 0.8, step = 0.1)
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
            ),
            sliderInput("ec_levels", "Levels", min = 1, max = 10, value = 5, step = 1),
            numericInput("ec_point_size", "Point size", value = 1, min = 0.1, step = 0.1),
            sliderInput("ec_point_alpha", "Point alpha", min = 0, max = 1, value = 0.5, step = 0.1),
            sliderInput("ec_contour_alpha", "Contour alpha", min = 0, max = 1, value = 0.1, step = 0.1),
            selectInput("ec_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C")
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
            textInput("ch_genes", "Genes (comma)", value = paste0("G", 1:12, collapse = ",")),
            numericInput(
              "ch_thr",
              "Correlation Threshold",
              value = 0.2,
              min = 0,
              max = 1,
              step = 0.05
            ),
            numericInput("ch_point_size", "Point size", value = 3, min = 0.1, step = 0.1),
            sliderInput("ch_point_alpha", "Point alpha", min = 0, max = 1, value = 0.8, step = 0.1),
            numericInput("ch_label_size", "Label size", value = 3, min = 0.1, step = 0.1),
            sliderInput("ch_curve_alpha", "Curve alpha", min = 0, max = 1, value = 0.5, step = 0.1),
            selectInput("ch_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C")
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
            textInput("sp_feats", "Features (comma)", value = "G1,G2,G3"),
            numericInput("sp_point_size", "Point size", value = 2, min = 0.1, step = 0.1),
            sliderInput("sp_alpha", "Alpha", min = 0, max = 1, value = 0.5, step = 0.1),
            selectInput("sp_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C")
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
            ),
            selectInput("hex_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C")
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
            width = 12,
            height = NULL,
            textInput("lisa_gene", "Gene", value = "G10"),
            sliderInput(
              "lisa_k",
              "Neighbors (k)",
              min = 5,
              max = 50,
              value = 15,
              step = 5
            ),
            selectInput("lisa_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C"),
            numericInput("lisa_point_size", "Point size", value = 2, min = 0.1, step = 0.1),
            sliderInput("lisa_point_alpha", "Point alpha", min = 0, max = 1, value = 0.8, step = 0.1)
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
            ),
            selectInput("flow_palette", "Palette", choices = c("A","B","C","D","E","F","G","H"), selected = "C"),
            numericInput("flow_point_size", "Point size", value = 7, min = 0.1, step = 0.1),
            sliderInput("flow_point_alpha", "Point alpha", min = 0, max = 1, value = 0.9, step = 0.1),
            numericInput("flow_label_size", "Label size", value = 5, min = 0.1, step = 0.1)
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
  parse_gene_input <- function(s) {
    if (is.null(s)) return(character())
    if (length(s) > 1) return(trimws(s))
    x <- trimws(s)
    if (grepl("^paste0\\(['\"]G['\"],\\s*\\d+:\\d+\\)$", x)) {
      a <- as.integer(sub(".*,(\\s*)(\\d+):(\\d+).*", "\\2", x))
      b <- as.integer(sub(".*,(\\s*)(\\d+):(\\d+).*", "\\3", x))
      return(paste0("G", a:b))
    }
    if (grepl("^c\\((?:\\s*['\"]G\\d+['\"]\\s*,?)+\\)$", x)) {
      g <- regmatches(x, gregexpr("G\\d+", x))[[1]]
      return(g)
    }
    v <- trimws(unlist(strsplit(x, ",")))
    v <- v[v != ""]
    v
  }
  output$qc_plot <- renderPlot({
    p <- VisQCPanel(
      obj,
      assay = input$qc_assay,
      genes_mt = input$mt,
      genes_ribo = input$ribo,
      group.by = input$group,
      interactive = FALSE,
      palette = input$qc_palette,
      violin_width = input$qc_violin_width,
      violin_alpha = input$qc_violin_alpha,
      box_width = input$qc_box_width,
      box_alpha = input$qc_box_alpha
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
    r <- VisMarkerAtlas(
      obj,
      markers_top = input$topn,
      logfc_threshold = input$logfc,
      min_percent = input$minpct,
      test_method = input$test,
      palette = input$markers_palette
    )
    print(r$plot)
  })
  output$marker_table <- renderDT({
    r <- VisMarkerAtlas(
      obj,
      markers_top = input$topn,
      logfc_threshold = input$logfc,
      min_percent = input$minpct,
      test_method = input$test,
      palette = input$markers_palette
    )
    datatable(r$markers, options = list(scrollY = "150px", pageLength = 10))
  })

  output$batch_plot <- renderPlot({
    r <- VisBatchAlign(
      obj,
      batch = input$batch,
      reduction = input$batch_reduction,
      dims = 1:10,
      k = input$batch_k,
      palette = input$batch_palette,
      violin_width = input$batch_violin_width,
      violin_alpha = input$batch_violin_alpha,
      box_width = input$batch_box_width,
      box_alpha = input$batch_box_alpha
    )
    print(r$plot)
  })
  output$batch_table <- renderDT({
    r <- VisBatchAlign(
      obj,
      batch = input$batch,
      reduction = input$batch_reduction,
      dims = 1:10,
      k = input$batch_k,
      palette = input$batch_palette,
      violin_width = input$batch_violin_width,
      violin_alpha = input$batch_violin_alpha,
      box_width = input$batch_box_width,
      box_alpha = input$batch_box_alpha
    )
    datatable(r$summary, options = list(scrollY = "150px", pageLength = 10))
  })

  output$trend_plot <- renderPlot({
    feats <- parse_gene_input(input$genes)
    p <- VisGeneTrend(
      obj,
      features = feats,
      by = input$trendBy,
      reduction = 'umap',
      dims = 1:2,
      smooth.method = 'loess',
      palette = input$trend_palette,
      point_size = input$trend_point_size,
      point_alpha = input$trend_point_alpha,
      smooth_alpha = input$trend_smooth_alpha,
      smooth_linewidth = input$trend_smooth_linewidth
    )
    print(p)
  })

  lr_demo <- reactiveVal(data.frame(
    ligand = paste0('G', 1:5),
    receptor = paste0('G', 6:10)
  ))
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
      assay = input$lr_assay,
      lr_table = lr,
      group.by = input$lr_group,
      palette = input$lr_palette,
      tile_alpha = input$lr_tile_alpha
    )
    print(r$plot)
  })
  output$lr_table <- renderDT({
    lr <- lr_demo()
    req(lr)
    r <- VisLigRec(
      obj,
      assay = input$lr_assay,
      lr_table = lr,
      group.by = input$lr_group,
      palette = input$lr_palette,
      tile_alpha = input$lr_tile_alpha
    )
    datatable(r$scores, options = list(scrollY = "150px", pageLength = 10))
  })

  output$ms_plot <- renderPlot({
    setA <- parse_gene_input(input$setA)
    setB <- parse_gene_input(input$setB)
    r <- VisMetaFeature(
      obj,
      feature_sets = list(SetA = setA, SetB = setB),
      group.by = input$ms_group,
      nbin = input$ms_nbin,
      min.size = input$ms_min_size,
      palette = input$ms_palette,
      violin_width = input$ms_violin_width,
      violin_alpha = input$ms_violin_alpha,
      box_width = input$ms_box_width,
      box_alpha = input$ms_box_alpha
    )
    print(r$plot)
  })

  output$cc_plot <- renderPlot({
    s.genes <- parse_gene_input(input$sGenes)
    g2m.genes <- parse_gene_input(input$g2mGenes)
    r <- VisCellCycle(
      obj,
      genes_s = s.genes,
      genes_g2m = g2m.genes,
      reduction = input$cc_reduction,
      dims = 1:input$cc_dims,
      palette = input$cc_palette,
      alpha = input$cc_alpha
    )
    print(r$plot)
  })

  output$ec_plot <- renderPlot({
    print(
      VisEmbeddingContour(
        obj,
        group.by = input$ec_group,
        reduction = 'umap',
        levels = input$ec_levels,
        palette = input$ec_palette,
        point_size = input$ec_point_size,
        point_alpha = input$ec_point_alpha,
        contour_alpha = input$ec_contour_alpha
      )
    )
  })



  output$ch_plot <- renderPlot({
    genes <- parse_gene_input(input$ch_genes)
    print(
      VisGeneCoexpHive(
        obj,
        genes = genes,
        reduction = 'pca',
        threshold = input$ch_thr,
        palette = input$ch_palette,
        point_size = input$ch_point_size,
        point_alpha = input$ch_point_alpha,
        label_size = input$ch_label_size,
        curve_alpha = input$ch_curve_alpha
      )
    )
  })

  output$sp_plot <- renderPlot({
    feats <- parse_gene_input(input$sp_feats)
    print(
      VisSpatialOverlay(
        obj,
        features = feats,
        image = NULL,
        coords_cols = c('x', 'y'),
        palette = input$sp_palette,
        point_size = input$sp_point_size,
        alpha = input$sp_alpha
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
        palette = input$hex_palette
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
        palette = input$lisa_palette,
        point_size = input$lisa_point_size,
        point_alpha = input$lisa_point_alpha
      )
    )
  })
  output$flow_plot <- renderPlot({
    print(
      VisClusterFlowGraph(
        obj,
        group.by = input$flow_group,
        reduction = 'umap',
        palette = input$flow_palette,
        point_size = input$flow_point_size,
        point_alpha = input$flow_point_alpha,
        label_size = input$flow_label_size
      )
    )
  })
}

# Run app
shinyApp(ui, server)
