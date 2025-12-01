#' @title Interactive QC panel for Seurat v5 objects
#' @description Visualize key QC metrics with customizable thresholds and optional interactivity.
#' @author benben-miao
#'
#' @return A assembled ggplot or a list of interactive plots.
#' @param object A `Seurat` object.
#' @param genes_mt Regex for mitochondrial genes. Default `^MT-` for human.
#' @param genes_ribo Regex for ribosomal genes (e.g., `^RPL|^RPS`).
#' @param group.by Column in `object@meta.data` used for grouping.
#' @param interactive If `TRUE`, returns interactive plots via `plotly::ggplotly`.
#' @param assay Assay to use; defaults to `DefaultAssay(object)`.
#' @param palette Viridis palette option.
#' @param violin_width Violin width.
#' @param violin_alpha Violin alpha.
#' @param box_width Box width.
#' @param box_alpha Box alpha.
#'
#' @export
#'
#' @examples
#' obj <- SeuratVisProExample(
#'     n_cells = 300,
#'     n_genes = 1000,
#'     n_clusters = 10,
#'     seed = 123,
#'     genes_mt = "^MT-",
#'     neighbor_dims = 10,
#'     cluster_res = 0.5,
#'     umap_dims = 10,
#'     spatial = FALSE)
#'
#' p <- VisQCPanel(
#'   obj,
#'   assay = NULL,
#'   genes_mt = "^MT-",
#'   genes_ribo = NULL,
#'   group.by = "seurat_clusters",
#'   interactive = FALSE,
#'   palette = "C",
#'   violin_width = 0.8,
#'   violin_alpha = 0.3,
#'   box_width = 0.3,
#'   box_alpha = 0.5)
#' p
#'
VisQCPanel <- function(object,
                       assay = NULL,
                       genes_mt = "^MT-",
                       genes_ribo = NULL,
                       group.by = NULL,
                       interactive = FALSE,
                       palette = "C",
                       violin_width = 0.8,
                       violin_alpha = 0.3,
                       box_width = 0.3,
                       box_alpha = 0.5) {
  # Check object class
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # DefaultAssay
  if (is.null(assay))
    assay <- Seurat::DefaultAssay(object)

  # Genes in assay
  assay_obj <- Seurat::GetAssay(object, assay = assay)
  genes <- rownames(assay_obj)

  # PercentageFeatureSet for mitochondrial genes
  mt_feats <- grep(genes_mt, genes, value = TRUE)
  if (length(mt_feats) == 0) {
    n_cells <- ncol(object)
    object$percent.mt <- rep(0, n_cells)
  } else {
    object <- Seurat::PercentageFeatureSet(object, features = mt_feats, col.name = "percent.mt")
  }

  # PercentageFeatureSet for ribosomal genes
  if (!is.null(genes_ribo)) {
    rb_feats <- grep(genes_ribo, genes, value = TRUE)
    if (length(rb_feats) == 0) {
      n_cells <- ncol(object)
      object$percent.ribo <- rep(0, n_cells)
    } else {
      object <- Seurat::PercentageFeatureSet(object, features = rb_feats, col.name = "percent.ribo")
    }
  }

  # Group in meta.data
  md <- object@meta.data
  if (!is.null(group.by) && group.by %in% colnames(md)) {
    md$.group <- md[[group.by]]
    group_col <- ".group"
  } else if (!is.null(Seurat::Idents(object))) {
    md$.group <- as.character(Seurat::Idents(object))
    group_col <- ".group"
  } else {
    md$.group <- factor("all")
    group_col <- ".group"
  }
  xlab <- ifelse(is.null(group.by), "group", group.by)

  # Plot
  v1 <- ggplot2::ggplot(
    md,
    ggplot2::aes(
      x = !!rlang::sym(group_col),
      y = nFeature_RNA,
      fill = !!rlang::sym(group_col),
      color = !!rlang::sym(group_col)
    )
  ) +
    ggplot2::geom_violin(
      width = violin_width,
      scale = "width",
      color = NA,
      alpha = violin_alpha,
      show.legend = TRUE
    ) +
    ggplot2::geom_boxplot(
      width = box_width,
      alpha = box_alpha,
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::labs(x = "Cell Cluster", y = "nFeature RNA") +
    ggplot2::guides() +
    svpp_theme()

  v2 <- ggplot2::ggplot(
    md,
    ggplot2::aes(
      x = !!rlang::sym(group_col),
      y = nCount_RNA,
      fill = !!rlang::sym(group_col),
      color = !!rlang::sym(group_col)
    )
  ) +
    ggplot2::geom_violin(
      width = violin_width,
      scale = "width",
      color = NA,
      alpha = violin_alpha,
      show.legend = TRUE
    ) +
    ggplot2::geom_boxplot(
      width = box_width,
      alpha = box_alpha,
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::labs(x = "Cell Cluster", y = "nCount RNA") +
    ggplot2::guides() +
    svpp_theme()

  v3 <- ggplot2::ggplot(md,
                        ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::scale_color_viridis_c(option = palette) +
    ggplot2::labs(x = "nCount RNA", y = "nFeature RNA", color = "MT Percent") +
    svpp_theme()

  if ("percent.ribo" %in% colnames(md)) {
    v4 <- ggplot2::ggplot(md, ggplot2::aes(x = percent.mt, y = percent.ribo)) +
      ggplot2::geom_bin2d(bins = 30) +
      ggplot2::labs(x = "%MT", y = "%Ribo") +
      svpp_theme()
  } else {
    v4 <- ggplot2::ggplot() +
      ggplot2::labs(title = "Ribosomal genes not provided") +
      svpp_theme()
  }
  patch <- patchwork::wrap_plots(v1, v2, v3, v4, ncol = 2)

  # Plotly
  if (interactive) {
    p1 <- plotly::ggplotly(v1)
    p2 <- plotly::ggplotly(v2)
    p3 <- plotly::ggplotly(v3)
    p4 <- plotly::ggplotly(v4)
    return(plotly::subplot(
      p1,
      p2,
      p3,
      p4,
      nrows = 2,
      shareX = FALSE,
      shareY = FALSE
    ))
  }
  patch
}
