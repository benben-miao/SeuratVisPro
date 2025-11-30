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
#'
#' @param palette Viridis palette option.
#' @examples
#' obj <- SeuratVisProExample()
#' p <- VisQCPanel(obj, genes_ribo = "^RPL|^RPS", group.by = "seurat_clusters")
#' p
#' @export
VisQCPanel <- function(object,
                       genes_mt = "^MT-",
                       genes_ribo = NULL,
                       group.by = NULL,
                       interactive = FALSE,
                       assay = NULL,
                       palette = "C") {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  mt_feats <- grep(genes_mt, rownames(object[[assay]]), value = TRUE)
  if (length(mt_feats) == 0) {
    object$percent.mt <- 0
  } else {
    object <- Seurat::PercentageFeatureSet(object, features = mt_feats, col.name = "percent.mt")
  }
  if (!is.null(genes_ribo)) {
    rb_feats <- grep(genes_ribo, rownames(object[[assay]]), value = TRUE)
    if (length(rb_feats) == 0) {
      object$percent.ribo <- 0
    } else {
      object <- Seurat::PercentageFeatureSet(object, features = rb_feats, col.name = "percent.ribo")
    }
  }
  md <- object@meta.data
  if (is.null(group.by) ||
      !(group.by %in% colnames(md)))
    group.by <- NULL

  xlab <- if (is.null(group.by))
    "seurat_clusters"
  else
    group.by
  v1 <- ggplot2::ggplot(md,
                        ggplot2::aes(
                          x = !!rlang::sym(xlab),
                          y = nFeature_RNA,
                          fill = !!rlang::sym(xlab)
                        )) +
    ggplot2::geom_violin(scale = "width", color = NA) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
    ggplot2::labs(x = xlab, y = "nFeature_RNA") +
    ggplot2::guides(fill = "none") + svpp_theme()
  v2 <- ggplot2::ggplot(md,
                        ggplot2::aes(
                          x = !!rlang::sym(xlab),
                          y = nCount_RNA,
                          fill = !!rlang::sym(xlab)
                        )) +
    ggplot2::geom_violin(scale = "width", color = NA) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
    ggplot2::labs(x = xlab, y = "nCount_RNA") +
    ggplot2::guides(fill = "none") + svpp_theme()
  v3 <- ggplot2::ggplot(md,
                        ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    ggplot2::geom_point(alpha = 0.6, size = 0.8) +
    ggplot2::scale_color_viridis_c(option = palette) +
    svpp_theme() +
    ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA", color = "%MT")
  if ("percent.ribo" %in% colnames(md)) {
    v4 <- ggplot2::ggplot(md, ggplot2::aes(x = percent.mt, y = percent.ribo)) +
      ggplot2::geom_bin2d(bins = 30) +
      svpp_theme() +
      ggplot2::labs(x = "%MT", y = "%Ribo")
  } else {
    v4 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Ribosomal genes not provided")
  }
  patch <- patchwork::wrap_plots(v1, v2, v3, v4, ncol = 2)
  if (interactive) {
    return(
      list(
        v1 = plotly::ggplotly(v1),
        v2 = plotly::ggplotly(v2),
        v3 = plotly::ggplotly(v3),
        v4 = plotly::ggplotly(v4)
      )
    )
  }
  patch
}
