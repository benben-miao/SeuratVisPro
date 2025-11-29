#' @title Embedding density contours per group
#' @description Overlay 2D density contours for groups on UMAP/PCA to highlight boundaries.
#' @param object A `Seurat` object with UMAP or PCA.
#' @param group.by Metadata column for grouping.
#' @param reduction Reduction name, default 'umap'.
#' @param levels Number of contour levels.
#' @param palette Optional manual palette for groups.
#' @return A ggplot with contours.
#' @examples
#' obj <- SeuratVisProExample()
#' VisEmbeddingContour(obj, group.by = 'seurat_clusters')
#' @export
VisEmbeddingContour <- function(object, group.by = 'seurat_clusters', reduction = 'umap', levels = 5, palette = NULL) {
  svpp_check_seurat_object(object)
  if (is.null(object@reductions[[reduction]])) {
    if (reduction == 'pca') object <- Seurat::RunPCA(object)
    if (reduction == 'umap') object <- Seurat::RunUMAP(object, dims = 1:10)
  }
  emb <- as.data.frame(Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE])
  colnames(emb) <- c('X1','X2')
  emb$group <- object@meta.data[[group.by]]
  p <- ggplot2::ggplot(emb, ggplot2::aes(x = X1, y = X2, color = group)) +
    ggplot2::geom_point(alpha = 0.2, size = 0.5) +
    ggplot2::stat_density2d(ggplot2::aes(fill = ..level..), geom = 'polygon', contour = TRUE, alpha = 0.25) +
    svpp_theme() + ggplot2::guides(fill = 'none') + ggplot2::labs(x = paste0(toupper(reduction),'1'), y = paste0(toupper(reduction),'2'))
  if (!is.null(palette)) p <- p + ggplot2::scale_color_manual(values = palette)
  p
}
