#' @title Embedding density contours per group
#' @description Overlay 2D density contours for groups on UMAP/PCA to highlight boundaries.
#' @author benben-miao
#'
#' @return A ggplot with contours.
#' @param object A `Seurat` object with UMAP or PCA.
#' @param group.by Metadata column for grouping.
#' @param reduction Reduction name, default 'umap'.
#' @param levels Number of contour levels.
#'
#' @param palette Optional manual palette for groups.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#' @param contour_alpha Contour alpha.
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
#' p <- VisEmbeddingContour(
#'   obj,
#'   group.by = "seurat_clusters",
#'   reduction = "umap",
#'   levels = 5,
#'   palette = "C",
#'   point_size = 1,
#'   point_alpha = 0.5,
#'   contour_alpha = 0.1)
#' p
#'
VisEmbeddingContour <- function(object,
                                group.by = "seurat_clusters",
                                reduction = "umap",
                                levels = 5,
                                palette = "C",
                                point_size = 1,
                                point_alpha = 0.5,
                                contour_alpha = 0.1) {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # PCA or UMAP in object@reductions
  if (is.null(object@reductions[[reduction]])) {
    if (reduction == 'pca')
      object <- Seurat::RunPCA(object)
    if (reduction == 'umap')
      object <- Seurat::RunUMAP(object, dims = 1:10)
  }

  # Embeddings
  emb <- as.data.frame(Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE])
  colnames(emb) <- c('X1', 'X2')
  emb$group <- object@meta.data[[group.by]]

  # Plot
  p <- ggplot2::ggplot(emb, ggplot2::aes(x = X1, y = X2, color = group)) +
    ggplot2::geom_point(alpha = point_alpha,
                        size = point_size,
                        show.legend = FALSE) +
    ggplot2::stat_density2d(
      ggplot2::aes(fill = ggplot2::after_stat(level)),
      geom = 'polygon',
      contour = TRUE,
      alpha = contour_alpha
    ) +
    ggplot2::guides(fill = 'none') +
    ggplot2::labs(x = paste0(toupper(reduction), '1'), y = paste0(toupper(reduction), '2')) +
    ggplot2::scale_color_viridis_d(option = palette) +
    ggplot2::scale_fill_viridis_c(option = palette) +
    svpp_theme()
  p
}
