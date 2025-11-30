#' @title Spatial feature overlay
#' @description Overlay selected features onto spatial images with improved color scales.
#' @param object A `Seurat` object; if not spatial, falls back to embedding overlay.
#' @param features Character vector of gene features.
#' @param image Image key; defaults to active image.
#' @param alpha Point alpha.
#' @param coords_cols Column names for fallback coordinates.
#' @param point_size Point size for fallback scatter.
#' @param palette Viridis palette option for color/fill.
#' @return A patchwork of spatial feature plots.
#' @examples
#' obj <- SeuratVisProExample(n_cells = 300, n_genes = 800, n_clusters = 3, spatial = TRUE)
#' VisSpatialOverlay(obj, features = c('G1','G2'))
#' @export
VisSpatialOverlay <- function(object, features, image = NULL, alpha = 0.8, coords_cols = c("x","y"), point_size = 0.8, palette = "C") {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!is.null(object@images) && length(object@images) > 0) {
    plots <- lapply(features, function(f) {
      p <- Seurat::SpatialFeaturePlot(object, features = f, image = image) +
        ggplot2::scale_fill_viridis_c(option = palette) + ggplot2::guides(fill = ggplot2::guide_colorbar(title = f)) + svpp_theme()
      p
    })
    return(patchwork::wrap_plots(plots, ncol = 2))
  }
  # fallback: embedding overlay using FeaturePlot
  md <- object@meta.data
  if (all(coords_cols %in% colnames(md))) {
    df <- Seurat::FetchData(object, vars = features)
    df[[coords_cols[1]]] <- md[[coords_cols[1]]]
    df[[coords_cols[2]]] <- md[[coords_cols[2]]]
    plots <- lapply(features, function(f) {
      ggplot2::ggplot(df, ggplot2::aes(x = !!rlang::sym(coords_cols[1]), y = !!rlang::sym(coords_cols[2]), color = !!rlang::sym(f))) +
        ggplot2::geom_point(alpha = alpha, size = point_size) + ggplot2::scale_color_viridis_c(option = palette) + svpp_theme() +
        ggplot2::labs(color = f, title = f)
    })
    return(patchwork::wrap_plots(plots, ncol = 2))
  }
  if (is.null(object@reductions$umap) && is.null(object@reductions$pca)) {
    object <- suppressMessages(Seurat::RunPCA(object))
    object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(object, dims = 1:10)))
  }
  red <- if (!is.null(object@reductions$umap)) "umap" else "pca"
  plots <- lapply(features, function(f) {
    Seurat::FeaturePlot(object, features = f, reduction = red, order = TRUE) +
      ggplot2::scale_color_viridis_c(option = palette) + svpp_theme() + ggplot2::labs(color = f)
  })
  patchwork::wrap_plots(plots, ncol = 2)
}
