#' @title Spatial feature overlay
#' @description Overlay selected features onto spatial images with improved color scales.
#' @author benben-miao
#'
#' @return A patchwork of spatial feature plots.
#' @param object A `Seurat` object; if not spatial, falls back to embedding overlay.
#' @param features Character vector of gene features.
#' @param image Image key; defaults to active image.
#' @param coords_cols Column names for fallback coordinates.
#'
#' @param palette Viridis palette option for color/fill.
#' @param point_size Point size for fallback scatter.
#' @param alpha Point alpha.
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
#'     spatial = TRUE)
#'
#' p <- VisSpatialOverlay(
#'   obj,
#'   features = c("G1", "G2", "G3", "G4"),
#'   image = NULL,
#'   coords_cols = c("x", "y"),
#'   palette = "C",
#'   point_size = 2,
#'   alpha = 0.5)
#' p
#'
VisSpatialOverlay <- function(object,
                              features,
                              image = NULL,
                              coords_cols = c("x", "y"),
                              palette = "C",
                              point_size = 2,
                              alpha = 0.5) {
  # Seurat object
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # SpatialFeaturePlot
  if (!is.null(object@images) && length(object@images) > 0) {
    plots <- lapply(features, function(f) {
      p <- Seurat::SpatialFeaturePlot(object, features = f, image = image) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(title = f)) +
        ggplot2::scale_fill_viridis_c(option = palette) +
        svpp_theme()
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
      ggplot2::ggplot(df,
                      ggplot2::aes(
                        x = !!rlang::sym(coords_cols[1]),
                        y = !!rlang::sym(coords_cols[2]),
                        color = !!rlang::sym(f)
                      )) +
        ggplot2::geom_point(alpha = alpha, size = point_size) +
        ggplot2::labs(color = f, title = f) +
        ggplot2::scale_color_viridis_c(option = palette) +
        svpp_theme()

    })
    return(patchwork::wrap_plots(plots, ncol = 2))
  }
  if (is.null(object@reductions$umap) &&
      is.null(object@reductions$pca)) {
    object <- suppressMessages(Seurat::RunPCA(object))
    object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(object, dims = 1:10)))
  }
  red <- if (!is.null(object@reductions$umap))
    "umap"
  else
    "pca"
  plots <- lapply(features, function(f) {
    Seurat::FeaturePlot(object,
                        features = f,
                        reduction = red,
                        order = TRUE) +
      ggplot2::labs(color = f) +
      ggplot2::scale_color_viridis_c(option = palette) +
      svpp_theme()

  })
  patchwork::wrap_plots(plots, ncol = 2)
}
