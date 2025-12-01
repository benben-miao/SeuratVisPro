#' @title Cell cycle visualization
#' @description Score cell cycle and visualize on embedding and distribution.
#' @author benben-miao
#'
#' @return A list with updated object and patchwork plot.
#' @param object A `Seurat` object.
#' @param genes_s Character vector for S phase genes.
#' @param genes_g2m Character vector for G2/M phase genes.
#' @param reduction Reduction name for overlay.
#' @param dims Dimensions for neighbor graph.
#'
#' @param palette Palette.
#' @param alpha Points and bars alpha.
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
#' genes_s <- paste0('G', 1:10); genes_g2m <- paste0('G', 11:20)
#'
#' res <- VisCellCycle(
#'   obj,
#'   genes_s,
#'   genes_g2m,
#'   reduction = "umap",
#'   dims = 1:10,
#'   palette = "C",
#'   alpha = 0.8)
#'
#' res$plot
#' res$cluster
#' res$bar
#'
VisCellCycle <- function(object,
                         genes_s,
                         genes_g2m,
                         reduction = "umap",
                         dims = 1:10,
                         palette = "C",
                         alpha = 0.8) {
  # Seurat object chek
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # CellCycleScoring
  object <- suppressMessages(
    Seurat::CellCycleScoring(
      object,
      s.features = genes_s,
      g2m.features = genes_g2m,
      set.ident = TRUE
    )
  )

  # PCA or UMAP in object@reductions
  if (is.null(object@reductions[[reduction]])) {
    if (reduction == "pca")
      object <- suppressMessages(Seurat::RunPCA(object))
    if (reduction == "umap")
      object <- suppressMessages(Seurat::RunUMAP(object, dims = dims))
  }

  # Plot
  p1 <- Seurat::DimPlot(object,
                        reduction = reduction,
                        group.by = "Phase",
                        alpha = alpha) +
    ggplot2::scale_color_viridis_d(option = palette) +
    ggplot2::scale_fill_viridis_d(option = palette) +
    svpp_theme()

  p2 <- ggplot2::ggplot(object@meta.data, ggplot2::aes(x = Phase, fill = Phase)) +
    ggplot2::geom_bar(width = 0.5, alpha = alpha) +
    ggplot2::scale_color_viridis_d(option = palette) +
    ggplot2::scale_fill_viridis_d(option = palette) +
    svpp_theme()

  list(
    object = object,
    plot = patchwork::wrap_plots(p1, p2, ncol = 2),
    cluster = p1,
    bar = p2
  )
}
