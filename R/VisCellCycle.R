#' @title Cell cycle visualization
#' @description Score cell cycle and visualize on embedding and distribution.
#' @param object A `Seurat` object.
#' @param s.genes Character vector for S phase genes.
#' @param g2m.genes Character vector for G2/M phase genes.
#' @param reduction Reduction name for overlay.
#' @param dims Dimensions for neighbor graph.
#' @return A list with updated object and patchwork plot.
#' @examples
#' obj <- SeuratVisProExample()
#' s.genes <- paste0('G', 1:10); g2m.genes <- paste0('G', 11:20)
#' res <- VisCellCycle(obj, s.genes, g2m.genes)
#' res$plot
#' @export
VisCellCycle <- function(object, s.genes, g2m.genes, reduction = "umap", dims = 1:10) {
  svpp_check_seurat_object(object)
  object <- suppressMessages(Seurat::CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE))
  if (is.null(object@reductions[[reduction]])) {
    if (reduction == "pca") object <- suppressMessages(Seurat::RunPCA(object))
    if (reduction == "umap") object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(object, dims = dims)))
  }
  p1 <- Seurat::DimPlot(object, reduction = reduction, group.by = "Phase") + svpp_theme()
  p2 <- ggplot2::ggplot(object@meta.data, ggplot2::aes(x = Phase)) + ggplot2::geom_bar() + svpp_theme()
  list(object = object, plot = patchwork::wrap_plots(p1, p2, ncol = 2))
}
