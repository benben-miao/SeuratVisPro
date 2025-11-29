#' @title Radial hive plot of gene co-expression
#' @description Visualize a radial network of selected genes; edge weight reflects co-expression across cells.
#' @param object A `Seurat` object.
#' @param genes Character vector of genes.
#' @param reduction Reduction to position genes radially by loading (PCA).
#' @param threshold Minimum absolute correlation to draw edges.
#' @param low_color Low color for correlation.
#' @param mid_color Mid color for correlation.
#' @param low_color Low color for correlation.
#' @param mid_color Mid color for correlation.
#' @param high_color High color for correlation.
#' @return A ggplot radial hive-like network.
#' @examples
#' obj <- SeuratVisProExample()
#' VisGeneCoexpHive(obj, genes = paste0('G', 1:12), threshold = 0.2)
#' @export
VisGeneCoexpHive <- function(object, genes, reduction = 'pca', threshold = 0.2, low_color = 'steelblue', mid_color = 'grey80', high_color = 'firebrick') {
  svpp_check_seurat_object(object)
  assay <- Seurat::DefaultAssay(object)
  data <- Seurat::GetAssayData(object, assay = assay, slot = 'data')
  genes <- intersect(genes, rownames(data))
  if (length(genes) < 3) stop('Need >=3 genes for hive plot')
  if (is.null(object@reductions$pca)) object <- Seurat::RunPCA(object)
  loadings <- Seurat::Loadings(object[['pca']])[, 1]
  theta <- seq(0, 2*pi, length.out = length(genes)+1)[1:length(genes)]
  pos <- data.frame(gene = genes, theta = theta, r = scales::rescale(abs(loadings[genes]), to = c(0.4, 1)))
  expr <- t(as.matrix(data[genes, , drop = FALSE]))
  cmat <- stats::cor(expr)
  edges <- which(abs(cmat) >= threshold & upper.tri(cmat), arr.ind = TRUE)
  if (nrow(edges) == 0) stop('No edges above threshold')
  edge_df <- data.frame(g1 = genes[edges[,1]], g2 = genes[edges[,2]], w = cmat[edges])
  pos2 <- pos; rownames(pos2) <- pos2$gene
  edge_df$x0 <- pos2[edge_df$g1, 'r'] * cos(pos2[edge_df$g1, 'theta'])
  edge_df$y0 <- pos2[edge_df$g1, 'r'] * sin(pos2[edge_df$g1, 'theta'])
  edge_df$x1 <- pos2[edge_df$g2, 'r'] * cos(pos2[edge_df$g2, 'theta'])
  edge_df$y1 <- pos2[edge_df$g2, 'r'] * sin(pos2[edge_df$g2, 'theta'])
  nodes <- transform(pos, x = r*cos(theta), y = r*sin(theta))
  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(data = edge_df, ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1, color = w, size = abs(w)), curvature = 0.2, alpha = 0.6) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y), size = 3) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y, label = gene), vjust = -0.8, size = 3) +
    ggplot2::scale_color_gradient2(low = low_color, mid = mid_color, high = high_color) + ggplot2::scale_size(range = c(0.2, 2)) +
    ggplot2::coord_equal() + svpp_theme() + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  p
}
