#' @title Radial hive plot of gene co-expression
#' @description Visualize a radial network of selected genes; edge weight reflects co-expression across cells.
#' @author benben-miao
#'
#' @return A ggplot radial hive-like network.
#' @param object A `Seurat` object; required.
#' @param genes Character vector of genes; must contain ≥3 present in the assay.
#' @param reduction Reduction used to position genes radially by loading (`'pca'`). Default: `'pca'`.
#' @param threshold Minimum absolute correlation to draw edges. Default: `0.2`.
#'
#' @param palette Viridis palette option for color/fill. Default: `"C"`.
#' @param point_size Node point size. Default: `3`.
#' @param point_alpha Node point alpha. Default: `0.8`.
#' @param label_size Label size. Default: `3`.
#' @param curve_alpha Edge curve alpha. Default: `0.5`.
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
#' p <- VisGeneCoexpHive(
#'   obj,
#'   genes = paste0("G", 1:100),
#'   reduction = "pca",
#'   threshold = 0.2,
#'   palette = "C",
#'   point_size = 3,
#'   point_alpha = 0.8,
#'   label_size = 3,
#'   curve_alpha = 0.5)
#' p
#'
VisGeneCoexpHive <- function(object,
                             genes,
                             reduction = "pca",
                             threshold = 0.2,
                             palette = "C",
                             point_size = 3,
                             point_alpha = 0.8,
                             label_size = 3,
                             curve_alpha = 0.5) {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # DefaultAssay
  assay <- Seurat::DefaultAssay(object)
  data <- Seurat::GetAssayData(object, assay = assay, slot = 'data')

  genes <- intersect(genes, rownames(data))
  if (length(genes) < 3)
    stop('Need >=3 genes for hive plot')

  # PCA
  if (is.null(object@reductions$pca))
    object <- Seurat::RunPCA(object)
  loadings <- Seurat::Loadings(object[['pca']])[, 1]
  theta <- seq(0, 2 * pi, length.out = length(genes) + 1)[1:length(genes)]
  pos <- data.frame(gene = genes,
                    theta = theta,
                    r = scales::rescale(abs(loadings[genes]), to = c(0.4, 1)))

  expr <- t(as.matrix(data[genes, , drop = FALSE]))
  cmat <- stats::cor(expr)
  edges <- which(abs(cmat) >= threshold &
                   upper.tri(cmat), arr.ind = TRUE)

  if (nrow(edges) == 0)
    stop('No edges above threshold')

  edge_df <- data.frame(g1 = genes[edges[, 1]], g2 = genes[edges[, 2]], w = cmat[edges])
  pos2 <- pos
  rownames(pos2) <- pos2$gene
  edge_df$x0 <- pos2[edge_df$g1, 'r'] * cos(pos2[edge_df$g1, 'theta'])
  edge_df$y0 <- pos2[edge_df$g1, 'r'] * sin(pos2[edge_df$g1, 'theta'])
  edge_df$x1 <- pos2[edge_df$g2, 'r'] * cos(pos2[edge_df$g2, 'theta'])
  edge_df$y1 <- pos2[edge_df$g2, 'r'] * sin(pos2[edge_df$g2, 'theta'])
  nodes <- transform(pos, x = r * cos(theta), y = r * sin(theta))

  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = edge_df,
      ggplot2::aes(
        x = x0,
        y = y0,
        xend = x1,
        yend = y1,
        color = w,
        size = abs(w)
      ),
      curvature = 0.5,
      alpha = curve_alpha
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y),
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = gene),
      vjust = -0.8,
      size = label_size
    ) +
    ggplot2::scale_size(range = c(0.5, 3)) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_viridis_c(option = palette) +
    ggplot2::scale_fill_viridis_c(option = palette) +
    svpp_theme()
  p
}
