#' @title Marker atlas visualization
#' @description Find top markers per cluster and visualize as a heatmap.
#' @author benben-miao
#'
#' @return A list containing marker table and ggplot heatmap.
#' @param object A `Seurat` object.
#' @param top_n Number of top markers per cluster.
#' @param logfc.threshold Log fold-change threshold for marker finding.
#' @param min.pct Minimum percent expressed.
#' @param test.use Differential test (e.g., 'wilcox').
#' @param palette Viridis palette option.
#'
#' @export
#'
#' @examples
#' obj <- SeuratVisProExample()
#'
#' res <- VisMarkerAtlas(obj, top_n = 5, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox", palette = "C")
#' res$plot
#' head(res$markers)
#'
VisMarkerAtlas <- function(object,
                           top_n = 5,
                           logfc.threshold = 0.25,
                           min.pct = 0.1,
                           test.use = "wilcox",
                           palette = "C") {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # Seurat clusters
  Seurat::Idents(object) <- object$seurat_clusters

  # FindAllMarkers
  markers <- Seurat::FindAllMarkers(
    object,
    only.pos = TRUE,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    test.use = test.use
  )

  # Markers top N
  top_markers <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()

  feats <- unique(top_markers$gene)

  # AverageExpression
  mat <- Seurat::AverageExpression(object, features = feats)$RNA
  mat <- mat[feats, , drop = FALSE]

  df <- mat |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "avg_exp")

  df$gene <- factor(df$gene, levels = rev(feats))

  # Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = gene, fill = avg_exp)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = palette) +
    ggplot2::labs(x = "Cluster", y = "Marker genes", fill = "AvgExpr") +
    svpp_theme()

  list(markers = top_markers, plot = p)
}
