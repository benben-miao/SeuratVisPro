#' @title Marker atlas visualization
#' @description Find top markers per cluster and visualize as a heatmap.
#' @param object A `Seurat` object.
#' @param top_n Number of top markers per cluster.
#' @param logfc.threshold Log fold-change threshold for marker finding.
#' @param min.pct Minimum percent expressed.
#' @param test.use Differential test (e.g., 'wilcox').
#' @return A list containing marker table and ggplot heatmap.
#' @param palette Viridis palette option.
#' @examples
#' obj <- SeuratVisProExample()
#' res <- VisMarkerAtlas(obj, top_n = 5)
#' res$plot
#' @export
VisMarkerAtlas <- function(object, top_n = 10, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox", palette = "C") {
  svpp_check_seurat_object(object)
  Seurat::Idents(object) <- object$seurat_clusters
  markers <- Seurat::FindAllMarkers(object, only.pos = TRUE, logfc.threshold = logfc.threshold, min.pct = min.pct, test.use = test.use)
  top_markers <- markers |> dplyr::group_by(cluster) |> dplyr::slice_head(n = top_n) |> dplyr::ungroup()
  feats <- unique(top_markers$gene)
  mat <- Seurat::AverageExpression(object, features = feats)$RNA
  mat <- mat[feats, , drop = FALSE]
  df <- mat |> as.data.frame() |> tibble::rownames_to_column("gene") |> tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "avg_exp")
  df$gene <- factor(df$gene, levels = rev(feats))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = gene, fill = avg_exp)) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = palette) + svpp_theme() +
    ggplot2::labs(x = "Cluster", y = "Marker genes", fill = "AvgExpr")
  list(markers = top_markers, plot = p)
}
