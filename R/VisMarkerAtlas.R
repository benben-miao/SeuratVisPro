#' @title Marker atlas visualization
#' @description Find top markers per cluster and visualize as a heatmap.
#' @author benben-miao
#'
#' @return A list containing marker table and ggplot heatmap.
#' @param object A `Seurat` object.
#' @param markers_top Number of top markers per cluster.
#' @param logfc_threshold Log fold-change threshold for marker finding.
#' @param min_percent Minimum percent expressed.
#' @param test_method Differential test (e.g., 'wilcox').
#'
#' @param palette Viridis palette option.
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
#' res <- VisMarkerAtlas(
#'   obj,
#'   markers_top = 5,
#'   logfc_threshold = 0.25,
#'   min_percent = 0.1,
#'   test_method = "wilcox",
#'   palette = "C")
#'
#' res$plot
#' head(res$markers)
#'
VisMarkerAtlas <- function(object,
                           markers_top = 5,
                           logfc_threshold = 0.25,
                           min_percent = 0.1,
                           test_method = "wilcox",
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
    logfc.threshold = logfc_threshold,
    min.pct = min_percent,
    test.use = test_method
  )

  # Markers top N
  top_markers <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = markers_top) |>
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
    ggplot2::labs(x = "Cluster", y = "Marker genes", fill = "AvgExpr") +
    ggplot2::scale_color_viridis_c(option = palette) +
    ggplot2::scale_fill_viridis_c(option = palette) +
    svpp_theme()

  list(markers = top_markers, plot = p)
}
