#' @title Cluster dendrogram visualization
#' @description Build a hierarchical tree of clusters using average expression profiles.
#' @author benben-miao
#'
#' @return A ggplot dendrogram (or patchwork with heatmap).
#' @param object A `Seurat` object; required.
#' @param group.by Metadata column for cluster identity. Default: `"seurat_clusters"`.
#' @param assay Assay name. Default: `Seurat::DefaultAssay(object)`.
#' @param dist.metric Distance metric: `'euclidean'` or `'correlation'`. Default: `'euclidean'`.
#' @param linkage Linkage method (e.g., `'complete'`, `'average'`, `'ward.D2'`). Default: `'complete'`.
#' @param show_heatmap Append pairwise similarity heatmap under the tree. Default: `TRUE`.
#'
#' @param palette Viridis palette option for color/fill. Default: `"C"`.
#' @param tile_alpha Tile alpha. Default: `0.8`.
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
#' p <- VisClusterTree(
#'   obj,
#'   assay = NULL,
#'   group.by = "seurat_clusters",
#'   dist.metric = "euclidean",
#'   linkage = "complete",
#'   show_heatmap = TRUE,
#'   palette = "C",
#'   tile_alpha = 0.8)
#' p
#'
VisClusterTree <- function(object,
                           assay = NULL,
                           group.by = "seurat_clusters",
                           dist.metric = "euclidean",
                           linkage = "complete",
                           show_heatmap = TRUE,
                           palette = "C",
                           tile_alpha = 0.8) {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # DefaultAssay
  if (is.null(assay))
    assay <- Seurat::DefaultAssay(object)

  # Groups in object@meta.data
  if (!(group.by %in% colnames(object@meta.data)))
    stop("group.by not found in meta.data")

  # AverageExpression
  avg <- Seurat::AverageExpression(object, assays = assay, group.by = group.by)
  mat <- as.matrix(avg[[assay]])
  if (dist.metric == "correlation") {
    cmat <- stats::cor(mat)
    dist_mat <- stats::as.dist(1 - cmat)
  } else {
    dist_mat <- stats::dist(t(mat))
  }

  # HClust
  hc <- stats::hclust(dist_mat, method = linkage)
  dd <- ggdendro::dendro_data(hc)

  # Plot
  p_tree <- ggplot2::ggplot(ggdendro::segment(dd), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = xend, yend = yend),
      stat = "identity",
      position = "identity",
      arrow = NULL,
      arrow.fill = NULL,
      lineend = "butt",
      linejoin = "round",
      na.rm = FALSE,
      show.legend = NA,
      inherit.aes = TRUE,
      linewidth = 1
    ) +
    # ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "Clusters", y = "Height") +
    ggplot2::scale_x_continuous(breaks = dd$labels$x,
                                labels = dd$labels$label) +
    ggplot2::scale_color_viridis_c(option = palette) +
    ggplot2::scale_fill_viridis_c(option = palette) +
    svpp_theme()

  if (!show_heatmap)
    return(p_tree)
  # heatmap of pairwise similarity
  sim_df <- as.data.frame(if (dist.metric == "correlation")
    stats::cor(mat)
    else
      as.matrix(dist_mat))
  if (dist.metric != "correlation") {
    # convert distance to similarity
    m <- max(sim_df)
    sim_df <- 1 - sim_df / m
  }

  sim_df <- sim_df |>
    tibble::rownames_to_column("cluster1") |>
    tidyr::pivot_longer(-cluster1, names_to = "cluster2", values_to = "sim")

  p_heat <- ggplot2::ggplot(sim_df, ggplot2::aes(x = cluster1, y = cluster2, fill = sim)) +
    ggplot2::geom_tile(alpha = tile_alpha) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Similarity") +
    ggplot2::scale_color_viridis_c(option = palette) +
    ggplot2::scale_fill_viridis_c(option = palette) +
    svpp_theme()

  patchwork::wrap_plots(p_tree, p_heat, ncol = 2)
}
