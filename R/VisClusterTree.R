#' @title Cluster dendrogram visualization
#' @description Build a hierarchical tree of clusters using average expression profiles.
#' @param object A `Seurat` object.
#' @param group.by Metadata column for cluster identity.
#' @param assay Assay name.
#' @param dist.metric Distance metric: 'euclidean' or 'correlation'.
#' @param linkage Linkage method, e.g., 'complete','average','ward.D2'.
#' @param label.size Label size for cluster names.
#' @param show_heatmap If TRUE, append a pairwise similarity heatmap under the tree.
#'
#' @return A ggplot dendrogram (or patchwork with heatmap).
#' @examples
#' obj <- SeuratVisProExample()
#' p <- VisClusterTree(obj, group.by = 'seurat_clusters')
#' p
#' @export
VisClusterTree <- function(object, group.by = "seurat_clusters", assay = NULL, dist.metric = "euclidean", linkage = "complete", label.size = 3, show_heatmap = TRUE) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  if (!(group.by %in% colnames(object@meta.data))) stop("group.by not found in meta.data")
  avg <- Seurat::AverageExpression(object, assays = assay, group.by = group.by)
  mat <- as.matrix(avg[[assay]])
  if (dist.metric == "correlation") {
    cmat <- stats::cor(mat)
    dist_mat <- stats::as.dist(1 - cmat)
  } else {
    dist_mat <- stats::dist(t(mat))
  }
  hc <- stats::hclust(dist_mat, method = linkage)
  dd <- ggdendro::dendro_data(hc)
  p_tree <- ggplot2::ggplot(ggdendro::segment(dd), ggplot2::aes(x = x, y = y)) + ggplot2::geom_segment(ggplot2::aes(xend = xend, yend = yend)) + svpp_theme() +
    ggplot2::scale_y_reverse() + ggplot2::labs(x = "Clusters", y = "Height") +
    ggplot2::scale_x_continuous(breaks = dd$labels$x, labels = dd$labels$label) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = label.size))
  if (!show_heatmap) return(p_tree)
  # heatmap of pairwise similarity
  sim_df <- as.data.frame(if (dist.metric == "correlation") stats::cor(mat) else as.matrix(dist_mat))
  if (dist.metric != "correlation") {
    # convert distance to similarity
    m <- max(sim_df)
    sim_df <- 1 - sim_df / m
  }
  sim_df <- sim_df |> tibble::rownames_to_column("cluster1") |> tidyr::pivot_longer(-cluster1, names_to = "cluster2", values_to = "sim")
  p_heat <- ggplot2::ggplot(sim_df, ggplot2::aes(x = cluster1, y = cluster2, fill = sim)) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c() + svpp_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + ggplot2::labs(x = NULL, y = NULL, fill = "Similarity")
  patchwork::wrap_plots(p_tree, p_heat, ncol = 1)
}
