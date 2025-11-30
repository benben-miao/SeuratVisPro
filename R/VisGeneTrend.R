#' @title Gene expression trends along pseudotime or groups
#' @description Plot smoothed trends for selected genes across pseudotime approximated from embeddings or grouped by clusters.
#' @param object A `Seurat` object with UMAP or PCA.
#' @param features Character vector of gene names.
#' @param by Either 'pseudotime' or a metadata column to group by.
#' @param reduction Reduction to use for pseudotime ('umap' or 'pca').
#' @param dims Dimensions used for pseudotime ranking.
#' @param smooth.method Smoothing method ('loess' or 'gam').
#' @param palette Viridis palette option.
#' @return A ggplot with smoothed curves.
#' @examples
#' obj <- SeuratVisProExample()
#' p <- VisGeneTrend(obj, features = c('G10','G20'), by = 'pseudotime')
#' p
#' @export
VisGeneTrend <- function(object, features, by = "pseudotime", reduction = "umap", dims = 1:2, smooth.method = "loess", palette = "C") {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  assay <- Seurat::DefaultAssay(object)
  if (by == "pseudotime") {
    if (is.null(object@reductions[[reduction]])) {
      if (reduction == "pca") object <- suppressMessages(Seurat::RunPCA(object))
      if (reduction == "umap") object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(object)))
    }
    emb <- Seurat::Embeddings(object[[reduction]])[, dims, drop = FALSE]
    pt <- scale(emb[, 1])[, 1]
    object$pseudotime <- rank(pt) / length(pt)
    by_col <- "pseudotime"
  } else {
    if (!(by %in% colnames(object@meta.data))) stop("Grouping column not found in meta.data")
    by_col <- by
  }
  df <- Seurat::FetchData(object, vars = c(features, by_col))
  names(df)[ncol(df)] <- "group"
  df <- df |> tidyr::pivot_longer(cols = all_of(features), names_to = "gene", values_to = "expr")
  if (smooth.method == "loess") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = expr, color = gene)) + ggplot2::geom_point(alpha = 0.2, size = 0.5) + ggplot2::geom_smooth(method = "loess") + svpp_theme() +
      ggplot2::labs(x = if (by == "pseudotime") "Pseudotime" else by, y = "Expression")
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = expr, color = gene)) + ggplot2::geom_point(alpha = 0.2, size = 0.5) + ggplot2::geom_smooth(method = "gam") + svpp_theme() +
      ggplot2::labs(x = if (by == "pseudotime") "Pseudotime" else by, y = "Expression")
  }
  p
}
