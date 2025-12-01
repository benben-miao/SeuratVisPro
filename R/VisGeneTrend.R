#' @title Gene expression trends along pseudotime or groups
#' @description Plot smoothed trends for selected genes across pseudotime approximated from embeddings or grouped by clusters.
#' @author benben-miao
#'
#' @return A ggplot with smoothed curves.
#' @param object A `Seurat` object with UMAP or PCA.
#' @param features Character vector of gene names.
#' @param by Either 'pseudotime' or a metadata column to group by.
#' @param reduction Reduction to use for pseudotime ('umap' or 'pca').
#' @param dims Dimensions used for pseudotime ranking.
#' @param smooth.method Smoothing method ('loess' or 'gam').
#' @param palette Viridis palette option.
#' @param point_size Point size.
#' @param smooth_alpha Smooth alpha.
#' @param smooth_linewidth Smooth line width.
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
#' p <- VisGeneTrend(
#'   obj,
#'   features = c("G10", "G20", "G30"),
#'   by = "pseudotime",
#'   reduction = "umap",
#'   dims = 1:2,
#'   smooth.method = "loess",
#'   palette = "C",
#'   point_size = 2,
#'   smooth_alpha = 0.3,
#'   smooth_linewidth = 1.5)
#' p
#'
VisGeneTrend <- function(object,
                         features,
                         by = "pseudotime",
                         reduction = "umap",
                         dims = 1:2,
                         smooth.method = "loess",
                         palette = "C",
                         point_size = 2,
                         smooth_alpha = 0.3,
                         smooth_linewidth = 1.5) {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # DefaultAssay
  assay <- Seurat::DefaultAssay(object)

  # PCA or UMAP in object@reductions
  if (by == "pseudotime") {
    if (is.null(object@reductions[[reduction]])) {
      if (reduction == "pca")
        object <- suppressMessages(Seurat::RunPCA(object))
      if (reduction == "umap")
        object <- suppressMessages(Seurat::RunUMAP(object))
    }

    # Embeddings
    emb <- Seurat::Embeddings(object[[reduction]])[, dims, drop = FALSE]
    pt <- scale(emb[, 1])[, 1]
    object$pseudotime <- rank(pt) / length(pt)
    by_col <- "pseudotime"
  } else {
    if (!(by %in% colnames(object@meta.data)))
      stop("Grouping column not found in meta.data")
    by_col <- by
  }

  # FetchData
  df <- Seurat::FetchData(object, vars = c(features, by_col))
  names(df)[ncol(df)] <- "group"
  df <- df |>
    tidyr::pivot_longer(cols = all_of(features),
                        names_to = "gene",
                        values_to = "expr")
  if (smooth.method == "loess") {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = group,
      y = expr,
      fill = gene,
      color = gene
    )) +
      ggplot2::geom_point(alpha = 0.5, size = point_size) +
      ggplot2::geom_smooth(
        method = "loess",
        alpha = smooth_alpha,
        linewidth = smooth_linewidth,
        show.legend = TRUE
      ) +
      ggplot2::labs(x = if (by == "pseudotime")
        "Pseudotime"
        else
          by, y = "Expression") +
      svpp_theme()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = group,
      y = expr,
      fill = gene,
      color = gene
    )) +
      ggplot2::geom_point(alpha = 0.5, size = point_size) +
      ggplot2::geom_smooth(
        method = "gam",
        alpha = smooth_alpha,
        linewidth = smooth_linewidth,
        show.legend = TRUE
      ) +
      ggplot2::labs(x = if (by == "pseudotime")
        "Pseudotime"
        else
          by, y = "Expression") +
      svpp_theme()
  }
  p
}
