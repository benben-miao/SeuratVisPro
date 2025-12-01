#' @title Hex-binned entropy map on embeddings
#' @description Compute per-bin Shannon entropy of group composition on UMAP/PCA and visualize as a tile heatmap.
#' @author benben-miao
#'
#' @return A ggplot heatmap of entropy.
#' @param object A `Seurat` object.
#' @param group.by Metadata column for groups.
#' @param reduction Reduction name, default 'umap'.
#' @param bins Number of bins per axis.
#'
#' @param palette Palette.
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
#' p <- VisHexEntropy(
#'   obj,
#'   group.by = "seurat_clusters",
#'   reduction = "umap",
#'   bins = 30,
#'   palette = "C")
#' p
#'
VisHexEntropy <- function(object,
                          group.by = "seurat_clusters",
                          reduction = "umap",
                          bins = 30,
                          palette = "C") {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # PCA or UMAP in object@reductions
  if (is.null(object@reductions[[reduction]])) {
    if (reduction == 'pca')
      object <- suppressMessages(Seurat::RunPCA(object))
    if (reduction == 'umap')
      object <- suppressMessages(Seurat::RunUMAP(object, dims = 1:10))
  }

  # Embeddings
  emb <- as.data.frame(Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE])
  colnames(emb) <- c('X1', 'X2')
  if (!(group.by %in% colnames(object@meta.data)))
    stop('group.by not found in meta.data')
  emb$group <- object@meta.data[[group.by]]
  if (!is.numeric(bins) ||
      bins < 2)
    stop('bins must be an integer >= 2')
  xr <- range(emb$X1, na.rm = TRUE)
  yr <- range(emb$X2, na.rm = TRUE)
  if (!is.finite(diff(xr)) ||
      diff(xr) == 0)
    xr <- scales::expand_range(xr, mul = 0.01)
  if (!is.finite(diff(yr)) ||
      diff(yr) == 0)
    yr <- scales::expand_range(yr, mul = 0.01)
  bx_breaks <- seq(xr[1], xr[2], length.out = bins + 1)
  by_breaks <- seq(yr[1], yr[2], length.out = bins + 1)
  emb$bx <- cut(
    emb$X1,
    breaks = bx_breaks,
    include.lowest = TRUE,
    right = TRUE
  )
  emb$by <- cut(
    emb$X2,
    breaks = by_breaks,
    include.lowest = TRUE,
    right = TRUE
  )
  agg <- emb |>
    dplyr::count(bx, by, group) |>
    tidyr::pivot_wider(
      names_from = group,
      values_from = n,
      values_fill = 0
    )
  probs <- as.matrix(agg[, -c(1, 2)])
  row_sums <- rowSums(probs)
  probs <- sweep(probs, 1, row_sums, FUN = '/')
  probs[is.na(probs)] <- 0
  entropy <- -rowSums(ifelse(probs > 0, probs * log(probs), 0))

  mx <- (bx_breaks[-1] + bx_breaks[-length(bx_breaks)]) / 2
  my <- (by_breaks[-1] + by_breaks[-length(by_breaks)]) / 2
  df <- agg[, c("bx", "by")]
  df$entropy <- entropy
  df$x <- mx[as.integer(df$bx)]
  df$y <- my[as.integer(df$by)]

  # Plot
  ggplot2::ggplot(df, ggplot2::aes(
    x = x,
    y = y,
    color = entropy,
    fill = entropy
  )) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      x = paste0(toupper(reduction), '1'),
      y = paste0(toupper(reduction), '2')
    ) +
    ggplot2::scale_color_viridis_c(option = palette) +
    ggplot2::scale_fill_viridis_c(option = palette) +
    svpp_theme()
}
