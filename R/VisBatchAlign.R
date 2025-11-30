#' @title Batch mixing diagnostics
#' @description Quantify and visualize batch mixing in a kNN graph before/after integration.
#' @param object A `Seurat` object.
#' @param batch Column name in `meta.data` indicating batch.
#' @param reduction Dimensional reduction to use for kNN (e.g., 'pca' or 'umap').
#' @param dims Dimensions to use.
#' @param k Number of neighbors.
#' @param palette Viridis palette option.
#' @return A list with summary data and plot comparing neighbor batch composition.
#' @examples
#' obj <- SeuratVisProExample()
#' obj$batch <- sample(c('A','B'), ncol(obj), replace = TRUE)
#' res <- VisBatchAlign(obj, batch = 'batch', reduction = 'pca')
#' res$plot
#' @export
VisBatchAlign <- function(object, batch, reduction = "pca", dims = 1:10, k = 20, palette = "C") {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!(batch %in% colnames(object@meta.data))) stop("batch column not found in meta.data")
  if (is.null(object@reductions[[reduction]])) {
    if (reduction == "pca") object <- suppressMessages(Seurat::RunPCA(object))
    if (reduction == "umap") object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(object, dims = dims)))
  }
  emb <- Seurat::Embeddings(object[[reduction]])[, dims, drop = FALSE]
  dmat <- as.matrix(stats::dist(emb))
  n <- nrow(dmat)
  nn_idx <- matrix(NA_integer_, nrow = n, ncol = k)
  for (i in seq_len(n)) {
    ord <- order(dmat[i, ])
    ord <- ord[ord != i]
    nn_idx[i, ] <- ord[seq_len(k)]
  }
  batches <- object@meta.data[[batch]]
  mix <- numeric(nrow(nn_idx))
  for (i in seq_len(nrow(nn_idx))) {
    nei <- nn_idx[i, ]
    mix[i] <- mean(batches[nei] != batches[i])
  }
  df <- dplyr::tibble(cell = colnames(object), mix_prop = mix, batch = batches)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = batch, y = mix_prop, fill = batch)) + ggplot2::geom_violin(scale = "width", color = NA) + ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
    svpp_theme() + ggplot2::guides(fill = "none") + ggplot2::labs(y = "Neighbor batch difference proportion")
  list(summary = df, plot = p)
}
