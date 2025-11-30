#' @title Cluster stability across resolutions
#' @description Assess clustering stability by resampling cells and varying resolution.
#' @author benben-miao
#'
#' @return A list with `summary` data.frame and a ggplot summarizing stability.
#' @param object A `Seurat` object with PCA computed.
#' @param resolution_range Numeric vector of resolutions, e.g., `seq(0.2, 1.2, by=0.2)`.
#' @param dims PCA dimensions to use.
#' @param reps Number of resampling repetitions.
#' @param prop Proportion of cells to sample per repetition.
#' @param palette Viridis palette option.
#'
#' @export
#'
#' @examples
#' obj <- SeuratVisProExample()
#' res <- VisClusterStability(obj, resolution_range = seq(0.2,1.2,0.2), dims = 1:10, reps = 5, prop = 0.8, palette = "C")
#' res$plot
#'
VisClusterStability <- function(object,
                                resolution_range = seq(0.2, 1.2, by = 0.2),
                                dims = 1:10,
                                reps = 5,
                                prop = 0.8,
                                palette = "C") {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # PCA in object@reductions$pca
  if (is.null(object@reductions$pca))
    object <- Seurat::RunPCA(object)

  # FindNeighbors
  object <- Seurat::FindNeighbors(object, dims = dims)

  # FindClusters
  full_clusters <- lapply(resolution_range, function(r) {
    as.character(Seurat::FindClusters(object, resolution = r)$seurat_clusters)
  })
  names(full_clusters) <- paste0("res_", resolution_range)

  # Stability
  stability <- dplyr::tibble(resolution = numeric(), agreement = numeric())
  cells <- colnames(object)
  for (r in resolution_range) {
    base_obj <- Seurat::FindClusters(object, resolution = r)
    base_assign <- as.character(base_obj$seurat_clusters)
    names(base_assign) <- colnames(base_obj)
    agree_vec <- numeric(reps)
    for (i in seq_len(reps)) {
      set.seed(i)
      sub_cells <- sample(cells, size = floor(length(cells) * prop))
      sub_obj <- subset(object, cells = sub_cells)
      sub_obj <- Seurat::FindNeighbors(sub_obj, dims = dims)
      sub_obj <- Seurat::FindClusters(sub_obj, resolution = r)
      sub_assign <- as.character(sub_obj$seurat_clusters)
      names(sub_assign) <- colnames(sub_obj)
      common <- intersect(names(base_assign), names(sub_assign))
      agree <- mean(base_assign[common] == sub_assign[common])
      agree_vec[i] <- agree
    }
    stability <- dplyr::bind_rows(stability, dplyr::tibble(resolution = r, agreement = mean(agree_vec)))
  }

  p <- ggplot2::ggplot(stability,
                       ggplot2::aes(
                         x = factor(resolution),
                         y = agreement,
                         fill = agreement
                       )) +
    ggplot2::geom_col(width = 0.5) +
    ggplot2::scale_color_viridis_c(option = palette) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Resolution", y = "Agreement (resampled)", fill = "Agreement") +
    svpp_theme()
  list(summary = stability, plot = p)
}
