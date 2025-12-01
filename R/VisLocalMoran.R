#' @title Local Moran's I hotspot for gene expression
#' @description Compute local Moran's I on UMAP/PCA neighborhoods for a given gene to detect spatial autocorrelation hotspots.
#' @author benben-miao
#'
#' @return A ggplot scatter of local Moran's I values.
#' @param object A `Seurat` object.
#' @param gene Gene name.
#' @param reduction Reduction name, default 'umap'.
#' @param k Number of nearest neighbors.
#'
#' @param palette Palette.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
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
#' p <- VisLocalMoran(
#'    obj,
#'    gene = 'G10',
#'    reduction = "umap",
#'    k = 15,
#'    palette = "C",
#'    point_size = 2,
#'    point_alpha = 0.8)
#' p
#'
VisLocalMoran <- function(object,
                          gene,
                          reduction = "umap",
                          k = 15,
                          palette = "C",
                          point_size = 2,
                          point_alpha = 0.8) {
   # Seurat object check
   if (!inherits(object, "Seurat"))
      stop("object must be a Seurat object")

   # PCA or UMAP in object@reductions
   if (is.null(object@reductions[[reduction]])) {
      if (reduction == 'pca')
         object <- suppressMessages(Seurat::RunPCA(object))
      if (reduction == 'umap')
         object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(object, dims = 1:10)))
   }

   # Embeddings
   emb <- Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE]
   assay <- Seurat::DefaultAssay(object)
   data <- Seurat::GetAssayData(object, assay = assay, slot = 'data')
   if (!(gene %in% rownames(data)))
      stop('gene not found in assay data')
   x <- Seurat::FetchData(object, vars = gene)[, 1]
   dmat <- as.matrix(stats::dist(emb))
   n <- nrow(dmat)
   W <- matrix(0, n, n)
   for (i in seq_len(n)) {
      ord <- order(dmat[i, ])
      ord <- ord[ord != i]
      nei <- ord[seq_len(min(k, n - 1))]
      W[i, nei] <- 1
   }
   W <- W / pmax(rowSums(W), 1)
   z <- (x - mean(x, na.rm = TRUE))
   s2 <- stats::var(x, na.rm = TRUE)
   local_I <- (z / s2) * (W %*% z)
   df <- as.data.frame(emb)
   colnames(df) <- c('X1', 'X2')
   df$I <- local_I[, 1]

   # Plot
   ggplot2::ggplot(df, ggplot2::aes(x = X1, y = X2, color = I)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::labs(
         x = paste0(toupper(reduction), '1'),
         y = paste0(toupper(reduction), '2'),
         color = "Local Moran's I"
      ) +
      ggplot2::scale_color_viridis_c(option = palette) +
      ggplot2::scale_fill_viridis_c(option = palette) +
      svpp_theme()
}
