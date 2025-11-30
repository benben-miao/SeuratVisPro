 #' @title Local Moran's I hotspot for gene expression
 #' @description Compute local Moran's I on UMAP/PCA neighborhoods for a given gene to detect spatial autocorrelation hotspots.
 #' @param object A `Seurat` object.
 #' @param gene Gene name.
 #' @param reduction Reduction name, default 'umap'.
 #' @param k Number of nearest neighbors.
#' @return A ggplot scatter of local Moran's I values.
#' @export
#' @examples
#' obj <- SeuratVisProExample(n_cells = 300, n_genes = 800, n_clusters = 3)
#' p <- VisLocalMoran(obj, gene = 'G10', k = 15)
#' p
VisLocalMoran <- function(object, gene, reduction = 'umap', k = 15) {
   if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (is.null(object@reductions[[reduction]])) {
    if (reduction == 'pca') object <- suppressMessages(Seurat::RunPCA(object))
    if (reduction == 'umap') object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(object, dims = 1:10)))
  }
  emb <- Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE]
  assay <- Seurat::DefaultAssay(object)
  data <- Seurat::GetAssayData(object, assay = assay, slot = 'data')
  if (!(gene %in% rownames(data))) stop('gene not found in assay data')
  x <- Seurat::FetchData(object, vars = gene)[,1]
   dmat <- as.matrix(stats::dist(emb))
   n <- nrow(dmat)
   W <- matrix(0, n, n)
   for (i in seq_len(n)) {
     ord <- order(dmat[i, ])
     ord <- ord[ord != i]
     nei <- ord[seq_len(min(k, n-1))]
     W[i, nei] <- 1
   }
   W <- W / pmax(rowSums(W), 1)
   z <- (x - mean(x, na.rm = TRUE))
   s2 <- stats::var(x, na.rm = TRUE)
   local_I <- (z / s2) * (W %*% z)
   df <- as.data.frame(emb)
   colnames(df) <- c('X1','X2')
   df$I <- local_I[,1]
  ggplot2::ggplot(df, ggplot2::aes(x = X1, y = X2, color = I)) + ggplot2::geom_point(size = 0.8) + ggplot2::scale_color_viridis_c() +
    svpp_theme() + ggplot2::labs(x = paste0(toupper(reduction),'1'), y = paste0(toupper(reduction),'2'), color = "Local Moran's I")
 }
