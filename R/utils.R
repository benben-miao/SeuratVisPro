#' @keywords internal
svpp_check_seurat_object <- function(object) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  invisible(TRUE)
}

#' @keywords internal
svpp_default_mt_pattern <- function() {
  "^MT-"
}

#' @keywords internal
svpp_compute_percent <- function(object, pattern, assay = NULL, name = "percent") {
  svpp_check_seurat_object(object)
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  feats <- grep(pattern, rownames(object[[assay]]), value = TRUE)
  if (length(feats) == 0) {
    object[[name]] <- 0
  } else {
    object <- Seurat::PercentageFeatureSet(object, features = feats, col.name = name)
  }
  object
}

#' @title Example Seurat object for demonstrations
#' @description Generate a small synthetic Seurat object suitable for examples and testing.
#' @param n_cells Number of cells.
#' @param n_genes Number of genes.
#' @param n_clusters Number of simulated clusters.
#' @param seed Random seed for reproducibility.
#' @param spatial If TRUE, adds synthetic spatial coordinates 'x','y' to meta.data.
#' @return A `Seurat` object with basic preprocessing (normalization, variable features, PCA, UMAP).
#' @examples
#' obj <- SeuratVisProExample(n_cells = 400, n_genes = 1000, n_clusters = 4)
#' Seurat::DimPlot(obj, group.by = "cluster")
#' @export
SeuratVisProExample <- function(n_cells = 300, n_genes = 1000, n_clusters = 3, seed = 123, spatial = FALSE) {
  set.seed(seed)
  clusters <- sample(paste0("C", seq_len(n_clusters)), n_cells, replace = TRUE)
  base_counts <- matrix(stats::rpois(n_cells * n_genes, lambda = 1), n_genes, n_cells)
  rownames(base_counts) <- paste0("G", seq_len(n_genes))
  colnames(base_counts) <- paste0("Cell", seq_len(n_cells))
  for (k in seq_len(n_clusters)) {
    idx <- which(clusters == paste0("C", k))
    gene_shift <- ((k - 1) * floor(n_genes / n_clusters) + 1):min(n_genes, k * floor(n_genes / n_clusters))
    base_counts[gene_shift, idx] <- base_counts[gene_shift, idx] + matrix(stats::rpois(length(gene_shift) * length(idx), lambda = 2), length(gene_shift), length(idx))
  }
  # convert to sparse to avoid coercion warning
  base_counts <- Matrix::Matrix(base_counts, sparse = TRUE)
  obj <- Seurat::CreateSeuratObject(base_counts)
  obj$cluster <- clusters
  obj <- svpp_compute_percent(obj, svpp_default_mt_pattern(), name = "percent.mt")
  obj <- suppressMessages(Seurat::NormalizeData(obj))
  obj <- suppressMessages(Seurat::FindVariableFeatures(obj))
  obj <- suppressMessages(Seurat::ScaleData(obj))
  obj <- suppressMessages(Seurat::RunPCA(obj))
  obj <- suppressMessages(Seurat::FindNeighbors(obj, dims = 1:10))
  obj <- suppressMessages(Seurat::FindClusters(obj, resolution = 0.5))
  obj <- suppressWarnings(suppressMessages(Seurat::RunUMAP(obj, dims = 1:10)))
  if (spatial) {
    coords <- matrix(stats::rnorm(n_cells * 2), ncol = 2)
    obj$x <- coords[, 1]
    obj$y <- coords[, 2]
  }
  obj
}

#' @keywords internal
svpp_theme <- function(grid = TRUE) {
  fam <- "Arial"
  psfonts <- names(grDevices::postscriptFonts())
  if (!is.null(psfonts) && !(fam %in% psfonts)) fam <- "sans"
  if (!interactive()) fam <- "sans"
  th <- ggplot2::theme_bw(base_size = 12, base_family = fam)
  if (!grid) th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  th + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(face = "bold")
  )
}
