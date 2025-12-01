#' @title Generate `Seurat` object for demonstrations
#' @description Generate a `Seurat` object with basic preprocessing (DefaultAssay, PercentageFeatureSet, NormalizeData, FindVariableFeatures, ScaleData, RunPCA, FindNeighbors, FindClusters, RunUMAP).
#' @author benben-miao
#'
#' @return A `Seurat` object with basic preprocessing steps applied.
#' @param n_cells Number of cells. Default: `300`.
#' @param n_genes Number of genes. Default: `1000`.
#' @param n_clusters Number of clusters. Default: `10`.
#' @param seed Random seed for reproducibility. Default: `123`.
#' @param genes_mt Regex for mitochondrial genes. Default: `"^MT-"`.
#' @param neighbor_dims Neighbor graph dimensions used in `FindNeighbors`. Default: `10`.
#' @param cluster_res Cluster resolution used in `FindClusters`. Default: `0.5`.
#' @param umap_dims UMAP dimensions used in `RunUMAP`. Default: `10`.
#' @param spatial Whether to add synthetic spatial coordinates `x`,`y` to `meta.data`. Default: `FALSE`.
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
#' Seurat::DimPlot(obj, group.by = "cluster")
#'
SeuratVisProExample <- function(n_cells = 300,
								n_genes = 1000,
								n_clusters = 10,
								seed = 123,
								genes_mt = "^MT-",
								neighbor_dims = 10,
								cluster_res = 0.5,
								umap_dims = 10,
								spatial = FALSE) {
	# Seed
	set.seed(seed)

	# Data
	clusters <- sample(paste0("C", seq_len(n_clusters)), n_cells, replace = TRUE)
	base_counts <- matrix(stats::rpois(n_cells * n_genes, lambda = 1), n_genes, n_cells)
	rownames(base_counts) <- paste0("G", seq_len(n_genes))
	colnames(base_counts) <- paste0("Cell", seq_len(n_cells))
	for (k in seq_len(n_clusters)) {
		idx <- which(clusters == paste0("C", k))
		gene_shift <- ((k - 1) * floor(n_genes / n_clusters) + 1):min(n_genes, k * floor(n_genes / n_clusters))
		base_counts[gene_shift, idx] <- base_counts[gene_shift, idx] + matrix(stats::rpois(length(gene_shift) * length(idx), lambda = 2),
																			  length(gene_shift),
																			  length(idx))
	}

	# Convert to sparse to avoid coercion warning
	base_counts <- Matrix::Matrix(base_counts, sparse = TRUE)
	obj <- Seurat::CreateSeuratObject(base_counts)
	obj$cluster <- clusters
	assay <- Seurat::DefaultAssay(obj)
	mt_feats <- grep(genes_mt, rownames(obj[[assay]]), value = TRUE)
	if (length(mt_feats) == 0) {
		obj$percent.mt <- 0
	} else {
		obj <- Seurat::PercentageFeatureSet(obj, features = mt_feats, col.name = "percent.mt")
	}

	# Pipeline with default parameters
	obj <- suppressMessages(Seurat::NormalizeData(obj))
	obj <- suppressMessages(Seurat::FindVariableFeatures(obj))
	obj <- suppressMessages(Seurat::ScaleData(obj))
	obj <- suppressMessages(Seurat::RunPCA(obj))
	obj <- suppressMessages(Seurat::FindNeighbors(obj, dims = 1:neighbor_dims))
	obj <- suppressMessages(Seurat::FindClusters(obj, resolution = cluster_res))
	obj <- suppressMessages(Seurat::RunUMAP(obj, dims = 1:umap_dims))
	if (spatial) {
		coords <- matrix(stats::rnorm(n_cells * 2), ncol = 2)
		obj$x <- coords[, 1]
		obj$y <- coords[, 2]
	}
	obj
}
