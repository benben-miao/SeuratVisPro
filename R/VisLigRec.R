#' @title Ligand-receptor co-expression scoring
#' @description Compute simple ligand-receptor co-expression scores across groups and visualize as a heatmap.
#' @param object A `Seurat` object.
#' @param lr_table Data frame with columns `ligand` and `receptor`.
#' @param group.by Metadata column for groups (e.g., clusters).
#' @param assay Assay name.
#' @return A list with score table and ggplot heatmap.
#' @param palette Viridis palette option.
#' @examples
#' obj <- SeuratVisProExample()
#' lr <- data.frame(ligand = paste0('G', 1:5), receptor = paste0('G', 6:10))
#' res <- VisLigRec(obj, lr_table = lr, group.by = 'seurat_clusters')
#' res$plot
#' @export
VisLigRec <- function(object, lr_table, group.by = "seurat_clusters", assay = NULL, palette = "C") {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  if (!(group.by %in% colnames(object@meta.data))) stop("group.by not found in meta.data")
  groups <- object@meta.data[[group.by]]
  avg <- Seurat::AverageExpression(object, assays = assay, group.by = group.by)
  mat <- avg[[assay]]
  lr_table <- lr_table |> dplyr::filter(ligand %in% rownames(mat) & receptor %in% rownames(mat))
  res <- dplyr::tibble()
  for (g1 in colnames(mat)) {
    for (g2 in colnames(mat)) {
      score <- mean(mat[lr_table$ligand, g1] * mat[lr_table$receptor, g2])
      res <- dplyr::bind_rows(res, dplyr::tibble(source = g1, target = g2, score = score))
    }
  }
  p <- ggplot2::ggplot(res, ggplot2::aes(x = source, y = target, fill = score)) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = palette) + svpp_theme() +
    ggplot2::labs(x = "Source group (ligand)", y = "Target group (receptor)", fill = "LR score")
  list(scores = res, plot = p)
}
