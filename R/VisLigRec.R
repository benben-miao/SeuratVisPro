#' @title Ligand-receptor co-expression scoring
#' @description Compute simple ligand-receptor co-expression scores across groups and visualize as a heatmap.
#' @author benben-miao
#'
#' @return A list with score table and ggplot heatmap.
#' @param object A `Seurat` object.
#' @param lr_table Data frame with columns `ligand` and `receptor`.
#' @param group.by Metadata column for groups (e.g., clusters).
#' @param assay Assay name.
#' @param palette Viridis palette option.
#' @param tile_alpha Tile alpha.
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
#' lr <- data.frame(
#'   ligand = paste0('G', 1:5),
#'   receptor = paste0('G', 6:10))
#'
#' res <- VisLigRec(
#'   obj,
#'   lr_table = lr,
#'   group.by = "seurat_clusters")
#'
#' res$plot
#' head(res$scores)
#'
VisLigRec <- function(object,
                      lr_table,
                      group.by = "seurat_clusters",
                      assay = NULL,
                      palette = "C",
                      tile_alpha = 0.8) {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # DefaultAssay
  if (is.null(assay))
    assay <- Seurat::DefaultAssay(object)

  # Groups in object@meta.data
  if (!(group.by %in% colnames(object@meta.data)))
    stop("group.by not found in meta.data")

  groups <- object@meta.data[[group.by]]

  # AverageExpression
  avg <- Seurat::AverageExpression(object, assays = assay, group.by = group.by)
  mat <- avg[[assay]]
  lr_table <- lr_table |>
    dplyr::filter(ligand %in% rownames(mat) &
                    receptor %in% rownames(mat))

  res <- dplyr::tibble()
  for (g1 in colnames(mat)) {
    for (g2 in colnames(mat)) {
      score <- mean(mat[lr_table$ligand, g1] * mat[lr_table$receptor, g2])
      res <- dplyr::bind_rows(res, dplyr::tibble(
        source = g1,
        target = g2,
        score = score
      ))
    }
  }

  # Plot
  p <- ggplot2::ggplot(res, ggplot2::aes(x = source, y = target, fill = score)) +
    ggplot2::geom_tile(
      stat = "identity",
      position = "identity",
      lineend = "butt",
      linejoin = "mitre",
      na.rm = FALSE,
      show.legend = NA,
      inherit.aes = TRUE,
      alpha = tile_alpha
    ) +
    ggplot2::scale_fill_viridis_c(option = palette) +
    ggplot2::labs(x = "Source group (ligand)", y = "Target group (receptor)", fill = "LR score") +
    svpp_theme()
  list(scores = res, plot = p)
}
