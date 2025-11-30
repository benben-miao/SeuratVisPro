#' @title Ranked dotmap for gene markers
#' @description A dotmap where genes are ranked per group by avg expression; dot size = pct, color = avg.
#' @param object A `Seurat` object.
#' @param group.by Metadata column for grouping.
#' @param features Optional gene list; if NULL, uses top markers.
#' @param top_n Number of top genes per group when features=NULL.
#' @param palette Optional color palette for avg expression.
#' @return A ggplot dotmap.
#' @examples
#' obj <- SeuratVisProExample()
#' VisRankedDotmap(obj, group.by = 'seurat_clusters', top_n = 5)
#' @export
VisRankedDotmap <- function(object, group.by = 'seurat_clusters', features = NULL, top_n = 8, palette = NULL) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!(group.by %in% colnames(object@meta.data))) stop('group.by not found')
  if (is.null(features)) {
    Seurat::Idents(object) <- object[[group.by]][, 1]
    mk <- try(Seurat::FindAllMarkers(object, only.pos = TRUE), silent = TRUE)
    if (inherits(mk, 'try-error') || nrow(mk) == 0) {
      avg <- Seurat::AverageExpression(object, group.by = group.by)
      mat0 <- avg[[Seurat::DefaultAssay(object)]]
      features <- unique(unlist(apply(mat0, 2, function(col) {
        utils::head(rownames(mat0)[order(col, decreasing = TRUE)], top_n)
      })))
    } else {
      features <- mk |> dplyr::group_by(cluster) |> dplyr::slice_head(n = top_n) |> dplyr::pull(gene) |> unique()
    }
  }
  avg <- Seurat::AverageExpression(object, group.by = group.by)
  mat <- avg[[Seurat::DefaultAssay(object)]]
  feats <- intersect(features, rownames(mat))
  if (length(feats) == 0) stop('No features available for dotmap after intersection with expression matrix')
  mat <- mat[feats, , drop = FALSE]
  groups <- unique(object@meta.data[[group.by]])
  fdat <- Seurat::FetchData(object, vars = feats)
  pct_mat <- sapply(groups, function(g) {
    cells <- rownames(object@meta.data)[object@meta.data[[group.by]] == g]
    if (length(cells) == 0) return(rep(NA_real_, length(feats)))
    colMeans(fdat[cells, , drop = FALSE] > 0, na.rm = TRUE)
  })
  colnames(pct_mat) <- groups
  rownames(pct_mat) <- feats
  df <- mat |> as.data.frame() |> tibble::rownames_to_column('gene') |> tidyr::pivot_longer(-gene, names_to = 'group', values_to = 'avg')
  df_pct <- pct_mat |> as.data.frame() |> tibble::rownames_to_column('gene') |> tidyr::pivot_longer(-gene, names_to = 'group', values_to = 'pct')
  df <- dplyr::left_join(df, df_pct, by = c('gene','group'))
  df <- df[!is.na(df$pct) & !is.na(df$avg), ]
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = gene, size = pct, color = avg)) + ggplot2::geom_point() + ggplot2::scale_size_area(max_size = 8) +
    svpp_theme() + ggplot2::labs(size = 'Pct', color = 'AvgExpr')
  if (!is.null(palette)) p <- p + ggplot2::scale_color_viridis_c(option = palette)
  p
}
