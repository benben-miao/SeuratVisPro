#' @title Module score visualization for gene sets
#' @description Compute module scores for multiple gene sets and visualize across groups.
#' @param object A `Seurat` object.
#' @param feature_sets A named list of character vectors of genes.
#' @param group.by Metadata column for grouping.
#' @param nbin Number of bins for `AddModuleScore`.
#' @return A list with updated object and ggplot boxplots.
#' @param min.size Minimum gene set size to keep.
#' @param palette Optional fill palette.
#' @examples
#' obj <- SeuratVisProExample()
#' sets <- list(SetA = paste0('G', 1:10), SetB = paste0('G', 11:20))
#' res <- VisMetaFeature(obj, feature_sets = sets, group.by = 'seurat_clusters')
#' res$plot
#' @export
VisMetaFeature <- function(object, feature_sets, group.by = "seurat_clusters", nbin = 24, min.size = 3, palette = NULL) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!(group.by %in% colnames(object@meta.data))) stop("group.by not found in meta.data")
  # filter sets to genes present
  assay <- Seurat::DefaultAssay(object)
  present_genes <- rownames(object[[assay]])
  feature_sets <- lapply(feature_sets, function(v) unique(intersect(v, present_genes)))
  # drop too-small sets
  feature_sets <- feature_sets[unlist(lapply(feature_sets, length)) >= min.size]
  if (length(feature_sets) == 0) stop("No valid feature sets with sufficient genes present")
  # try AddModuleScore, fallback to scaled mean expression
  add_ok <- TRUE
  try({ object <- Seurat::AddModuleScore(object, features = feature_sets, name = names(feature_sets), nbin = nbin) }, silent = TRUE)
  if (!any(grepl(paste0("^", paste(names(feature_sets), collapse = "|^")), colnames(object@meta.data)))) {
    add_ok <- FALSE
    data <- Seurat::GetAssayData(object, assay = assay, slot = "data")
    scores <- lapply(feature_sets, function(genes) {
      g <- intersect(genes, rownames(data))
      if (length(g) < min.size) return(rep(NA_real_, ncol(data)))
      colMeans(scale(as.matrix(data[g, , drop = FALSE])))
    })
    for (nm in names(scores)) {
      object[[nm]] <- scores[[nm]]
    }
  }
  md <- object@meta.data
  # collect module score columns created by AddModuleScore (prefixes of feature set names)
  score_cols <- grep(paste0("^", paste(names(feature_sets), collapse = "|^")), colnames(md), value = TRUE)
  df <- tibble::tibble(group = md[[group.by]])
  df[score_cols] <- md[score_cols]
  long <- df |> tidyr::pivot_longer(cols = dplyr::all_of(score_cols), names_to = "set", values_to = "score")
  p <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = score, fill = group)) + ggplot2::geom_violin(scale = "width", color = NA) + ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
    ggplot2::facet_wrap(~set, scales = "free_y") + ggplot2::theme_bw() + ggplot2::guides(fill = "none")
  if (!is.null(palette)) p <- p + ggplot2::scale_fill_manual(values = palette)
  list(object = object, plot = p)
}
