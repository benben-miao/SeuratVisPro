#' @title Module score visualization for gene sets
#' @description Compute module scores for multiple gene sets and visualize across groups.
#' @author benben-miao
#'
#' @return A list with `object` and `plot` (violin+box per gene set).
#' @param object A `Seurat` object; required.
#' @param feature_sets A named list of character vectors of genes; required.
#' @param group.by Metadata column for grouping. Default: `"seurat_clusters"`.
#' @param nbin Number of bins for `Seurat::AddModuleScore`. Default: `24`.
#' @param min.size Minimum gene set size to keep. Default: `3`.
#'
#' @param palette Viridis palette option for color/fill. Default: `"C"`.
#' @param violin_width Violin width. Default: `0.8`.
#' @param violin_alpha Violin alpha. Default: `0.3`.
#' @param box_width Box width. Default: `0.3`.
#' @param box_alpha Box alpha. Default: `0.5`.
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
#' sets <- list(
#'   SetA = paste0('G', 1:10),
#'   SetB = paste0('G', 11:20))
#'
#' res <- VisMetaFeature(
#'   obj,
#'   feature_sets = sets,
#'   group.by = "seurat_clusters",
#'   nbin = 24,
#'   min.size = 3,
#'   palette = "C",
#'   violin_width = 0.8,
#'   violin_alpha = 0.3,
#'   box_width = 0.3,
#'   box_alpha = 0.5)
#'
#' res$plot
#'
VisMetaFeature <- function(object,
                           feature_sets,
                           group.by = "seurat_clusters",
                           nbin = 24,
                           min.size = 3,
                           palette = "C",
                           violin_width = 0.8,
                           violin_alpha = 0.3,
                           box_width = 0.3,
                           box_alpha = 0.5) {
  # Seurat object check
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")

  # Groups in object@meta.data
  if (!(group.by %in% colnames(object@meta.data)))
    stop("group.by not found in meta.data")

  # Filter sets to genes present
  assay <- Seurat::DefaultAssay(object)
  present_genes <- rownames(object[[assay]])
  feature_sets <- lapply(feature_sets, function(v)
    unique(intersect(v, present_genes)))

  # Drop too-small sets
  feature_sets <- feature_sets[unlist(lapply(feature_sets, length)) >= min.size]
  if (length(feature_sets) == 0)
    stop("No valid feature sets with sufficient genes present")

  # Try AddModuleScore, fallback to scaled mean expression
  add_ok <- TRUE
  try({
    object <- Seurat::AddModuleScore(
      object,
      features = feature_sets,
      name = names(feature_sets),
      nbin = nbin
    )
  }, silent = TRUE)
  if (!any(grepl(paste0("^", paste(
    names(feature_sets), collapse = "|^"
  )), colnames(object@meta.data)))) {
    add_ok <- FALSE
    data <- Seurat::GetAssayData(object, assay = assay, layer = "data")
    scores <- lapply(feature_sets, function(genes) {
      g <- intersect(genes, rownames(data))
      if (length(g) < min.size)
        return(rep(NA_real_, ncol(data)))
      colMeans(scale(as.matrix(data[g, , drop = FALSE])))
    })
    for (nm in names(scores)) {
      object[[nm]] <- scores[[nm]]
    }
  }
  md <- object@meta.data

  # Collect module score columns created by AddModuleScore (prefixes of feature set names)
  score_cols <- grep(paste0("^", paste(names(feature_sets), collapse = "|^")), colnames(md), value = TRUE)
  df <- tibble::tibble(group = md[[group.by]])
  df[score_cols] <- md[score_cols]
  long <- df |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(score_cols),
      names_to = "set",
      values_to = "score"
    )

  # Plot
  p <- ggplot2::ggplot(long,
                       ggplot2::aes(
                         x = group,
                         y = score,
                         fill = group,
                         color = group
                       )) +
    ggplot2::geom_violin(
      width = violin_width,
      scale = "width",
      alpha = violin_alpha,
      color = NA,
      show.legend = TRUE
    ) +
    ggplot2::geom_boxplot(
      width = box_width,
      alpha = box_alpha,
      linewidth = 1,
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.2e", x)
    ) +
    ggplot2::facet_wrap(~ set, scales = "free_y") +
    ggplot2::guides() +
    ggplot2::labs(x = "Group", y = "Score") +
    ggplot2::scale_color_viridis_d(option = palette) +
    ggplot2::scale_fill_viridis_d(option = palette) +
    svpp_theme()

  list(object = object, plot = p)
}
