#' @title Cluster centroid flow graph on embeddings
#' @description Build a minimum spanning tree over cluster centroids on UMAP/PCA and draw arrows to visualize putative transitions.
#' @author benben-miao
#'
#' @return A ggplot with centroid nodes and MST edges.
#' @param object A `Seurat` object.
#' @param group.by Metadata column for cluster identity.
#' @param reduction Reduction name, default 'umap'.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#' @param label_size Label size.
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
#' p <- VisClusterFlowGraph(
#'    obj,
#'    group.by = "seurat_clusters",
#'    reduction = "umap",
#'    point_size = 7,
#'    point_alpha = 0.9,
#'    label_size = 5)
#' p
#'
VisClusterFlowGraph <- function(object,
                                group.by = "seurat_clusters",
                                reduction = "umap",
                                point_size = 7,
                                point_alpha = 0.9,
                                label_size = 7) {
   # Seurat object check
   if (!inherits(object, "Seurat"))
      stop("object must be a Seurat object")

   # Groups in object@meta.data
   if (!(group.by %in% colnames(object@meta.data)))
      stop('group.by not found')

   # PCA or UMAP in object@reductions
   if (is.null(object@reductions[[reduction]])) {
      if (reduction == 'pca')
         object <- suppressMessages(Seurat::RunPCA(object))
      if (reduction == 'umap')
         object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(object, dims = 1:10)))
   }

   # Embeddings
   emb <- as.data.frame(Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE])
   colnames(emb) <- c('X1', 'X2')
   emb$group <- object@meta.data[[group.by]]
   cents <- emb |>
      dplyr::group_by(group) |>
      dplyr::summarize(x = mean(X1),
                       y = mean(X2),
                       .groups = 'drop')
   C <- as.matrix(cents[, c('x', 'y')])
   if (nrow(C) < 2) {
      return(
         ggplot2::ggplot(cents, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(size = 3) + ggplot2::geom_text(
            ggplot2::aes(label = group),
            vjust = -0.8,
            size = 3
         ) +
            svpp_theme() + ggplot2::labs(
               x = paste0(toupper(reduction), '1'),
               y = paste0(toupper(reduction), '2')
            )
      )
   }
   D <- as.matrix(stats::dist(C))

   # Prim's algorithm for MST
   k <- nrow(D)
   selected <- rep(FALSE, k)
   selected[1] <- TRUE
   edges <- data.frame(from = integer(), to = integer())
   while (sum(selected) < k) {
      best <- c(NA_integer_, NA_real_, NA_integer_)
      for (i in which(selected)) {
         for (j in which(!selected)) {
            w <- D[i, j]
            if (is.na(best[2]) || w < best[2])
               best <- c(i, w, j)
         }
      }
      edges <- rbind(edges, data.frame(from = best[1], to = best[3]))
      selected[best[3]] <- TRUE
   }
   edges$x0 <- cents$x[edges$from]
   edges$y0 <- cents$y[edges$from]
   edges$x1 <- cents$x[edges$to]
   edges$y1 <- cents$y[edges$to]

   # Plot
   ggplot2::ggplot() +
      ggplot2::geom_segment(
         data = edges,
         ggplot2::aes(
            x = x0,
            y = y0,
            xend = x1,
            yend = y1
         ),
         linewidth = 1,
         color = "#333333",
         arrow = ggplot2::arrow(length = ggplot2::unit(0.5, 'cm')),
      ) +
      ggplot2::geom_point(
         data = cents,
         ggplot2::aes(x = x, y = y, color = group),
         size = point_size,
         alpha = point_alpha
      ) +
      ggplot2::geom_text(
         data = cents,
         ggplot2::aes(
            x = x,
            y = y,
            label = paste0("Cluster", group),
            color = group
         ),
         vjust = -0.8,
         size = label_size
      ) +
      ggplot2::labs(x = paste0(toupper(reduction), '1'), y = paste0(toupper(reduction), '2')) +
      svpp_theme()
}
