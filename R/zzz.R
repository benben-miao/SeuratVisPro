#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("SeuratVisPro 0.1.0 - Innovative Seurat v5 visualization toolkit")
}

utils::globalVariables(c(
  "Phase","agreement","avg_exp","cluster","gene","group","ligand","mix_prop",
  "nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","receptor","score","target",
  "x","xend","y","yend","X1","X2","x0","y0","x1","y1","pct","sim","w","cluster1","cluster2","level","r",
  "entropy","bx","by","I","n","resolution","..level..","expr"
))
