# SeuratVisPro

Innovative visualization and analysis toolkit for Seurat v5. Focused on real-world workflows, it provides cluster stability diagnostics, batch mixing assessment, marker atlas, ligand–receptor directional scoring, expression trend curves, cluster tree with similarity heatmap, spatial overlay with non-spatial fallback, publication-ready dotmap, and a radial gene co-expression network. A bundled Shiny app (optional `bs4Dash`) enables interactive exploration.

## Rationale & Benefits

- Cluster stability: resolution × resampling agreement curves to select hyperparameters quantitatively — adding a missing stability metric to standard workflows.
- Batch mixing diagnostics: neighborhood batch difference proportion on embeddings to quickly locate integration issues.
- Markers & similarity: marker atlas plus cluster tree + similarity heatmap for multi-view understanding of group differences and proximity.
- Signaling & trends: LR co-expression heatmap across groups and pseudotime/group-based smoothed trends for mechanism and temporal analysis.
- Spatial fallback: if no images are available, overlay via `x,y` or embedding coordinates to keep workflows unblocked.
- Publication styling: unified `svpp_theme()` and `viridis` palettes for consistent, controllable visuals.

## Method Details

- Stability: `VisClusterStability` resamples per resolution, aligns cluster assignments by cell names, computes mean agreement, and plots curves.
- Batch: `VisBatchAlign` builds kNN on embeddings (distance-based), measures neighbor batch difference vs self, and compares across groups.
- Tree & heatmap: `VisClusterTree` uses group-average expression, supports Euclidean/correlation distances and multiple linkages, with a similarity heatmap.
- Signaling: `VisLigRec` computes LR source/target products from averages to form a directional score matrix.
- Trends: `VisGeneTrend` derives pseudotime from UMAP/PCA or uses metadata grouping; plots loess/GAM smoothed curves.
- Spatial: `VisSpatialOverlay` prefers `SpatialFeaturePlot`; falls back to `x,y` or embeddings otherwise.
- Dotmap & network: `VisRankedDotmap` combines mean expression and percent; `VisGeneCoexpHive` builds a radial network from correlations.

## 快速开始

```r
library(SeuratVisPro)
obj <- SeuratVisProExample(n_cells = 300, n_genes = 800, n_clusters = 3, spatial = TRUE)
VisClusterStability(obj)$plot
VisClusterTree(obj, dist.metric = 'correlation', show_heatmap = TRUE)
VisSpatialOverlay(obj, features = c('G1','G2'))
VisRankedDotmap(obj, group.by = 'seurat_clusters', top_n = 6)
VisGeneCoexpHive(obj, genes = paste0('G',1:12), threshold = 0.2)
```

## Full Function Examples

```r
library(SeuratVisPro)
obj <- SeuratVisProExample(n_cells=400, n_genes=800, n_clusters=4, spatial=TRUE)
obj$batch <- sample(c('A','B'), ncol(obj), replace=TRUE)

# 1) QC panel (static)
VisQCPanel(obj, genes_mt='^MT-', genes_ribo='^RPL|^RPS', group.by='seurat_clusters')

# 2) Cluster stability
stab <- VisClusterStability(obj, resolution_range=seq(0.2,1.2,0.2), reps=3)
stab$plot; head(stab$summary)

# 3) Marker atlas
ma <- VisMarkerAtlas(obj, top_n=5)
ma$plot; head(ma$markers)

# 4) Batch mixing diagnostics
ba <- VisBatchAlign(obj, batch='batch', reduction='pca')
ba$plot; head(ba$summary)

# 5) Gene trends (pseudotime or group)
VisGeneTrend(obj, features=c('G10','G20','G30'), by='pseudotime')

# 6) Ligand–receptor directional scores
lr <- data.frame(ligand=paste0('G',1:5), receptor=paste0('G',6:10))
lrres <- VisLigRec(obj, lr_table=lr, group.by='seurat_clusters')
lrres$plot; head(lrres$scores)

# 7) Module scores (meta features)
ms <- VisMetaFeature(obj, feature_sets=list(SetA=paste0('G',1:10), SetB=paste0('G',11:20)), group.by='seurat_clusters')
ms$plot

# 8) Cell cycle view
cc <- VisCellCycle(obj, s.genes=paste0('G',1:10), g2m.genes=paste0('G',11:20))
cc$plot

# 9) Embedding contours (innovative)
VisEmbeddingContour(obj, group.by='seurat_clusters', levels=5)

# 10) Ranked dotmap (innovative)
VisRankedDotmap(obj, group.by='seurat_clusters', top_n=6)

# 11) Radial co-expression network (innovative)
VisGeneCoexpHive(obj, genes=paste0('G',1:12), threshold=0.2)

# 12) Spatial overlay (with non-spatial fallback)
VisSpatialOverlay(obj, features=c('G1','G2','G3'))

# 13) Hex-binned entropy heatmap (innovative; group mixing)
VisHexEntropy(obj, group.by='seurat_clusters', bins=30)

# 14) Local Moran's I hotspots (innovative; spatial autocorrelation)
VisLocalMoran(obj, gene='G10', k=15)

# 15) Cluster centroid flow graph (innovative; MST transitions)
VisClusterFlowGraph(obj, group.by='seurat_clusters')

# 13) Launch Shiny (bs4Dash dashboard)
launchSeuratVisPro()  # requires bs4Dash installed
```

## Shiny

```r
# Install bs4Dash if needed
install.packages('bs4Dash')
launchSeuratVisPro() # uses bs4Dash dashboard layout
```

## pkgdown

With Pandoc available:

```r
pkgdown::build_site() # generates docs for GitHub Pages
```
Help pages for new functions are generated from roxygen comments and will appear under `docs/reference` after building the site.

## Dependencies

`Imports`: ggplot2, dplyr, tidyr, patchwork, shiny, DT, stats, methods, tibble, ggdendro, plotly, rlang, Seurat, scales, Matrix

Optional: bs4Dash (Shiny dashboard enhancements; will gracefully fallback if not installed).
