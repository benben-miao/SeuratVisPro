# SeuratVisPro

面向Seurat v5的创新型单细胞可视化与分析工具包。聚焦真实项目痛点，提供聚类稳健性诊断、批次混合评估、群体标志物图谱、配体-受体方向评分、趋势曲线、群体树与相似性热图、空间叠色回退、发表级点图以及基因协同结构径向网络等功能，并集成Shiny（支持bs4Dash可选）便于交互探索。

## 设计理由与优势

- 聚类稳健性：分辨率×重采样一致性曲线定量选择超参，补足官方工作流的“稳健性度量”。
- 批次混合诊断：嵌入空间邻域批次差异比例快速定位整合问题群体。
- 标志物与相似性：Marker atlas与群体树+热图组合，多视角呈现群体差异与近似程度。
- 通信与趋势：群体×群体LR共表达热图与伪时间/分组平滑曲线，便于机制与时序分析。
- 空间回退：无真实图像也可依赖坐标/嵌入进行“空间式”叠色，保障工作流不中断。
- 发表级美化：统一`svpp_theme()`与`viridis`色系参数化，图形风格一致可控。

## 方法细节

- 稳健性：`VisClusterStability` 对每个分辨率进行重复抽样，按细胞名对齐完整与子样本的簇分配，计算一致性均值并绘制曲线。
- 批次：`VisBatchAlign` 在嵌入上构建k近邻（基于距离排序），统计邻居批次与自身差异比例，分组比较混合程度。
- 树与热图：`VisClusterTree` 基于群体平均表达，支持欧氏或相关性度量与多种联接法，树下叠加相似性热图。
- 通信：`VisLigRec` 基于平均表达计算LR对在源/靶群体的乘积均值，形成方向性评分矩阵。
- 趋势：`VisGeneTrend` 从UMAP/PCA近似伪时间或任意元数据分组，绘制loess/gam平滑曲线。
- 空间：`VisSpatialOverlay` 优先使用`SpatialFeaturePlot`，无图像时回退至`x,y`或嵌入坐标点染色。
- 点图与网络：`VisRankedDotmap` 合并平均表达与表达比例，`VisGeneCoexpHive` 使用相关性构建径向网络。

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

## 全函数示例

```r
library(SeuratVisPro)
obj <- SeuratVisProExample(n_cells=400, n_genes=800, n_clusters=4, spatial=TRUE)
obj$batch <- sample(c('A','B'), ncol(obj), replace=TRUE)

# 1) 质控面板（静态）
VisQCPanel(obj, genes_mt='^MT-', genes_ribo='^RPL|^RPS', group.by='seurat_clusters')

# 2) 聚类稳健性
stab <- VisClusterStability(obj, resolution_range=seq(0.2,1.2,0.2), reps=3)
stab$plot; head(stab$summary)

# 3) 标志物图谱
ma <- VisMarkerAtlas(obj, top_n=5)
ma$plot; head(ma$markers)

# 4) 批次混合诊断
ba <- VisBatchAlign(obj, batch='batch', reduction='pca')
ba$plot; head(ba$summary)

# 5) 基因趋势（伪时间或分组）
VisGeneTrend(obj, features=c('G10','G20','G30'), by='pseudotime')

# 6) 配体-受体方向评分
lr <- data.frame(ligand=paste0('G',1:5), receptor=paste0('G',6:10))
lrres <- VisLigRec(obj, lr_table=lr, group.by='seurat_clusters')
lrres$plot; head(lrres$scores)

# 7) 模块分数（元特征）
ms <- VisMetaFeature(obj, feature_sets=list(SetA=paste0('G',1:10), SetB=paste0('G',11:20)), group.by='seurat_clusters')
ms$plot

# 8) 细胞周期视图
cc <- VisCellCycle(obj, s.genes=paste0('G',1:10), g2m.genes=paste0('G',11:20))
cc$plot

# 9) 嵌入等高线（原创）
VisEmbeddingContour(obj, group.by='seurat_clusters', levels=5)

# 10) 排序点图（原创）
VisRankedDotmap(obj, group.by='seurat_clusters', top_n=6)

# 11) 协同网络径向图（原创）
VisGeneCoexpHive(obj, genes=paste0('G',1:12), threshold=0.2)

# 12) 空间叠色（含非空间回退）
VisSpatialOverlay(obj, features=c('G1','G2','G3'))

# 13) Hex-binned 熵热图（原创，群体混合度）
VisHexEntropy(obj, group.by='seurat_clusters', bins=30)

# 14) 局部Moran's I热点（原创，空间自相关）
VisLocalMoran(obj, gene='G10', k=15)

# 15) 群体质心流图（原创，MST指示潜在转移）
VisClusterFlowGraph(obj, group.by='seurat_clusters')

# 13) 启动Shiny（bs4Dash仪表盘）
launchSeuratVisPro()  # 需已安装 bs4Dash
```

## Shiny

```r
# 需要安装 bs4Dash
install.packages('bs4Dash')
launchSeuratVisPro() # 使用bs4Dash仪表盘布局
```

## pkgdown

在具备Pandoc的环境中：

```r
pkgdown::build_site() # 生成docs目录用于GitHub Pages
```
新函数的帮助文档由roxygen注释自动生成，构建站点后将出现在`docs/reference`。

## 依赖

`Imports`: ggplot2, dplyr, tidyr, patchwork, shiny, DT, stats, methods, tibble, ggdendro, plotly, rlang, Seurat, scales, Matrix

可选：bs4Dash（Shiny界面增强，未安装则自动回退）。
