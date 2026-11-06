rm(list = ls())

PATH_CONFIG <- Sys.getenv("AML_PATH_CONFIG", unset = file.path("config", "paths.R"))
if (!file.exists(PATH_CONFIG)) {
  PATH_CONFIG <- file.path("config", "paths_example.R")
}
if (!file.exists(PATH_CONFIG)) {
  stop("Path configuration not found. Run scripts from the repository root or set AML_PATH_CONFIG.")
}
source(PATH_CONFIG, local = FALSE)

# 


 
 
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(parallel)
library(scCustomize)
library(doParallel)
library(qs)


ANNO_DIR <- SCRNA_ANNOTATION_DIR
hdWGCNA_DIR <- SCRNA_HDWGCNA_DIR
if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)
if (!dir.exists(hdWGCNA_DIR)) dir.create(hdWGCNA_DIR, recursive = TRUE)
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

.fig_id <- 0L
fig_id <- function() { .fig_id <<- .fig_id + 1L; .fig_id }
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}


setwd(ANNO_DIR)
seurat <- qread(SCRNA_HDWGCNA_INPUT_FILE)

seurat_obj <- seurat
setwd(hdWGCNA_DIR)

seurat_obj$celltype=Idents(seurat_obj)

seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.03,
  wgcna_name = "hdWGCNA"
)

head(seurat_obj@meta.data)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype", "type"),   # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony',   # select the dimensionality reduction to perform KNN on
  k = 15,   # nearest-neighbors parameter
  max_shared = 10,   # maximum number of shared cells between two metacells
  ident.group = 'type',   # set the Idents of the metacell seurat object
  min_cells = 50
)


seurat_obj <- NormalizeMetacells(seurat_obj)

cellName <- "AML_Blast"
#MonoBlast   AML_Blast   AML_Progenitor  Baso_Mast   ProMono   iMono   ncMono


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = cellName,   # the name of the group of interest in the group.by column
  group.by='celltype',   # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA',   # using RNA assay
  layer = 'data'   # using normalized data
)


seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'   # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(seurat_obj)
pdf(file="Fig01_power.pdf", width=8, height=7)
wrap_plots(plot_list, ncol=2)
dev.off()
power_table <- GetPowerTable(seurat_obj)
head(power_table)

# seurat_obj <- TestSoftPowers(
#   seurat_obj,
#   networkType = 'signed hybrid'
# )
# 
# plot_list <- PlotSoftPowers(seurat_obj)
# pdf(file="Fig01_power.pdf", width=8, height=7)
# wrap_plots(plot_list, ncol=2)
# dev.off()
# 
# power_table <- GetPowerTable(seurat_obj)
# head(power_table)

selectPower=9
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power = selectPower,
  tom_name = cellName,    # name of the topoligical overlap matrix written to disk
  deepSplit = 4,
  minModuleSize = 20,
  mergeCutHeight = 0.15
)
pdf(file="Fig02_cluster.pdf", width=8, height=6)
PlotDendrogram(seurat_obj, main=paste0(cellName, ' hdWGCNA Dendrogram'))
dev.off()


TOM <- GetTOM(seurat_obj)
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="celltype"
)
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'celltype', group_name = cellName
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(cellName, "-M")
)

pdf(file="Fig03_kME.pdf",width = 15,height = 15)
PlotKMEs(seurat_obj, ncol=5, n_hubs = 30)
dev.off()


modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
# show the first 6 columns:
head(modules[,1:6])
write.csv(modules, file="Tab01_moduleGene.csv", row.names=F)
hub_df <- GetHubGenes(seurat_obj, n_hubs = 30)
head(hub_df)
write.csv(hub_df, file="Tab02_hubGene.csv", row.names=F)
#saveRDS(seurat_obj, file='hdWGCNA_object.rds')


seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 50,
  method='Seurat'
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "tsne",
  features='hMEs',    # plot the hMEs
  order=TRUE     # order so the points with highest hMEs are on top
)
Ncol=5
Height=ifelse(length(mods)/5*3>6, length(mods)/5*3, 6)
pdf(file="Fig04_FeaturePlot.pdf",width = 15,height = Height)
wrap_plots(plot_list, ncol=Ncol)
dev.off()


seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$orig.ident), ' '))[,1]
pdf(file = "Fig05_ModuleRadarPlot.pdf", width = 8, height = 8)
ModuleRadarPlot(
  seurat_obj,
  group.by = 'cluster',
  barcodes = seurat_obj@meta.data %>% 
    subset(celltype == cellName) %>% 
    rownames(),
  axis.label.size = 4,
  grid.label.size = 4
)
dev.off()


pdf(file="Fig06_ModuleCorrelogram.pdf",width = 8,height = 7)
ModuleCorrelogram(seurat_obj)
dev.off()


library(dplyr)
library(tidyr)
library(openxlsx)
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
p6 <- DotPlot(seurat_obj, features = mods, group.by = "celltype")
dotplot_data <- p6$data
dotplot_wide_avg <- dotplot_data %>%
  select(features.plot, id, avg.exp) %>%
  pivot_wider(names_from = id, values_from = avg.exp)
dotplot_wide_pct <- dotplot_data %>%
  select(features.plot, id, pct.exp) %>%
  pivot_wider(names_from = id, values_from = pct.exp)
write.xlsx(
  list(
    "AverageExpression" = dotplot_wide_avg,
    "PercentExpression" = dotplot_wide_pct
  ),
  file = "Tab03__cellModuleCor_wide.xlsx",
  rowNames = FALSE
)
p6 <- p6 + RotatedAxis() +     
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")    
pdf(file="Fig07_cellModuleCor.pdf", width = 9, height = 6)
print(p6)
dev.off()


p7 <- VlnPlot(
  seurat_obj,
  features = mods,
  group.by = "celltype",
  ncol=Ncol,
  sort="increasing",
  pt.size = 0
)
p7 <- p7 + geom_boxplot(width = .25, fill = "white")
p7 <- p7 + xlab("") + ylab("hME") + NoLegend()
pdf(file="Fig08_VlnPlot.pdf", width =16, height =Height)
print(p7)
dev.off()


ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "all")


library(future)
options(future.globals.maxSize = 20 * 1024^3)  
pdf(file="Fig09_hubNetwork.pdf", width =7, height =7)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3,     # The number of hub genes to plot for each module
  n_other = 5,     # The number of non-hub genes to sample from each module
  edge_prop = 0.75,    # The proportion of edges in the graph to sample
  mods = "all"
)
dev.off()


seurat_obj <- RunModuleUMAP(
  seurat_obj, 
  n_hubs = 10,     #number of hub genes to use in the UMAP computation
  n_neighbors = 15,   #The size of local neighborhood used for manifold approximation
  min_dist = 0.1    #The effective minimum distance between embedded points
)
umap_df <- GetModuleUMAP(seurat_obj)

pdf(file = "Fig10_geneUmap.pdf", width = 9, height = 7)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha = 0.25,
  sample_edges = TRUE,
  edge_prop = 0.1,
  label_hubs = 2,
  keep_grey_edges = FALSE
)
dev.off()
