rm(list = ls())

PATH_CONFIG <- Sys.getenv("AML_PATH_CONFIG", unset = file.path("config", "paths.R"))
if (!file.exists(PATH_CONFIG)) {
  PATH_CONFIG <- file.path("config", "paths_example.R")
}
if (!file.exists(PATH_CONFIG)) {
  stop("Path configuration not found. Run scripts from the repository root or set AML_PATH_CONFIG.")
}
source(PATH_CONFIG, local = FALSE)

library(qs)
library(Seurat)
library(ggplot2)
library(data.table)
library(SeuratDisk)
library(decontX)
library(celda)
library(DoubletFinder)
library(harmony)
library(ggsci)
library(clustree)

############################################################-
############################################################-

AML_TXT_DIR <- SCRNA_AML_TXT_DIR
Normal_TXT_DIR <- SCRNA_NORMAL_TXT_DIR

OUT_QC_DIR <- SCRNA_QC_DIR

dir.create(OUT_QC_DIR, recursive = TRUE, showWarnings = FALSE)

randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

.fig_id <- 0L
fig_id <- function() {.fig_id <<- .fig_id + 1L; .fig_id}
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}

############################################################-
############################################################-
# setwd(RAW_10X_TREAT)
#Seurat_R
#----treat
#--------P1
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#--------P2
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#--------P3
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#----normal
#--------N1
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#--------N2
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#--------N3
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz


read_txt_to_seurat <- function(file) {
  message("📥 Reading: ", basename(file))
  dat <- fread(file, data.table = FALSE)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat <- as.matrix(dat)
  storage.mode(dat) <- "integer"
  
  seu <- CreateSeuratObject(
    counts = dat,
    project = "AML",
    min.cells = 3,
    min.features = 200
  )
  
  sid <- tools::file_path_sans_ext(basename(file))
  seu$type <- sid
  seu$orig.ident <- sid
  return(seu)
}

files <- list.files(AML_TXT_DIR, pattern = "\\.txt$", full.names = TRUE)
seurat_list <- lapply(files, read_txt_to_seurat)
names(seurat_list) <- tools::file_path_sans_ext(basename(files))


seurat_aml <- Reduce(function(x, y) merge(x, y), seurat_list)

head(seurat_aml@meta.data)
dim(seurat_aml)


seurat_aml[["percent.mt"]] <- PercentageFeatureSet(seurat_aml, pattern = "^MT-")
seurat_aml[["percent.rb"]] <- PercentageFeatureSet(seurat_aml, pattern = "^RP")
seurat_aml[["percent.hb"]] <- PercentageFeatureSet(seurat_aml, pattern = "^HB[AB]")

randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}
randomColors <- replicate(100,randomColor())

p1 <- VlnPlot(
  seurat_aml, 
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb", "percent.hb"), 
  ncol = 5
) + scale_fill_manual(values = randomColors)

pdf(file.path(OUT_QC_DIR, "Fig01_QC_violin_AML.pdf"), width = 20, height = 8)
print(p1)
dev.off()


min_nCount <- 0
max_nCount <- 2500
min_nFeature <- 0
max_nFeature <- 1000
thr_mt <- 18
min_thr_rb <- 5
max_thr_rb <- 45
thr_hb <- 18

dim_before <- dim(seurat_aml)

mask1 <- min_nCount <= seurat_aml$nCount_RNA & seurat_aml$nCount_RNA <= max_nCount
mask2 <- min_nFeature <= seurat_aml$nFeature_RNA &  seurat_aml$nFeature_RNA <= max_nFeature
mask3 <- seurat_aml$percent.mt <= thr_mt
mask4 <- min_thr_rb <= seurat_aml$percent.rb & seurat_aml$percent.rb <= max_thr_rb
# mask5 <- seurat_pnh$percent.mt <= thr_hb
seurat_aml <- seurat_aml[, mask1 & mask2 & mask3]

dim_after <- dim(seurat_aml)

log_file <- file.path(OUT_QC_DIR, "Fig01_QC_pnh_log.txt")

cat(
  " 时间:\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "阈值:\n",
  "  nCount_RNA: ", min_nCount, " - ", max_nCount, "\n",
  "  nFeature_RNA: ", min_nFeature, " - ", max_nFeature, "\n",
  "  percent.mt <= ", thr_mt, "\n",
  if (exists("thr_rb")) {
    paste0("  percent.rb <= ", thr_rb, "\n")
  } else {
    "  percent.rb: 未筛选\n"
  },
  if (exists("thr_hb")) {
    paste0("  percent.hb <= ", thr_hb, "\n")
  } else {
    "  percent.hb: 未筛选\n"
  },
  "\n",
  "筛选前 dim: ", dim_before[1], " features × ", dim_before[2], " cells\n",
  "筛选后 dim: ", dim_after[1], " features × ", dim_after[2], " cells\n",
  "说明: dim 的第一个数字是基因/特征数，第二个是细胞数\n",
  "----------------------------------------\n",
  file = log_file, append = TRUE
)

qsave(seurat_aml, file.path(OUT_QC_DIR, "1.QC_AML_filtered.qs"))


setwd(Normal_TXT_DIR)

read_txt_to_seurat <- function(file) {
  message("📥 Reading: ", basename(file))
  dat <- fread(file, data.table = FALSE)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat <- as.matrix(dat)
  storage.mode(dat) <- "integer"
  
  seu <- CreateSeuratObject(
    counts = dat,
    project = "Normal",
    min.cells = 3,
    min.features = 200
  )
  
  sid <- tools::file_path_sans_ext(basename(file))
  seu$type <- sid
  seu$orig.ident <- sid
  return(seu)
}

files <- list.files(Normal_TXT_DIR, pattern = "\\.txt$", full.names = TRUE)
seurat_list <- lapply(files, read_txt_to_seurat)
names(seurat_list) <- tools::file_path_sans_ext(basename(files))


seurat_con <- Reduce(function(x, y) merge(x, y), seurat_list)

head(seurat_con@meta.data)
dim(seurat_con)

seurat_con[["percent.mt"]] <- PercentageFeatureSet(seurat_con, pattern = "^MT-")
seurat_con[["percent.rb"]] <- PercentageFeatureSet(seurat_con, pattern = "^RP")
seurat_con[["percent.hb"]] <- PercentageFeatureSet(seurat_con, pattern = "^HB[AB]")

randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}
randomColors <- replicate(100,randomColor())

p2 <- VlnPlot(
  seurat_con,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb", "percent.hb"),
  ncol = 5
) + scale_fill_manual(values = randomColors)

pdf(file.path(OUT_QC_DIR, "Fig02_QC_violin_col.pdf"), width = 20, height = 8)
print(p2)
dev.off()


min_nCount <- 0
max_nCount <- 1000
min_nFeature <- 0
max_nFeature <- 400
thr_mt <- 8

dim_before <- dim(seurat_con)

mask1 <- min_nCount <= seurat_con$nCount_RNA & seurat_con$nCount_RNA <= max_nCount
mask2 <- min_nFeature <= seurat_con$nFeature_RNA &  seurat_con$nFeature_RNA <= max_nFeature
mask3 <- seurat_con$percent.mt <= thr_mt
# mask4 <- seurat_con$nFeature_RNA <= thr_rb
# mask5 <- seurat_con$percent.mt <= thr_hb
seurat_con <- seurat_con[, mask1 & mask2 & mask3]

dim_after <- dim(seurat_con)

log_file <- file.path(OUT_QC_DIR, "Fig02_QC_col_log.txt")

cat(
  " 时间:\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "阈值:\n",
  "  nCount_RNA: ", min_nCount, " - ", max_nCount, "\n",
  "  nFeature_RNA: ", min_nFeature, " - ", max_nFeature, "\n",
  "  percent.mt <= ", thr_mt, "\n",
  if (exists("thr_rb")) {
    paste0("  percent.rb <= ", thr_rb, "\n")
  } else {
    "  percent.rb: 未筛选\n"
  },
  if (exists("thr_hb")) {
    paste0("  percent.hb <= ", thr_hb, "\n")
  } else {
    "  percent.hb: 未筛选\n"
  },
  "\n",
  "筛选前 dim: ", dim_before[1], " features × ", dim_before[2], " cells\n",
  "筛选后 dim: ", dim_after[1], " features × ", dim_after[2], " cells\n",
  "说明: dim 的第一个数字是基因/特征数，第二个是细胞数\n",
  "----------------------------------------\n",
  file = log_file, append = TRUE
)

qsave(seurat_con, file.path(OUT_QC_DIR, "2.QC_COL_filtered.qs"))

seurat_obj <- merge(seurat_con, seurat_aml, add.cell.ids = c("Control", "AML"))

head(seurat_obj@meta.data)
dim(seurat_obj)


pdf(file.path(OUT_QC_DIR, "Fig03_Con_AML_VlnPlot_QC.pdf"))
VlnPlot(seurat_obj, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.hb"),ncol = 5)+
  scale_fill_manual(values =randomColors)
dev.off()

LayerData(seurat_obj, assay = "RNA", layer = "counts")
seurat_obj <- JoinLayers(seurat_obj)
dim(seurat_obj[["RNA"]]$counts )

qsave(seurat_obj, file.path(OUT_QC_DIR, "3.Con_AML_Seurat_obj.qs"))


############################################################-
############################################################-

# Normalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)

plot1 <- VariableFeaturePlot(seurat_obj) + ggtitle(fig_title("Variable features (all)"))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + ggtitle(fig_title("Top variable features (top10)"))
save_pdf(file.path(OUT_QC_DIR, "Fig04_Variable_features.pdf"), plot2, w = 7, h = 6)

# # Scale
# all.genes <- rownames(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj, features = all.genes)

seurat_obj <- ScaleData(
  seurat_obj,
  features = VariableFeatures(seurat_obj)
)


# PCA / UMAP
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

p_umap <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.5) +
  theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5)) +
  ggtitle(fig_title("UMAP (post-normalization)"))
save_pdf(file.path(OUT_QC_DIR, "Fig05_UMAP_post_normalization.pdf"), p_umap, w = 7.5, h = 5.5)

p_pca <- DimPlot(seurat_obj, reduction = "pca", pt.size = 0.5) +
  theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5)) +
  ggtitle(fig_title("PCA (post-normalization)"))
save_pdf(file.path(OUT_QC_DIR, "Fig06_PCA_post_normalization.pdf"), p_pca, w = 7.5, h = 5.5)


qsave(seurat_obj, file.path(OUT_QC_DIR, "4.Normalized_scaled.qs"))
############################################################
############################################################
suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(celda)
  library(ggplot2)
})


seurat_obj_tmp <- qread(file.path(OUT_QC_DIR, "4.Normalized_scaled.qs"))

seurat_obj_tmp <- qread(file.path(OUT_QC_DIR, "7.Doublet_removed.qs"))
DefaultAssay(seurat_obj_tmp) <- "RNA"
# counts <- GetAssayData(
#   object = seurat_obj_tmp,
#   assay  = "RNA",
#   layer  = "counts"
# )
# 
# 
# message(">>> Running decontX (ambient RNA decontamination) ...")
# contam <- dxc$contamination
# seurat_obj_tmp$contamination <- contam
# 
# p_contam_before <- FeaturePlot(
#   seurat_obj_tmp, features = "contamination", raster = FALSE
# ) + ggtitle(fig_title("Ambient contamination (before filtering)")) +
#   theme_bw() + 
#   xlab("UMAP_1") + 
#   ylab("UMAP_2") +
# 
# save_pdf(file.path(OUT_QC_DIR, "Fig07_Contamination_before.pdf"), p_contam_before, w = 8, h = 6)
# 
# # thr_contam <- 0.2 
# # 
# # sce_filt <- seurat_obj_tmp[, seurat_obj_tmp$contamination < thr_contam]
# 
# 
# thr_contam <- quantile(seurat_obj_tmp$contamination, 0.95)
# 
# sce_filt <- subset(
#   seurat_obj_tmp,
#   subset =
#     contamination <= thr_contam &
#     nCount_RNA >= 500 & nFeature_RNA >= 200 & percent.mt <= 10
# )
# 
# 
# p_contam_after <- FeaturePlot(
#   sce_filt, features = "contamination", raster = FALSE
# ) + ggtitle(fig_title("Ambient contamination (after filtering)")) +
#   theme_bw() + 
#   xlab("UMAP_1") + 
#   ylab("UMAP_2") +
# 
# save_pdf(file.path(OUT_QC_DIR, "Fig08_Contamination_after.pdf"), p_contam_after, w = 8, h = 6)
# 
# dim(sce_filt)
# qsave(sce_filt, file.path(OUT_QC_DIR, "5.Decontaminated.qs"))
# 
# ----------------------------
# ----------------------------
rna_assay <- seurat_obj_tmp[["RNA"]]
all_layers <- Layers(rna_assay)
count_layers <- grep("^counts", all_layers, value = TRUE)

if (length(count_layers) == 0) {
  stop("❌ 没有找到 counts* layers。请先检查：Layers(seurat_obj_tmp[['RNA']])")
}
message(">>> Found count layers: ", paste(count_layers, collapse = ", "))

seurat_obj_tmp$contamination <- NA_real_

# ----------------------------
# ----------------------------
has_cluster <- "seurat_clusters" %in% colnames(seurat_obj_tmp@meta.data)

for (ly in count_layers) {
  message("\n>>> Running decontX on layer: ", ly)
  
  mat <- LayerData(seurat_obj_tmp, assay = "RNA", layer = ly)
  
  if (is.null(mat) || ncol(mat) == 0 || nrow(mat) == 0) {
    message("... skip empty layer: ", ly)
    next
  }
  
  cells <- colnames(mat)
  
  z <- NULL
  if (has_cluster) {
    z <- seurat_obj_tmp@meta.data[cells, "seurat_clusters", drop = TRUE]
    z <- as.factor(z)
  }
  
  sce <- SingleCellExperiment(list(counts = mat))
  
  sce <- celda::decontX(sce, z = z, verbose = TRUE)
  
  contam_vec <- NULL
  if ("decontX_contamination" %in% colnames(colData(sce))) {
    contam_vec <- colData(sce)$decontX_contamination
  } else if ("contamination" %in% colnames(colData(sce))) {
    contam_vec <- colData(sce)$contamination
  } else {
    stop("❌ decontX 输出里找不到 contamination 字段（请检查 colnames(colData(sce))）")
  }
  
  seurat_obj_tmp$contamination[cells] <- as.numeric(contam_vec)
  
  rm(sce, mat, contam_vec); gc()
}

if (all(is.na(seurat_obj_tmp$contamination))) {
  stop("❌ contamination 全是 NA：说明 decontX 没有成功写回。请检查 layer 与细胞名是否匹配。")
}

message(">>> decontX finished. contamination summary:")
print(summary(seurat_obj_tmp$contamination))

# ----------------------------
# ----------------------------
p_contam_before <- FeaturePlot(
  seurat_obj_tmp,
  features = "contamination",
  reduction = "umap",
  raster = TRUE
) +
  ggtitle("Ambient contamination (before filtering)") +
  theme_bw() +
  xlab("UMAP_1") + ylab("UMAP_2") +
  scale_color_gradient(low = "lightblue", high = "darkblue", limits = c(0, 1))

ggsave(
  filename = file.path(OUT_QC_DIR, "Fig07_Contamination_before.pdf"),
  plot = p_contam_before,
  width = 8, height = 6
)

# ----------------------------
# ----------------------------

thr_contam <- 0.2

md <- seurat_obj_tmp@meta.data

cells_keep <- rownames(md)[
  !is.na(md$contamination) &
    md$contamination <= thr_contam &
    md$nCount_RNA >= 200 &
    md$nFeature_RNA >= 100 
   
]

sce_filt <- subset(seurat_obj_tmp, cells = cells_keep)

message(">>> Filtered dims: ")
print(dim(sce_filt))


# ----------------------------
# ----------------------------
p_contam_after <- FeaturePlot(
  sce_filt,
  features = "contamination",
  reduction = "umap",
  raster = TRUE
) +
  ggtitle("Ambient contamination (after filtering)") +
  theme_bw() +
  xlab("UMAP_1") + ylab("UMAP_2") +
  scale_color_gradient(low = "lightblue", high = "darkblue", limits = c(0, 1))

ggsave(
  filename = file.path(OUT_QC_DIR, "Fig08_Contamination_after.pdf"),
  plot = p_contam_after,
  width = 8, height = 6
)

# ----------------------------
# ----------------------------
qsave(sce_filt, file.path(OUT_QC_DIR, "5.Decontaminated.qs"))
message("✅ Saved: ", file.path(OUT_QC_DIR, "5.Decontaminated.qs"))


############################################################-
############################################################-

seurat_cc <- qread(file.path(OUT_QC_DIR, "5.Decontaminated.qs"))

# s.genes   <- intersect(cc.genes$s.genes, rownames(seurat_cc))
# g2m.genes <- intersect(cc.genes$g2m.genes, rownames(seurat_cc))
# 
# seurat_cc <- CellCycleScoring(seurat_cc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# # S/G2M violin
# p_cc_s  <- VlnPlot(seurat_cc, features = "S.Score", pt.size = 0, group.by = 'type') + NoLegend() +
#   ggtitle(fig_title("Cell cycle S.Score (by Type)"))
# p_cc_g2m <- VlnPlot(seurat_cc, features = "G2M.Score", pt.size = 0, group.by = 'type') + NoLegend() +
#   ggtitle(fig_title("Cell cycle G2M.Score (by Type)"))
# save_pdf(file.path(OUT_QC_DIR, "Fig09_CellCycle_S_and_G2M.pdf"), p_cc_s + p_cc_g2m, w = 10, h = 5)
# 
# seurat_cc <- ScaleData(seurat_cc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_cc))
# p_umap_before_cc <- DimPlot(seurat_cc, reduction = "umap", group.by = "type") +
#   ggtitle(fig_title("UMAP before re-run (with CC regressed)"))
# seurat_cc_new <- RunPCA(seurat_cc, features = VariableFeatures(seurat_cc))
# seurat_cc_new <- RunUMAP(seurat_cc_new, dims = 1:15, reduction = "pca")
# p_umap_after_cc <- DimPlot(seurat_cc_new, reduction = "umap", group.by = "type") +
#   ggtitle(fig_title("UMAP after regressing cell cycle"))
# save_pdf(file.path(OUT_QC_DIR, "Fig10_CC_effect_before_after.pdf"), p_umap_before_cc + p_umap_after_cc, w = 10, h = 5)
# print(dim(seurat_cc_new))
# 
# 
# qsave(seurat_cc_new, file.path(OUT_QC_DIR, "6.CC_regressed.qs"))


if (length(VariableFeatures(seurat_cc)) == 0) {
  seurat_cc <- FindVariableFeatures(seurat_cc, selection.method = "vst", nfeatures = 2000)
}

s.genes   <- intersect(cc.genes$s.genes, rownames(seurat_cc))
g2m.genes <- intersect(cc.genes$g2m.genes, rownames(seurat_cc))

seurat_cc <- CellCycleScoring(
  seurat_cc,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# S/G2M violin
p_cc_s <- VlnPlot(seurat_cc, features = "S.Score", pt.size = 0, group.by = "type") +
  NoLegend() + ggtitle(fig_title("Cell cycle S.Score (by Type)"))
p_cc_g2m <- VlnPlot(seurat_cc, features = "G2M.Score", pt.size = 0, group.by = "type") +
  NoLegend() + ggtitle(fig_title("Cell cycle G2M.Score (by Type)"))

save_pdf(
  file.path(OUT_QC_DIR, "Fig09_CellCycle_S_and_G2M.pdf"),
  p_cc_s + p_cc_g2m,
  w = 10, h = 5
)

seurat_cc <- ScaleData(
  seurat_cc,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = VariableFeatures(seurat_cc)
)

p_umap_before_cc <- DimPlot(seurat_cc, reduction = "umap", group.by = "type") +
  ggtitle(fig_title("UMAP before re-run (with CC regressed)"))

seurat_cc_new <- RunPCA(
  seurat_cc,
  features = VariableFeatures(seurat_cc)
)

seurat_cc_new <- RunUMAP(
  seurat_cc_new,
  dims = 1:15,
  reduction = "pca"
)

p_umap_after_cc <- DimPlot(seurat_cc_new, reduction = "umap", group.by = "type") +
  ggtitle(fig_title("UMAP after regressing cell cycle"))

save_pdf(
  file.path(OUT_QC_DIR, "Fig10_CC_effect_before_after.pdf"),
  p_umap_before_cc + p_umap_after_cc,
  w = 10, h = 5
)


print(dim(seurat_cc_new))

qsave(seurat_cc_new, file.path(OUT_QC_DIR, "6.CC_regressed.qs"))


############################################################-
############################################################-
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)

# seurat_df <- qread(file.path(OUT_QC_DIR, "6.CC_regressed.qs"))

# seurat_df <- qread(file.path(OUT_QC_DIR, "4.Normalized_scaled.qs"))

head(seurat_df@meta.data)
dim(seurat_df)

packageVersion("Seurat")
names(seurat_df@reductions)
table(seurat_df$type)

############################################################
############################################################

stopifnot(exists("seurat_df"))
stopifnot(exists("OUT_QC_DIR"))
stopifnot("type" %in% colnames(seurat_df@meta.data))
stopifnot(all(c("pca","umap") %in% names(seurat_df@reductions)))

DefaultAssay(seurat_df) <- "RNA"

# ----------------------------------------------------------
# ----------------------------------------------------------
pca_embed_all  <- Embeddings(seurat_df, "pca")
umap_embed_all <- Embeddings(seurat_df, "umap")


# ----------------------------------------------------------
# ----------------------------------------------------------
seu_list <- SplitObject(seurat_df, split.by = "type")

rm(seurat_df); gc()

# ----------------------------------------------------------
# ----------------------------------------------------------
run_scdblfinder_group <- function(seu, nfeatures = 2000, verbose = TRUE) {
  log_msg <- function(...) if (verbose) cat(..., "\n")
  
  seu@graphs <- list()
  seu@neighbors <- list()
  seu@reductions <- list()
  
  log_msg("🔄 as.SingleCellExperiment")
  sce <- as.SingleCellExperiment(seu)
  
  if ("counts" %in% assayNames(sce)) {
    assays(sce) <- assays(sce)["counts"]
  }
  
  log_msg("🧪 scDblFinder")
  sce <- scDblFinder(sce, nfeatures = nfeatures)
  
  seu$scDblFinder.class <- colData(sce)$scDblFinder.class
  seu$scDblFinder.score <- colData(sce)$scDblFinder.score
  
  rm(sce); gc()
  
  keep <- colnames(seu)[seu$scDblFinder.class == "singlet"]
  seu  <- subset(seu, cells = keep)
  
  seu
}

failed_groups <- character()
singlet_list <- list()
dbl_stat <- data.frame(
  type = names(seu_list),
  n_before = as.integer(table(names(seu_list))[names(seu_list)]),
  stringsAsFactors = FALSE
)

dbl_stat$n_before <- as.integer(sapply(seu_list, ncol))

for (nm in names(seu_list)) {
  cat("\n=== 🔁 scDblFinder:", nm, "(", ncol(seu_list[[nm]]), "cells ) ===\n")
  
  res <- tryCatch(
    run_scdblfinder_group(seu_list[[nm]], nfeatures = 2000, verbose = TRUE),
    error = function(e) { message("❗ ", nm, " failed: ", conditionMessage(e)); NULL }
  )
  
  seu_list[[nm]] <- NULL
  gc()
  
  if (is.null(res)) {
    failed_groups <- c(failed_groups, nm)
  } else {
    singlet_list[[nm]] <- res
  }
}

cat("\n✅ scDblFinder 完成。\n")
if (length(failed_groups) > 0) {
  cat("❌ 失败分组：\n")
  print(failed_groups)
}

stopifnot(length(singlet_list) > 0)

dbl_stat$n_after <- as.integer(sapply(dbl_stat$type, function(x) {
  if (x %in% names(singlet_list)) ncol(singlet_list[[x]]) else NA_integer_
}))

dbl_stat$kept_ratio <- dbl_stat$n_after / dbl_stat$n_before
print(dbl_stat)

# ----------------------------------------------------------
# ----------------------------------------------------------
singlet_list <- lapply(singlet_list, function(x) {
  x@graphs <- list()
  x@neighbors <- list()
  x
})
gc()

seurat_singlet <- singlet_list[[1]]
if (length(singlet_list) > 1) {
  seurat_singlet <- merge(seurat_singlet, y = singlet_list[-1], merge.data = FALSE)
}
rm(singlet_list); gc()

DefaultAssay(seurat_singlet) <- "RNA"

# ----------------------------------------------------------
# ----------------------------------------------------------
common_cells <- intersect(colnames(seurat_singlet), rownames(pca_embed_all))
if (length(common_cells) == 0) stop("❌ 合并后的对象与原 PCA cells 无交集，检查细胞名是否中途被改动。")

if (length(common_cells) < ncol(seurat_singlet)) {
  seurat_singlet <- subset(seurat_singlet, cells = common_cells)
}

pca_embed  <- pca_embed_all[colnames(seurat_singlet), , drop = FALSE]
umap_embed <- umap_embed_all[colnames(seurat_singlet), , drop = FALSE]

seurat_singlet[["pca"]] <- CreateDimReducObject(
  embeddings = pca_embed,
  key = "PC_",
  assay = DefaultAssay(seurat_singlet)
)
seurat_singlet[["umap"]] <- CreateDimReducObject(
  embeddings = umap_embed,
  key = "UMAP_",
  assay = DefaultAssay(seurat_singlet)
)

rm(pca_embed_all, umap_embed_all, pca_embed, umap_embed); gc()

# ----------------------------------------------------------
# ----------------------------------------------------------
#df_colors <- c("doublet" = "#F8766D", "singlet" = "#00BFC4")
p_after <- DimPlot(
  seurat_singlet,
  reduction = "umap",
  group.by = "scDblFinder.class",
  raster = TRUE
) +
  scale_color_manual(values = df_colors, drop = FALSE) +
  ggtitle("Doublet classification (after filtering)")

ggsave(
  filename = file.path(OUT_QC_DIR, "Fig11_Doublet_after_filter.pdf"),
  plot = p_after,
  width = 7.5,
  height = 5.5
)


p_type <- DimPlot(seurat_singlet, reduction = "umap", group.by = "type") +
  ggtitle("UMAP after scDblFinder (singlets only) - colored by type")
ggsave(file.path(OUT_QC_DIR, "Fig11_UMAP_singlets_by_type.pdf"), p_type, width = 7.5, height = 5.5)

p_score <- FeaturePlot(seurat_singlet, reduction = "umap", features = "scDblFinder.score") +
  ggtitle("UMAP after scDblFinder (singlets only) - scDblFinder.score")
ggsave(file.path(OUT_QC_DIR, "Fig12_UMAP_singlets_score.pdf"), p_score, width = 7.5, height = 5.5)

cat("\n✅ 合并后 singlet 对象维度：\n")
print(dim(seurat_singlet))

# ----------------------------------------------------------
# ----------------------------------------------------------
qsave(seurat_singlet, file.path(OUT_QC_DIR, "7.Doublet_removed.qs"))

############################################################-
############################################################-

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Matrix)
  library(harmony)
  library(ggplot2)
  library(future)
})

suppressWarnings({
  plan("multisession", workers = 6)
})
options(future.globals.maxSize = 60 * 1024^3)

# ----------------------------
# ----------------------------
seurat_hm_raw <- qread(file.path(OUT_QC_DIR, "6.CC_regressed.qs"))
DefaultAssay(seurat_hm_raw) <- "RNA"
stopifnot("orig.ident" %in% colnames(seurat_hm_raw@meta.data))

# ----------------------------
# ----------------------------
count_layers <- grep("^counts\\.", Layers(seurat_hm_raw[["RNA"]]), value = TRUE)

if (length(count_layers) > 0) {
  message(">>> Detected multi-layer RNA assay: merging ", length(count_layers), " counts layers...")
  
  counts_list <- lapply(count_layers, function(ly) {
    GetAssayData(seurat_hm_raw, assay = "RNA", layer = ly)
  })
  
  counts_merged <- Reduce(base::cbind, counts_list)
  rm(counts_list); gc()
  
  md <- seurat_hm_raw@meta.data[colnames(counts_merged), , drop = FALSE]
  seurat_hm <- CreateSeuratObject(counts = counts_merged, meta.data = md)
  DefaultAssay(seurat_hm) <- "RNA"
  
  rm(counts_merged, md, seurat_hm_raw); gc()
} else {
  message(">>> RNA assay already single-layer. Proceeding directly.")
  seurat_hm <- seurat_hm_raw
  rm(seurat_hm_raw); gc()
}

# ----------------------------
# ----------------------------
seurat_hm <- NormalizeData(seurat_hm, verbose = FALSE)
seurat_hm <- FindVariableFeatures(seurat_hm, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

top10_hm <- head(VariableFeatures(seurat_hm), 10)
plot_hm1 <- VariableFeaturePlot(seurat_hm) + ggtitle(fig_title("Variable features (before Harmony)"))
plot_hm2 <- LabelPoints(plot = plot_hm1, points = top10_hm, repel = TRUE) +
  ggtitle(fig_title("Top variable features (before Harmony)"))
save_pdf(file.path(OUT_QC_DIR, "Fig13_Variable_features_before_Harmony.pdf"), plot_hm2, w = 7, h = 6)

seurat_hm <- ScaleData(seurat_hm, features = VariableFeatures(seurat_hm), verbose = FALSE)
seurat_hm <- RunPCA(seurat_hm, features = VariableFeatures(seurat_hm), npcs = 50, verbose = FALSE)

# ----------------------------
# ----------------------------
suppressPackageStartupMessages(library(harmony))

pca_emb <- Embeddings(seurat_hm, reduction = "pca")[, 1:30, drop = FALSE]

meta_use <- seurat_hm@meta.data[colnames(seurat_hm), , drop = FALSE]
stopifnot("orig.ident" %in% colnames(meta_use))

set.seed(1)
harmony_emb <- harmony::HarmonyMatrix(
  data_mat = pca_emb,
  meta_data = meta_use,
  vars_use = "orig.ident",
  do_pca = FALSE,
  verbose = TRUE
)

seurat_hm[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_emb,
  key = "harmony_",
  assay = DefaultAssay(seurat_hm)
)

rm(pca_emb, meta_use, harmony_emb); gc()

# ----------------------------
# ----------------------------
seurat_hm <- FindNeighbors(seurat_hm, reduction = "harmony", dims = 1:30, verbose = FALSE)
seurat_hm <- FindClusters(seurat_hm, resolution = 0.5, verbose = FALSE)
seurat_hm <- RunUMAP(seurat_hm, reduction = "harmony", dims = 1:30, verbose = FALSE)

p_harmony <- DimPlot(seurat_hm, reduction = "harmony", pt.size = 0.5) +
  theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5)) +
  ggtitle(fig_title("Harmony embedding"))
save_pdf(file.path(OUT_QC_DIR, "Fig14_Harmony_embedding.pdf"), p_harmony, w = 7.5, h = 5.5)

p_umap_hm <- DimPlot(seurat_hm, reduction = "umap", pt.size = 0.5) +
  theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5)) +
  ggtitle(fig_title("UMAP after Harmony"))
save_pdf(file.path(OUT_QC_DIR, "Fig15_UMAP_after_Harmony.pdf"), p_umap_hm, w = 7.5, h = 5.5)

save_pdf(
  file.path(OUT_QC_DIR, "Fig16_UMAP_split_by_orig.ident.pdf"),
  DimPlot(seurat_hm, reduction = "umap", split.by = "orig.ident", pt.size = 0.5) +
    ggtitle(fig_title("UMAP split by sample (orig.ident)")),
  w = 12, h = 7.5
)

# ----------------------------
# ----------------------------
qsave(seurat_hm, file.path(OUT_QC_DIR, "8.Harmony_corrected.qs"))

message(">>> Done. dim(seurat_hm) = ", paste(dim(seurat_hm), collapse = " x "))
message(">>> RNA layers now: ", paste(Layers(seurat_hm[["RNA"]]), collapse = ", "))

print(dim(seurat_hm))


############################################################-
############################################################-
seurat_clu <- seurat_hm
seurat_clu <- qread(file.path(OUT_QC_DIR, "8.Harmony_corrected.qs"))

save_pdf(file.path(OUT_QC_DIR, "Fig17_ElbowPlot.pdf"),
         ElbowPlot(seurat_clu, ndims = 50) + ggtitle(fig_title("Elbow plot")), w = 5, h = 4)

save_pdf(file.path(OUT_QC_DIR, "Fig18_DimHeatmap_PC1_20.pdf"),
         DimHeatmap(seurat_clu, dims = 1:20, cells = 1000, balanced = TRUE) +
           ggtitle(fig_title("PC1-20 heatmap")), w = 7.5, h = 11)

seuratPC <- 30
seurat_clu <- FindNeighbors(seurat_clu, dims = 1:seuratPC, reduction = "pca")

for (res in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)) {
  seurat_clu <- FindClusters(seurat_clu, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
p_tree <- clustree(seurat_clu@meta.data, prefix = "RNA_snn_res.")
save_pdf(file.path(OUT_QC_DIR, "Fig19_Clustree_resolutions.pdf"), p_tree + ggtitle(fig_title("Clustree across resolutions")),
         w = 25, h = 20)

Idents(seurat_clu) <- seurat_clu@meta.data[["RNA_snn_res.0.5"]]
p_umap_final <- DimPlot(seurat_clu, reduction = "umap", pt.size = 0.1, label = TRUE) +
  ggtitle(fig_title("Final UMAP (resolution = 0.5)")) +
  theme_classic()
save_pdf(file.path(OUT_QC_DIR, "Fig20_UMAP_final_clusters_res0.5.pdf"), p_umap_final, w = 7.5, h = 5.5)

qsave(seurat_clu, file.path(OUT_QC_DIR, "9.Clustered_final.qs"))

############################################################-
############################################################-
#
message(">>> Pipeline finished. Outputs saved to: ", OUT_QC_DIR)


# message("
# Pipeline finished successfully.
# 
# Data: GSE130756 (15 M5+ 4 M4 AML samples  + 3 NORMAL, TXT)
# Steps:
# TXT → QC → Normalize → HVG → PCA/UMAP → Doublet →
# decontX → Cell Cycle → Harmony → Clustering
# 
# All outputs saved in: ", OUT_QC_DIR
# )
