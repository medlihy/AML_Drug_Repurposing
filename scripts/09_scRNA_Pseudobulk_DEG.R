rm(list = ls())

PATH_CONFIG <- Sys.getenv("AML_PATH_CONFIG", unset = file.path("config", "paths.R"))
if (!file.exists(PATH_CONFIG)) {
  PATH_CONFIG <- file.path("config", "paths_example.R")
}
if (!file.exists(PATH_CONFIG)) {
  stop("Path configuration not found. Run scripts from the repository root or set AML_PATH_CONFIG.")
}
source(PATH_CONFIG, local = FALSE)

library(Seurat)
library(qs)
library(Matrix)
library(dplyr)
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

# =========================
# =========================

ANNO_DIR <- SCRNA_ANNOTATION_DIR
OUT_DIR <- SCRNA_PSEUDOBULK_DIR

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

setwd(ANNO_DIR)

seurat_obj <- qread(SCRNA_ANNOTATED_OBJECT_FILE)

setwd(OUT_DIR)

DefaultAssay(seurat_obj) <- "RNA"

# =========================
# =========================

seurat_obj$sample_id <- seurat_obj$orig.ident

seurat_obj$group <- ifelse(
  grepl("^AML", seurat_obj$sample_id),
  "AML",
  "Normal"
)

print(table(seurat_obj$group))
print(table(seurat_obj$sample_id))
print(sort(unique(seurat_obj$celltype)))

celltype_group_tab <- as.data.frame(table(seurat_obj$group, seurat_obj$celltype))
colnames(celltype_group_tab) <- c("group", "celltype", "n")

fwrite(
  celltype_group_tab,
  "Tab00_celltype_count_by_group.tsv",
  sep = "\t"
)

sample_celltype_tab <- as.data.frame(table(seurat_obj$sample_id, seurat_obj$celltype))
colnames(sample_celltype_tab) <- c("sample_id", "celltype", "n")

fwrite(
  sample_celltype_tab,
  "Tab00_celltype_count_by_sample.tsv",
  sep = "\t"
)

# =========================
# =========================

aml_malignant_celltypes <- c(
  "AML_Blast",
  "AML_Progenitor",
  "MonoBlast"
)

missing_ct <- setdiff(aml_malignant_celltypes, unique(seurat_obj$celltype))
if (length(missing_ct) > 0) {
  stop("这些细胞类型在对象中不存在：", paste(missing_ct, collapse = ", "))
}

meta <- seurat_obj@meta.data

# AML malignant compartment
aml_cells <- rownames(meta)[
  meta$group == "AML" &
    meta$celltype %in% aml_malignant_celltypes
]

# Normal bone marrow total cells
normal_cells <- rownames(meta)[
  meta$group == "Normal"
]

length(aml_cells)
length(normal_cells)

obj_deg <- subset(
  seurat_obj,
  cells = c(aml_cells, normal_cells)
)

obj_deg$compare_group <- ifelse(
  obj_deg$group == "AML",
  "AML_malignant",
  "Normal_BM"
)

print(table(obj_deg$compare_group))
print(table(obj_deg$sample_id, obj_deg$compare_group))


# =========================
# =========================

cell_count_tab <- obj_deg@meta.data %>%
  dplyr::count(sample_id, group, compare_group, name = "n") %>%
  dplyr::arrange(compare_group, sample_id)

print(cell_count_tab)

fwrite(
  cell_count_tab,
  "Tab01_cells_used_for_pseudobulk_by_sample.tsv",
  sep = "\t"
)

MIN_CELLS_PER_SAMPLE <- 50

valid_samples <- cell_count_tab %>%
  dplyr::filter(n >= MIN_CELLS_PER_SAMPLE) %>%
  dplyr::pull(sample_id)

obj_deg <- subset(
  obj_deg,
  cells = rownames(obj_deg@meta.data)[obj_deg$sample_id %in% valid_samples]
)

print(table(obj_deg$sample_id, obj_deg$compare_group))
print(table(obj_deg$compare_group))

sample_group_check <- obj_deg@meta.data %>%
  dplyr::distinct(sample_id, compare_group) %>%
  dplyr::count(compare_group, name = "sample_n")

print(sample_group_check)

fwrite(
  sample_group_check,
  "Tab01_sample_number_after_filter.tsv",
  sep = "\t"
)

if (any(sample_group_check$sample_n < 3)) {
  warning("某一组样本数 < 3，pseudo-bulk DEG 结果需要谨慎解释。")
}
# =========================
# =========================

DefaultAssay(obj_deg) <- "RNA"
Layers(obj_deg[["RNA"]])
counts <- GetAssayData(obj_deg, assay = "RNA", layer = "counts")

meta_deg <- obj_deg@meta.data

samples <- unique(meta_deg$sample_id)

pb_counts <- do.call(cbind, lapply(samples, function(sid) {
  cells_use <- rownames(meta_deg)[meta_deg$sample_id == sid]
  Matrix::rowSums(counts[, cells_use, drop = FALSE])
}))

colnames(pb_counts) <- samples
pb_counts <- as.matrix(pb_counts)

pb_meta <- meta_deg %>%
  dplyr::distinct(sample_id, group, compare_group) %>%
  as.data.frame()

rownames(pb_meta) <- pb_meta$sample_id
pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]

stopifnot(identical(rownames(pb_meta), colnames(pb_counts)))

fwrite(
  data.frame(gene = rownames(pb_counts), pb_counts, check.names = FALSE),
  "Tab02_pseudobulk_raw_counts.tsv",
  sep = "\t"
)

fwrite(
  data.frame(sample_id = rownames(pb_meta), pb_meta, check.names = FALSE),
  "Tab03_pseudobulk_metadata.tsv",
  sep = "\t"
)
# =========================
# 6. limma-voom DEG
# =========================

pb_meta$compare_group <- factor(
  pb_meta$compare_group,
  levels = c("Normal_BM", "AML_malignant")
)

print(table(pb_meta$compare_group))

dge <- DGEList(
  counts = pb_counts,
  group = pb_meta$compare_group
)

keep <- filterByExpr(
  dge,
  group = pb_meta$compare_group
)

dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + compare_group, data = pb_meta)
colnames(design) <- levels(pb_meta$compare_group)

print(design)

# voom
pdf("Fig00_voom_mean_variance_trend.pdf", width = 6, height = 5)
v <- voom(dge, design, plot = TRUE)
dev.off()

fit <- lmFit(v, design)

contrast.matrix <- makeContrasts(
  AML_malignant_vs_Normal_BM = AML_malignant - Normal_BM,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, robust = TRUE)

deg <- topTable(
  fit2,
  coef = "AML_malignant_vs_Normal_BM",
  number = Inf,
  adjust.method = "BH"
)

deg$gene <- rownames(deg)

deg <- deg %>%
  dplyr::select(gene, everything())

# =========================
# =========================

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

deg <- deg %>%
  dplyr::mutate(
    regulated = dplyr::case_when(
      adj.P.Val < FDR_CUTOFF & logFC > LOGFC_CUTOFF  ~ "Up_in_AML_malignant",
      adj.P.Val < FDR_CUTOFF & logFC < -LOGFC_CUTOFF ~ "Down_in_AML_malignant",
      TRUE ~ "Not_sig"
    )
  )

print(table(deg$regulated))

deg_sig <- deg %>%
  dplyr::filter(regulated != "Not_sig")

deg_up <- deg %>%
  dplyr::filter(regulated == "Up_in_AML_malignant")

deg_down <- deg %>%
  dplyr::filter(regulated == "Down_in_AML_malignant")

fwrite(
  deg,
  "Tab04_AML_malignant_vs_NormalBM_pseudobulk_DEG_all.tsv",
  sep = "\t"
)

fwrite(
  deg_sig,
  "Tab05_AML_malignant_vs_NormalBM_pseudobulk_DEG_sig.tsv",
  sep = "\t"
)

fwrite(
  deg_up,
  "Tab06_AML_malignant_vs_NormalBM_pseudobulk_DEG_up.tsv",
  sep = "\t"
)

fwrite(
  deg_down,
  "Tab07_AML_malignant_vs_NormalBM_pseudobulk_DEG_down.tsv",
  sep = "\t"
)

# =========================
# =========================

CMAP_N <- 150

cmap_up <- deg %>%
  dplyr::filter(regulated == "Up_in_AML_malignant") %>%
  dplyr::arrange(adj.P.Val, dplyr::desc(logFC)) %>%
  dplyr::slice_head(n = CMAP_N) %>%
  dplyr::pull(gene) %>%
  unique()

cmap_down <- deg %>%
  dplyr::filter(regulated == "Down_in_AML_malignant") %>%
  dplyr::arrange(adj.P.Val, logFC) %>%
  dplyr::slice_head(n = CMAP_N) %>%
  dplyr::pull(gene) %>%
  unique()

writeLines(
  cmap_up,
  "CMap_AML_malignant_vs_NormalBM_UP_top150.txt"
)

writeLines(
  cmap_down,
  "CMap_AML_malignant_vs_NormalBM_DOWN_top150.txt"
)

length(cmap_up)
length(cmap_down)

# =========================
# =========================

deg$plot_p <- deg$P.Value
deg$plot_p[deg$plot_p == 0] <- 1e-300

deg$regulated <- factor(
  deg$regulated,
  levels = c(
    "Up_in_AML_malignant",
    "Not_sig",
    "Down_in_AML_malignant"
  )
)

top_label <- deg %>%
  dplyr::filter(regulated != "Not_sig") %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = 20)

p_volcano <- ggplot(deg, aes(x = logFC, y = -log10(plot_p))) +
  geom_point(
    aes(color = regulated),
    alpha = 0.75,
    size = 1.1
  ) +
  scale_color_manual(
    values = c(
      "Up_in_AML_malignant" = "#D73027",
      "Not_sig" = "grey75",
      "Down_in_AML_malignant" = "#4575B4"
    ),
    breaks = c(
      "Up_in_AML_malignant",
      "Not_sig",
      "Down_in_AML_malignant"
    ),
    labels = c(
      "Up in AML",
      "Not sig.",
      "Down in AML"
    )
  ) +
  geom_vline(
    xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  ggrepel::geom_text_repel(
    data = top_label,
    aes(label = gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.25
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    legend.position = "right",
    legend.box.spacing = grid::unit(0.02, "cm"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.45, "cm"),
    panel.grid = element_blank()
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(size = 3, alpha = 1)
    )
  ) +
  labs(
    title = "AML malignant vs Normal BM",
    x = "log2FC",
    y = "-log10(P value)",
    color = NULL
  )

ggsave(
  "Fig01_AML_malignant_vs_NormalBM_pseudobulk_volcano.pdf",
  p_volcano,
  width = 5.8,
  height = 5.8
)

# =========================
# 10. pseudo-bulk PCA
# =========================

logCPM <- cpm(dge, log = TRUE, prior.count = 1)

pca <- prcomp(t(logCPM), scale. = TRUE)

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  pb_meta[rownames(pca$x), , drop = FALSE]
)

pc_var <- summary(pca)$importance[2, 1:2] * 100

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = compare_group, label = sample_id)
) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  theme_bw(base_size = 13) +
  labs(
    title = "Pseudo-bulk PCA",
    x = paste0("PC1: ", round(pc_var[1], 1), "%"),
    y = paste0("PC2: ", round(pc_var[2], 1), "%"),
    color = NULL
  )

ggsave(
  "Fig02_pseudobulk_PCA.pdf",
  p_pca,
  width = 7,
  height = 6
)

# =========================
# =========================

marker_genes <- c(
  "MPO", "ELANE", "AZU1", "PRTN3",
  "LYZ", "S100A8", "S100A9",
  "CD34", "KIT", "CEBPA",
  "MKI67", "TOP2A", "PCNA", "TYMS",
  "H2AFZ", "STMN1", "MIF"
)

marker_check <- deg %>%
  dplyr::filter(gene %in% marker_genes) %>%
  dplyr::arrange(dplyr::desc(logFC))

fwrite(
  marker_check,
  "Tab08_marker_gene_DEG_check.tsv",
  sep = "\t"
)

marker_check
