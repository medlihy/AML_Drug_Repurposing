############################################################
## SCI-level CellChat analysis for AML scRNA-seq
## Author: revised for publication-quality figures
############################################################

rm(list = ls())

PATH_CONFIG <- Sys.getenv("AML_PATH_CONFIG", unset = file.path("config", "paths.R"))
if (!file.exists(PATH_CONFIG)) {
  PATH_CONFIG <- file.path("config", "paths_example.R")
}
if (!file.exists(PATH_CONFIG)) {
  stop("Path configuration not found. Run scripts from the repository root or set AML_PATH_CONFIG.")
}
source(PATH_CONFIG, local = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(qs)
  library(openxlsx)
  library(patchwork)
  library(grid)
  library(ComplexHeatmap)
  library(NMF)
  library(ggalluvial)
})

############################################################
## 0. User settings
############################################################

ANNO_DIR <- SCRNA_ANNOTATION_DIR
COMMUNI_DIR <- SCRNA_CELLCHAT_DIR

FIG_DIR <- file.path(COMMUNI_DIR, "figures")
TAB_DIR <- file.path(COMMUNI_DIR, "tables")
OBJ_DIR <- file.path(COMMUNI_DIR, "objects")

dir.create(COMMUNI_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_CELLS <- 20
USE_POPULATION_SIZE <- TRUE

CIRCLE_TOP <- 0.35

K_OUTGOING <- 6
K_INCOMING <- 3

TOP_N_PATHWAYS <- 12

FONT_FAMILY <- "sans"

############################################################
## 1. Plot helper functions
############################################################

open_pdf <- function(file, width = 8, height = 6) {
  grDevices::cairo_pdf(filename = file, width = width, height = height)
}

theme_sci <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(color = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.35, "cm"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(size = base_size, face = "bold", color = "black")
    )
}

save_gg <- function(plot, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(FIG_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

############################################################
## 2. Cell type order and SCI color palette
############################################################

celltype_order <- c(
  "AML_Blast",
  "AML_Progenitor",
  "MonoBlast",
  "Promyelocyte",
  "Myelocyte",
  "IFN-Myelocyte",
  "Neutrophil",
  "Neu&Ret",
  "Eosinophil",
  "ProEosinophil",
  "ProMono",
  "iMono",
  "ncMono",
  "Macrophage",
  "NK",
  "NK_T",
  "ProB",
  "Plasma",
  "EBP",
  "ProEry",
  "Erythroid"
)

celltype_colors <- c(
  "AML_Blast"       = "#5E3C99",
  "AML_Progenitor"  = "#8E63B0",
  "MonoBlast"       = "#B2182B",
  
  "Promyelocyte"    = "#D6604D",
  "Myelocyte"       = "#E6550D",
  "IFN-Myelocyte"   = "#3182BD",
  "Neutrophil"      = "#FDAE6B",
  "Neu&Ret"         = "#BDB76B",
  "Eosinophil"      = "#E6AB02",
  "ProEosinophil"   = "#A6761D",
  
  "ProMono"         = "#1B9E77",
  "iMono"           = "#66C2A5",
  "ncMono"          = "#01665E",
  "Macrophage"      = "#008080",
  
  "NK"              = "#1F78B4",
  "NK_T"            = "#6BAED6",
  "ProB"            = "#33A02C",
  "Plasma"          = "#238B45",
  
  "EBP"             = "#C51B7D",
  "ProEry"          = "#F768A1",
  "Erythroid"       = "#DD1C77"
)

disease_colors <- c(
  "Normal" = "#4D9221",
  "AML"    = "#B2182B",
  "Other"  = "#999999"
)

############################################################
## 3. Load Seurat object and clean metadata
############################################################

# Input is supplied through SCRNA_ANNOTATED_OBJECT_FILE.

sce <- qread(SCRNA_ANNOTATED_OBJECT_FILE)

DefaultAssay(sce) <- "RNA"

stopifnot("celltype" %in% colnames(sce@meta.data))
stopifnot("orig.ident" %in% colnames(sce@meta.data))

sce$celltype <- as.character(sce$celltype)
sce <- subset(sce, subset = !is.na(celltype) & celltype != "")

sce$Sample <- as.character(sce$orig.ident)

sce$Disease <- dplyr::case_when(
  str_detect(sce$Sample, "^AML") ~ "AML",
  str_detect(sce$Sample, "^N")   ~ "Normal",
  TRUE ~ "Other"
)

sce$Disease <- factor(sce$Disease, levels = c("Normal", "AML", "Other"))

real_celltypes <- sort(unique(as.character(sce$celltype)))

cat("Detected cell types:\n")
print(real_celltypes)

missing_from_order <- setdiff(real_celltypes, celltype_order)
missing_from_data  <- setdiff(celltype_order, real_celltypes)

if (length(missing_from_order) > 0) {
  message("以下细胞类型存在于数据中，但不在 celltype_order 中：")
  print(missing_from_order)
  message("它们会被自动追加到最后。")
}

if (length(missing_from_data) > 0) {
  message("以下细胞类型写在 celltype_order 中，但数据中没有：")
  print(missing_from_data)
}

celltype_order_final <- c(
  celltype_order[celltype_order %in% real_celltypes],
  setdiff(real_celltypes, celltype_order)
)

missing_color <- setdiff(celltype_order_final, names(celltype_colors))

if (length(missing_color) > 0) {
  extra_colors <- grDevices::hcl.colors(
    n = length(missing_color),
    palette = "Dark 3"
  )
  names(extra_colors) <- missing_color
  celltype_colors <- c(celltype_colors, extra_colors)
}

celltype_colors <- celltype_colors[celltype_order_final]

sce$celltype <- factor(as.character(sce$celltype), levels = celltype_order_final)
Idents(sce) <- sce$celltype

meta_check <- sce@meta.data %>%
  dplyr::select(orig.ident, Sample, Disease, celltype) %>%
  dplyr::count(Sample, Disease, celltype, name = "CellNumber") %>%
  dplyr::arrange(Disease, Sample, celltype)

openxlsx::write.xlsx(
  meta_check,
  file = file.path(TAB_DIR, "Tab00_celltype_sample_disease_check.xlsx"),
  overwrite = TRUE
)

############################################################
## 4. Normalize and extract data
############################################################

future::plan("sequential")
options(future.globals.maxSize = 20 * 1024^3)

sce <- NormalizeData(
  sce,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

data.input <- tryCatch(
  {
    GetAssayData(sce, assay = "RNA", layer = "data")
  },
  error = function(e) {
    GetAssayData(sce, assay = "RNA", slot = "data")
  }
)

meta <- sce@meta.data
meta <- meta[match(colnames(data.input), rownames(meta)), , drop = FALSE]

stopifnot(all(colnames(data.input) == rownames(meta)))

cat("Genes x Cells:", nrow(data.input), "x", ncol(data.input), "\n")
cat("Cell types:", length(levels(sce$celltype)), "\n")

############################################################
## 5. Create CellChat object
############################################################

cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "celltype"
)

cellchat <- setIdent(cellchat, ident.use = "celltype")

cellchat@idents <- factor(
  as.character(cellchat@idents),
  levels = celltype_order_final
)

celltypes_use <- levels(cellchat@idents)

groupSize <- as.numeric(table(cellchat@idents)[celltypes_use])
names(groupSize) <- celltypes_use

############################################################
## 6. Load CellChat database
############################################################

CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

stopifnot(nrow(cellchat@DB$interaction) > 0)

open_pdf(file.path(FIG_DIR, "Fig01_CellChat_database_category.pdf"), width = 7, height = 5)
showDatabaseCategory(CellChatDB)
dev.off()

############################################################
## 7. Preprocessing for CellChat
############################################################

cellchat <- subsetData(cellchat)

cat("data.signaling dimension:\n")
print(dim(cellchat@data.signaling))

stopifnot(nrow(cellchat@data.signaling) > 0)

future::plan("sequential")
gc()

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

############################################################
## 8. Compute communication probability
############################################################

cellchat <- computeCommunProb(
  cellchat,
  type = "truncatedMean",
  trim = 0.1,
  population.size = USE_POPULATION_SIZE,
  raw.use = TRUE
)

cellchat <- filterCommunication(
  cellchat,
  min.cells = MIN_CELLS
)

df.net <- subsetCommunication(cellchat)

if ("prob" %in% colnames(df.net)) {
  df.net <- df.net %>% arrange(desc(prob))
}

openxlsx::write.xlsx(
  df.net,
  file = file.path(TAB_DIR, "Tab01_CellChat_ligand_receptor_network.xlsx"),
  overwrite = TRUE
)

############################################################
## 9. Pathway-level communication and aggregated network
############################################################

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

qsave(
  cellchat,
  file = file.path(OBJ_DIR, "cellchat_all_after_aggregateNet.qs")
)

############################################################
## 10. Export count and weight matrices
############################################################

mat_count  <- cellchat@net$count
mat_weight <- cellchat@net$weight

mat_count  <- mat_count[celltypes_use, celltypes_use, drop = FALSE]
mat_weight <- mat_weight[celltypes_use, celltypes_use, drop = FALSE]

count_df <- as.data.frame(mat_count) %>%
  tibble::rownames_to_column("Source")

weight_df <- as.data.frame(mat_weight) %>%
  tibble::rownames_to_column("Source")

openxlsx::write.xlsx(
  count_df,
  file = file.path(TAB_DIR, "Tab02_CellChat_interaction_count_matrix.xlsx"),
  overwrite = TRUE
)

openxlsx::write.xlsx(
  weight_df,
  file = file.path(TAB_DIR, "Tab03_CellChat_interaction_weight_matrix.xlsx"),
  overwrite = TRUE
)

############################################################
## 11. Global circle plots
############################################################

open_pdf(file.path(FIG_DIR, "Fig02_Global_network_count_top_edges.pdf"), width = 8.5, height = 8.5)
netVisual_circle(
  mat_count,
  color.use = celltype_colors[celltypes_use],
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions",
  top = CIRCLE_TOP,
  edge.weight.max = max(mat_count, na.rm = TRUE),
  edge.width.max = 6,
  vertex.size.max = 18,
  vertex.label.cex = 0.75,
  alpha.edge = 0.55,
  margin = 0.12
)
dev.off()

open_pdf(file.path(FIG_DIR, "Fig03_Global_network_weight_top_edges.pdf"), width = 8.5, height = 8.5)
netVisual_circle(
  mat_weight,
  color.use = celltype_colors[celltypes_use],
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction strength",
  top = CIRCLE_TOP,
  edge.weight.max = max(mat_weight, na.rm = TRUE),
  edge.width.max = 6,
  vertex.size.max = 18,
  vertex.label.cex = 0.75,
  alpha.edge = 0.55,
  margin = 0.12
)
dev.off()

open_pdf(file.path(FIG_DIR, "FigS01_Global_network_weight_all_edges.pdf"), width = 9, height = 9)
netVisual_circle(
  mat_weight,
  color.use = celltype_colors[celltypes_use],
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction strength - all edges",
  top = 1,
  edge.weight.max = max(mat_weight, na.rm = TRUE),
  edge.width.max = 5,
  vertex.size.max = 17,
  vertex.label.cex = 0.65,
  alpha.edge = 0.35,
  margin = 0.12
)
dev.off()

############################################################
## 12. Communication heatmaps
############################################################

open_pdf(file.path(FIG_DIR, "Fig04_Communication_heatmap_weight.pdf"), width = 7.5, height = 7)
netVisual_heatmap(
  cellchat,
  measure = "weight",
  color.heatmap = "Reds"
)
dev.off()

open_pdf(file.path(FIG_DIR, "Fig05_Communication_heatmap_count.pdf"), width = 7.5, height = 7)
netVisual_heatmap(
  cellchat,
  measure = "count",
  color.heatmap = "Blues"
)
dev.off()

############################################################
## 13. Focused malignant-myeloid-immune network
############################################################

focus_cells <- c(
  "AML_Blast",
  "AML_Progenitor",
  "MonoBlast",
  "ProMono",
  "iMono",
  "ncMono",
  "Macrophage",
  "Neutrophil",
  "NK",
  "NK_T",
  "ProB",
  "Plasma"
)

focus_cells <- focus_cells[focus_cells %in% celltypes_use]

mat_focus_weight <- mat_weight[focus_cells, focus_cells, drop = FALSE]
mat_focus_count  <- mat_count[focus_cells, focus_cells, drop = FALSE]

groupSize_focus <- groupSize[focus_cells]

open_pdf(file.path(FIG_DIR, "Fig06_Focused_AML_myeloid_immune_network_weight.pdf"), width = 7.5, height = 7.5)
netVisual_circle(
  mat_focus_weight,
  color.use = celltype_colors[focus_cells],
  vertex.weight = groupSize_focus,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Focused AML-myeloid-immune network",
  top = 0.55,
  edge.weight.max = max(mat_focus_weight, na.rm = TRUE),
  edge.width.max = 7,
  vertex.size.max = 20,
  vertex.label.cex = 0.85,
  alpha.edge = 0.60,
  margin = 0.12
)
dev.off()

############################################################
## 14. Single source outgoing networks
############################################################

open_pdf(file.path(FIG_DIR, "FigS02_Single_source_outgoing_networks.pdf"), width = 12, height = 12)

par(mfrow = c(3, 3), xpd = TRUE, mar = c(1, 1, 2, 1))

for (cel in celltypes_use) {
  cir_mat <- matrix(
    0,
    nrow = nrow(mat_weight),
    ncol = ncol(mat_weight),
    dimnames = dimnames(mat_weight)
  )
  
  cir_mat[cel, ] <- mat_weight[cel, ]
  
  netVisual_circle(
    cir_mat,
    color.use = celltype_colors[celltypes_use],
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_weight, na.rm = TRUE),
    edge.width.max = 6,
    vertex.size.max = 16,
    vertex.label.cex = 0.65,
    label.edge = FALSE,
    title.name = paste0(cel, " outgoing"),
    alpha.edge = 0.60,
    margin = 0.08
  )
}

dev.off()

############################################################
## 15. Network centrality and signaling role heatmaps
############################################################

cellchat <- netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)

qsave(
  cellchat,
  file = file.path(OBJ_DIR, "cellchat_all_after_centrality.qs")
)

open_pdf(file.path(FIG_DIR, "Fig07_Signaling_role_heatmap_outgoing_incoming.pdf"), width = 10, height = 7)

ht_out <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing",
  color.heatmap = "Reds",
  width = unit(4.2, "cm"),
  height = unit(12, "cm")
)

ht_in <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming",
  color.heatmap = "Blues",
  width = unit(4.2, "cm"),
  height = unit(12, "cm")
)

ComplexHeatmap::draw(
  ht_out + ht_in,
  ht_gap = unit(0.5, "cm")
)

dev.off()

############################################################
## 16. Signaling role scatter
############################################################

gg_role <- netAnalysis_signalingRole_scatter(
  cellchat,
  color.use = celltype_colors[celltypes_use]
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey30"
  ) +
  labs(
    x = "Outgoing interaction strength",
    y = "Incoming interaction strength"
  ) +
  theme_sci(base_size = 10) +
  theme(
    legend.position = "right"
  )

save_gg(
  gg_role,
  filename = "Fig08_Signaling_role_scatter.pdf",
  width = 7,
  height = 5.5
)

############################################################
## 17. Top pathway ranking
############################################################

pathway_strength <- apply(cellchat@netP$prob, 3, sum)
pathway_strength <- sort(pathway_strength, decreasing = TRUE)

pathway_rank <- data.frame(
  Pathway = names(pathway_strength),
  TotalStrength = as.numeric(pathway_strength)
)

openxlsx::write.xlsx(
  pathway_rank,
  file = file.path(TAB_DIR, "Tab04_CellChat_pathway_strength_rank.xlsx"),
  overwrite = TRUE
)

top_pathways <- names(pathway_strength)[1:min(TOP_N_PATHWAYS, length(pathway_strength))]

cat("Top pathways:\n")
print(top_pathways)

############################################################
## 18. Top pathway bubble plot
############################################################

gg_bubble <- netVisual_bubble(
  cellchat,
  signaling = top_pathways,
  remove.isolate = FALSE,
  color.heatmap = "viridis",
  angle.x = 45,
  font.size = 9
) +
  theme_sci(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )

save_gg(
  gg_bubble,
  filename = "Fig09_Top_pathway_ligand_receptor_bubbleplot.pdf",
  width = 12,
  height = 7
)

############################################################
## 19. Pathway-specific network plots
############################################################

pathways_to_plot <- top_pathways[1:min(6, length(top_pathways))]

for (pw in pathways_to_plot) {
  out_file <- file.path(
    FIG_DIR,
    paste0("FigS03_Pathway_network_", make.names(pw), ".pdf")
  )
  
  open_pdf(out_file, width = 7.2, height = 7.2)
  
  netVisual_aggregate(
    cellchat,
    signaling = pw,
    layout = "circle",
    color.use = celltype_colors[celltypes_use],
    vertex.weight = groupSize,
    vertex.size.max = 18,
    vertex.label.cex = 0.75,
    edge.width.max = 6,
    alpha.edge = 0.55,
    margin = 0.12
  )
  
  dev.off()
}

############################################################
## 20. Select K for outgoing communication patterns
############################################################

open_pdf(file.path(FIG_DIR, "Fig10_SelectK_outgoing.pdf"), width = 8, height = 6)
selectK(cellchat, pattern = "outgoing")
dev.off()

open_pdf(file.path(FIG_DIR, paste0("Fig11_Outgoing_communication_patterns_k", K_OUTGOING, ".pdf")),
         width = 10, height = 5.5)

cellchat <- identifyCommunicationPatterns(
  cellchat,
  pattern = "outgoing",
  k = K_OUTGOING
)

dev.off()

open_pdf(file.path(FIG_DIR, paste0("Fig12_Outgoing_river_plot_k", K_OUTGOING, ".pdf")),
         width = 11, height = 7)

netAnalysis_river(
  cellchat,
  pattern = "outgoing"
)

dev.off()

open_pdf(file.path(FIG_DIR, paste0("Fig13_Outgoing_dot_plot_k", K_OUTGOING, ".pdf")),
         width = 10, height = 6.5)

netAnalysis_dot(
  cellchat,
  pattern = "outgoing"
)

dev.off()

############################################################
## 21. Select K for incoming communication patterns
############################################################

open_pdf(file.path(FIG_DIR, "Fig14_SelectK_incoming.pdf"), width = 8, height = 6)
selectK(cellchat, pattern = "incoming")
dev.off()

open_pdf(file.path(FIG_DIR, paste0("Fig15_Incoming_communication_patterns_k", K_INCOMING, ".pdf")),
         width = 10, height = 5.5)

cellchat <- identifyCommunicationPatterns(
  cellchat,
  pattern = "incoming",
  k = K_INCOMING
)

dev.off()

open_pdf(file.path(FIG_DIR, paste0("Fig16_Incoming_river_plot_k", K_INCOMING, ".pdf")),
         width = 11, height = 7)

netAnalysis_river(
  cellchat,
  pattern = "incoming"
)

dev.off()

open_pdf(file.path(FIG_DIR, paste0("Fig17_Incoming_dot_plot_k", K_INCOMING, ".pdf")),
         width = 10, height = 6.5)

netAnalysis_dot(
  cellchat,
  pattern = "incoming"
)

dev.off()

############################################################
## 22. Save final CellChat object
############################################################

qsave(
  cellchat,
  file = file.path(OBJ_DIR, "cellchat_all_final.qs")
)

cat("\nCellChat SCI-level analysis finished.\n")
cat("Figures saved to:", FIG_DIR, "\n")
cat("Tables saved to:", TAB_DIR, "\n")
cat("Objects saved to:", OBJ_DIR, "\n")
