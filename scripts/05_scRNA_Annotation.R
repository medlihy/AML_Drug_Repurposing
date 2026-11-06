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
library(ggsci)
library(dplyr)
library(openxlsx)
library(patchwork)
library(readr)
library(stringr)

QC_DIR <- SCRNA_QC_DIR
ANNO_DIR <- SCRNA_ANNOTATION_DIR
if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)

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
seurat_obj <- qread(file.path(QC_DIR, "9.Clustered_final.qs"))

head(seurat_obj@meta.data)

px=0.5
seurat_obj <- FindClusters(seurat_obj, resolution = px)

DefaultAssay(seurat_obj)
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

seurat_obj.markers$difference <- seurat_obj.markers$pct.1 - seurat_obj.markers$pct.2

seurat_obj.markers_sorted <- seurat_obj.markers %>%
  group_by(cluster) %>%
  arrange(desc(difference)) %>%
  ungroup()


file_name <- sprintf("Tab01_resolution_%.1f_marker_gene.xlsx", px)
write.xlsx(seurat_obj.markers_sorted, file = file_name)


top_genes_per_cluster <- lapply(unique(seurat_obj.markers_sorted$cluster), function(cluster) {
  top_genes <- seurat_obj.markers_sorted[seurat_obj.markers_sorted$cluster == cluster,]
  top_genes <- head(top_genes$gene, 10)
  paste(top_genes, collapse = ",")
})
names(top_genes_per_cluster) <- unique(seurat_obj.markers_sorted$cluster)

print(head(seurat_obj.markers_sorted, 3))

txt_content <- sapply(names(top_genes_per_cluster), function(cluster) {
  paste0("cluster_", cluster, ": ", top_genes_per_cluster[[cluster]])
})
txt_content <- as.character(txt_content)
txt_filename <- sprintf("Tab02_resolution_%.1f_top10_marker_gene.txt", px)
writeLines(txt_content, con = txt_filename)


top5seurat_obj.markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = difference)

col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))

pdf_filename <- sprintf("Fig01_resolution_%.1f_heatmap_plot.pdf", px)

pdf(file = pdf_filename, width = 22, height = 16)

DoHeatmap(
  seurat_obj,
  features = top5seurat_obj.markers$gene,
  group.colors = col
) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(
    low = '#0099CC', mid = 'white', high = '#CC0033',
    name = 'Z-score'
  )

dev.off()


seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, reduction = "pca")
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20, reduction = "pca")


combined_pdf <- sprintf("Fig02_resolution_%.1f_umap_tsne_plot_1.pdf", px)
pdf(file = combined_pdf, width = 12, height = 5.5)
p_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  label.size = 3.5,
  pt.size = 0.1
) +
  ggtitle("UMAP") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.position = "right"
  )

p_tsne <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  label = TRUE,
  label.size = 3.5,
  pt.size = 0.1
) +
  ggtitle("tSNE") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.position = "right"
  )
print(p_umap + p_tsne)
dev.off()


head(seurat_obj@meta.data)
res_col <- if (px == 0.5) {
  "RNA_snn_res.0.5"
} else {
  sprintf("RNA_snn_res.%.1f", px)
}


seurat_obj@meta.data <- seurat_obj@meta.data[, c("orig.ident", "type", res_col)]

rds_file <- sprintf("1.resolution_%.1f_seuratUMAP.qs", px)

qsave(seurat_obj, rds_file)


seurat_obj <- qread("1.resolution_0.5_seuratUMAP.qs")

file_path <- SCRNA_ANNOTATION_REFERENCE_FILE

df <- read.delim(
  file(file_path, encoding = "GB18030"),
  stringsAsFactors = FALSE,
  header = TRUE
)
df$ClusterNum <- as.numeric(gsub("cluster_", "", df$Cluster))

df$Cell.Type.Clean <- str_trim(gsub("[^a-zA-Z0-9 ,]", "", df$Cell.Type))

df$GeneList <- strsplit(df$whg.overlap, ",")

df$GeneListTop3 <- lapply(df$GeneList, function(x) {
  if (is.null(x)) return(NULL)
  x <- lapply(x, function(gene) gsub(" ", "", gene))
  x <- unlist(x)
  head(x, 3)
})

df <- df[order(df$ClusterNum), ]
marker_lines <- c("marker_cluster <- c(")

na_messages <- c()

for (i in seq_len(nrow(df))) {
  cluster_id <- df$ClusterNum[i]
  
  genes <- df$GeneListTop3[[i]]
  if (is.null(genes) || all(is.na(genes))) {
    genes <- c("NA", "NA", "NA")
    na_messages <- c(na_messages, sprintf("第 %s 簇的 whg.overlap 是 NA", cluster_id))
  } else {
    genes <- trimws(genes)
    genes[genes == ""] <- "NA"
    genes[is.na(genes)] <- "NA"
  }
  
  celltype_clean <- df$Cell.Type.Clean[i]
  if (is.na(celltype_clean) || celltype_clean == "") {
    celltype_clean <- "NA"
    na_messages <- c(na_messages, sprintf("第 %s 簇的 Cell.Type 是 NA", cluster_id))
  }
  
  gene_line <- paste0('  "', paste(genes, collapse = '", "'), '"')
  if (i < nrow(df)) {
    gene_line <- paste0(gene_line, ",")
  }
  comment <- sprintf(" #%s. %s", cluster_id, celltype_clean)
  marker_lines <- c(marker_lines, paste0(gene_line, comment))
}

marker_lines <- c(marker_lines, ")")

writeLines(marker_lines, "1_marker_cluster_output.R")
cat(paste(marker_lines, collapse = "\n"))

if (length(na_messages) > 0) {
  cat("\n# Missing annotation information:\n")
  for (msg in na_messages) cat(msg, "\n")
}


# FCER1A/CLEC10A/CLEC9A/IL3RA/GZMB/LST1/TYROBP/FCN1
# p6
# ggsave("gene_feature.pdf", plot = p6, device = "pdf", width = 15, height = 10)


source("1_marker_cluster_output.R")
marker_cluster <- unique(marker_cluster)

dotplot <- DotPlot(seurat_obj, features = marker_cluster, assay = 'RNA') + coord_flip()
dotplot_file <- sprintf("Fig03_resolution_%.1f_dotplot.pdf", px)
ggsave(dotplot_file, plot = dotplot, width = 15, height = 12)


heatmap <- DoHeatmap(seurat_obj, features = marker_cluster) + NoLegend()
heatmap_file <- sprintf("Fig04_resolution_%.1f_marker_heatmap.pdf", px)
ggsave(heatmap_file, plot = heatmap, width = 15, height = 12)


df$ClusterNum <- as.numeric(gsub("cluster_", "", df$Cluster))

df$Cell.Type.Clean <- str_trim(gsub("[^a-zA-Z0-9 ,]", "", df$Cell.Type))

df$Cell.Type.Clean[df$Cell.Type.Clean == ""] <- "NA"
df$Cell.Type.Clean[is.na(df$Cell.Type.Clean)] <- "NA"

df <- df[order(df$ClusterNum), ]

cluster_lines <- c("new.cluster.ids <- c(")

for (i in seq_len(nrow(df))) {
  cell_type <- df$Cell.Type.Clean[i]
  cluster_id <- df$ClusterNum[i]
  
  line <- paste0('  "', cell_type, '"')
  
  if (i < nrow(df)) {
    line <- paste0(line, ",")
  }
  
  comment <- sprintf(" #%d", cluster_id)
  cluster_lines <- c(cluster_lines, paste0(line, comment))
}
cluster_lines <- c(cluster_lines, ")")
writeLines(cluster_lines, "1_new_cluster_ids_output.R")
cat(paste(cluster_lines, collapse = "\n"))


source("1_new_cluster_ids_output.R")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

cluster_col <- paste0("RNA_snn_res.", px)
if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
  stop(paste("找不到聚类列:", cluster_col))
}
cluster_ids <- as.integer(as.character(seurat_obj@meta.data[[cluster_col]]))
seurat_obj$celltype <- ifelse(
  !is.na(cluster_ids) & cluster_ids < length(new.cluster.ids),
  new.cluster.ids[cluster_ids + 1],
  "Unknown"
)
head(seurat_obj@meta.data)


library(ggplot2)
library(ggrepel)
library(ggunchull)
library(dplyr)
library(stringr)
library(viridis)
meta <- cbind(seurat_obj@meta.data, as.data.frame(seurat_obj@reductions$tsne@cell.embeddings))
colnames(meta)[(ncol(meta)-1):ncol(meta)] <- c("tSNE_1", "tSNE_2")

meta$cell_type <- as.character(meta$celltype)
meta$cell_type[is.na(meta$celltype)] <- "NA"
meta$main_type_wrapped <- str_wrap(meta$cell_type, width = 20)

main_type_med <- meta %>%
  group_by(main_type_wrapped) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2),
    .groups = "drop"
  )

cell_types <- unique(meta$cell_type)
non_na_types <- setdiff(cell_types, "NA")
cols <- viridis(length(non_na_types), option = "D")
names(cols) <- non_na_types
if ("NA" %in% cell_types) cols <- c(cols, "NA" = "grey80")

x_limits <- c(min(meta$tSNE_1) - 1, max(meta$tSNE_1) + 1)
y_limits <- c(min(meta$tSNE_2) - 1, max(meta$tSNE_2) + 1)

annotated_tSNE <- ggplot(meta, aes(x = tSNE_1, y = tSNE_2)) +
  stat_unchull(aes(fill = cell_type, color = cell_type),
               alpha = 0.15, linewidth = 0.8, linetype = "dashed", delta = 0.3,
               show.legend = FALSE) +
  geom_point(aes(color = cell_type), size = 0.4, show.legend = FALSE) +
  geom_text_repel(data = main_type_med,
                  aes(label = main_type_wrapped),
                  fontface = "bold",
                  size = 5,
                  box.padding = 0.5,
                  segment.size = 0.3,
                  max.overlaps = 100,
                  color = "black",
                  show.legend = FALSE) +
  scale_x_continuous(name = "tSNE 1", breaks = NULL, limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(name = "tSNE 2", breaks = NULL, limits = y_limits, expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(linewidth = 0.6,
                             arrow = arrow(type = "closed", length = unit(0.25, "cm")))
  ) +
  ggtitle("")

umap_file <- sprintf("Fig05_resolution_%.1f_celltype_tSNE.pdf", px)
ggsave(umap_file, plot = annotated_tSNE, width = 10, height = 8)


meta <- cbind(seurat_obj@meta.data, as.data.frame(seurat_obj@reductions$umap@cell.embeddings))
colnames(meta)[(ncol(meta)-1):ncol(meta)] <- c("UMAP_1", "UMAP_2")

meta$cell_type <- as.character(meta$celltype)
meta$cell_type[is.na(meta$celltype)] <- "NA"
meta$main_type_wrapped <- str_wrap(meta$cell_type, width = 20)

main_type_med <- meta %>%
  group_by(main_type_wrapped) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

cell_types <- unique(meta$cell_type)
non_na_types <- setdiff(cell_types, "NA")
cols <- viridis(length(non_na_types), option = "D")
names(cols) <- non_na_types
if ("NA" %in% cell_types) cols <- c(cols, "NA" = "grey80")

x_limits <- c(min(meta$UMAP_1) - 1, max(meta$UMAP_1) + 1)
y_limits <- c(min(meta$UMAP_2) - 1, max(meta$UMAP_2) + 1)

annotated_umap <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_unchull(aes(fill = cell_type, color = cell_type),
               alpha = 0.15, linewidth = 0.8, linetype = "dashed", delta = 0.3,
               show.legend = FALSE) +
  geom_point(aes(color = cell_type), size = 0.4, show.legend = FALSE) +
  geom_text_repel(data = main_type_med,
                  aes(label = main_type_wrapped),
                  fontface = "bold",
                  size = 5,
                  box.padding = 0.5,
                  segment.size = 0.3,
                  max.overlaps = 100,
                  color = "black",
                  show.legend = FALSE) +
  scale_x_continuous(name = "UMAP 1", breaks = NULL, limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(name = "UMAP 2", breaks = NULL, limits = y_limits, expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(linewidth = 0.6,
                             arrow = arrow(type = "closed", length = unit(0.25, "cm")))
  ) +
  ggtitle("")

umap_file <- sprintf("Fig05_resolution_%.1f_celltype_umap.pdf", px)
ggsave(umap_file, plot = annotated_umap, width = 10, height = 8)


qsave(seurat_obj,"2.resolution_0.5_annotation.qs")


# 
# library(Seurat)
# library(wordcloud)
# library(RColorBrewer)
# 
# 
# seurat_object<-qread("2.resolution_0.1_annotation.qs")
# seurat_object$celltype<-seurat_object@active.ident
# 
# colnames(seurat_object@meta.data)
# circ_data <- prepare_circlize_data(seurat_object, scale = 0.8 )
# library(colorspace)
# cluster_colors <- rainbow_hcl(length(levels(seurat_object)))
# group_colors <- rainbow_hcl(length(names(table(seurat_object$celltype))))
# rep_colors <- rainbow_hcl(length(names(table(seurat_object$orig.ident))))
# # png(filename =  '01.circlize_plot.png', width = 8, height = 8,units = 'in', res = 300)
# # plot_circlize(circ_data, do.label = T, pt.size = 0.01, 
# #               col.use = cluster_colors, bg.color = 'white',
# #               kde2d.n = 500, repel = TRUE, label.cex = 1,
# # add_track(circ_data, group = "celltype", colors = group_colors, track_num = 4)
# # add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num =5)
# # dev.off()
# library(circlize)
# pdf('Fig06_circlize_plot.pdf', width = 8, height = 8)
# plot_circlize(circ_data, do.label = T, pt.size = 0.01, 
#               col.use = cluster_colors, bg.color = 'white',
#               kde2d.n = 500, repel = TRUE, label.cex = 1,
# add_track(circ_data, group = "celltype", colors = group_colors, track_num = 4)
# add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num =5)
# dev.off()


# seurat_obj$celltype <- gsub(
#   "Naive/Memory T cells",
#   "Naive or Memory T cells",
#   seurat_obj$celltype
# )
# Idents(seurat_obj) <- seurat_obj$celltype
# levels(seurat_obj)

head(seurat_obj@meta.data)

cell_types <- levels(seurat_obj)

for (cell_type in cell_types) {
  cell_type_seurat <- subset(seurat_obj, subset = celltype == cell_type)
  
  file_name <- paste0("3.", gsub(" ", "_", cell_type), "_seurat.qs")
  
  qsave(cell_type_seurat, file_name)
}
