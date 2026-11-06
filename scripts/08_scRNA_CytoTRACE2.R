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
library(CytoTRACE2)
library(qs)
library(ggplot2)
ANNO_DIR <- SCRNA_ANNOTATION_DIR
TRACE_DIR <- SCRNA_CYTOTRACE_DIR


if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)
if (!dir.exists(TRACE_DIR)) dir.create(TRACE_DIR, recursive = TRUE)
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
seurat_obj <- qread(SCRNA_ANNOTATED_OBJECT_FILE)


setwd(TRACE_DIR)
cytotrace2_result <- cytotrace2(seurat_obj, is_seurat = TRUE,
                                slot_type = "counts",
                                species = "human",
                                ncores = 5)


phe.anno = as.data.frame(cytotrace2_result$celltype)
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  expression_data = NULL, 
                  is_seurat = TRUE,
                  annotation = phe.anno, 
                  slot_type = "counts",
                  pc_dims = 50)


p2 <- plots$CytoTRACE2_Boxplot_byPheno + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 15),
        ...)
pdf("CytoTRACE2_Boxplot.pdf", width = 8, height = 5)
p2
dev.off()

library(ggplot2)
library(ggsignif)
library(patchwork)
library(dplyr)
pub_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1),
    legend.position = "none"
  )

p_umap <- plots$CytoTRACE2_Potency_UMAP +
  ggtitle("CytoTRACE2 Differentiation Potency") +
  pub_theme +
  theme(legend.position = "right")

# p_box <- plots$CytoTRACE2_Boxplot_byPheno +
#   ggtitle("CytoTRACE2 Score by Cell Type") +
#   pub_theme +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "none"
#   ) +
#   ylab("CytoTRACE2 Potency Score") + xlab("Cell Type")

  p_box <- plots$CytoTRACE2_Boxplot_byPheno +
  ggtitle("CytoTRACE2 Potency Score by Cell Type") +
  xlab("Cell Type") +
  ylab("CytoTRACE2 Potency Score") +
  pub_theme +
  
  geom_signif(
    test = "wilcox.test",
    map_signif_level = TRUE,
    textsize = 3.8,
    vjust = 0.4
  )


# if(!is.null(plots$CytoTRACE2_Gene_UMAP)){
#   p_gene <- plots$CytoTRACE2_Gene_UMAP +
#     ggtitle("CytoTRACE2 Gene Signature") +
#     pub_theme
#   print(p_gene)
# }
if(!is.null(plots$CytoTRACE2_Density)){
  p_density <- plots$CytoTRACE2_Density +
    ggtitle("Potency Score Distribution") +
    pub_theme
  print(p_density)
}
ggsave("CytoTRACE2_UMAP.png", p_umap, width = 8, height = 6, dpi = 300)
ggsave("CytoTRACE2_Boxplot2.pdf", p_box, width = 10, height = 6, dpi = 300)
