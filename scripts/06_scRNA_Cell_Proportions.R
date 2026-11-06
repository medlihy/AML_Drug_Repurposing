############################################
## SCI-level single-cell cell proportion plots
############################################

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
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(readr)
library(tibble)
library(forcats)
library(grid)

############################################
############################################

ANNO_DIR <- SCRNA_ANNOTATION_DIR
RATIO_DIR <- SCRNA_CELL_PROPORTION_DIR

if (!dir.exists(RATIO_DIR)) dir.create(RATIO_DIR, recursive = TRUE)

setwd(ANNO_DIR)
sce <- qread(SCRNA_ANNOTATED_OBJECT_FILE)

setwd(RATIO_DIR)

############################################
############################################

meta <- sce@meta.data %>%
  rownames_to_column("CellID") %>%
  mutate(
    Sample   = as.character(orig.ident),
    Type_old = as.character(type),
    CellType = as.character(celltype)
  )

if (all(meta$Sample == meta$Type_old)) {
  message("注意：orig.ident 和 type 完全一致。后续使用 orig.ident 作为 Sample，并根据样本名前缀生成 Disease。")
} else {
  message("注意：orig.ident 和 type 不完全一致，请人工确认哪个字段代表样本，哪个字段代表分组。")
}

meta <- meta %>%
  mutate(
    Disease = case_when(
      str_detect(Sample, "^AML") ~ "AML",
      str_detect(Sample, "^N")   ~ "Normal",
      TRUE ~ "Other"
    )
  )

meta$Disease <- factor(meta$Disease, levels = c("Normal", "AML", "Other"))

############################################
############################################

sample_order <- c(
  str_sort(unique(meta$Sample[meta$Disease == "Normal"]), numeric = TRUE),
  str_sort(unique(meta$Sample[meta$Disease == "AML"]), numeric = TRUE),
  str_sort(unique(meta$Sample[meta$Disease == "Other"]), numeric = TRUE)
)

meta$Sample <- factor(meta$Sample, levels = sample_order)

############################################
############################################

celltype_order_auto <- meta %>%
  count(CellType, sort = TRUE) %>%
  pull(CellType)

celltype_order_manual <- c(
  "AML_Blast",
  "AML_Progenitor",
  "MonoBlast",
  "Neutrophil",
  "Plasma",
  "Eosinophil",
  "Erythroid",
  "IFN-Myelocyte",
  "iMono",
  "Macrophage",
  "Myelocyte",
  "ncMono",
  "Neu&Ret",
  "NK",
  "NK_T",
  "ProB",
  "EBP",
  "ProEosinophil",
  "ProEry",
  "ProMono",
  "Promyelocyte"
)

celltype_order <- celltype_order_manual[celltype_order_manual %in% unique(meta$CellType)]

if (length(celltype_order) < 2) {
  celltype_order <- celltype_order_auto
}

celltype_order <- c(
  celltype_order,
  setdiff(celltype_order_auto, celltype_order)
)

meta$CellType <- factor(meta$CellType, levels = celltype_order)

############################################
############################################

manual_colors <- c(
  "AML_Blast"      = "#815e99",
  "AML_Progenitor" = "#7B3294",
  "MonoBlast"      = "#d4a6a8",
  "Neutrophil"     = "#F4A582",
  "Plasma"         = "#E08214",
  "Eosinophil"     = "#FDB863",
  "Erythroid"      = "#DFC27D",
  "IFN-Myelocyte"  = "#2166AC",
  "iMono"          = "#4393C3",
  "Macrophage"     = "#f9b769",
  "Myelocyte"      = "#ec5051",
  "ncMono"         = "#C51B7D",
  "Neu&Ret"        = "#9970AB",
  "NK"             = "#f6f28f",
  "NK_T"           = "#a4cde1",
  "ProB"           = "#67a4cc",
  "EBP"            = "#277fb8",
  "ProEosinophil"  = "#549da3",
  "ProEry"         = "#96cb8f",
  "ProMono"        = "#5c9e43",
  "Promyelocyte"   = "#4dae47"
)

missing_celltypes <- setdiff(celltype_order, names(manual_colors))

auto_colors <- grDevices::hcl.colors(
  n = length(missing_celltypes),
  palette = "Dark 3"
)

names(auto_colors) <- missing_celltypes

celltype_colors <- c(
  manual_colors[names(manual_colors) %in% celltype_order],
  auto_colors
)

celltype_colors <- celltype_colors[celltype_order]

disease_colors <- c(
  Normal = "#4D9221",
  AML    = "#B2182B",
  Other  = "#999999"
)

############################################
############################################

theme_sci <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
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
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(size = base_size, face = "bold", color = "black")
    )
}

############################################
############################################

cell_count <- meta %>%
  count(Sample, Disease, CellType, name = "CellNumber") %>%
  group_by(Sample) %>%
  mutate(
    TotalCell = sum(CellNumber),
    Proportion = CellNumber / TotalCell,
    Percent = Proportion * 100
  ) %>%
  ungroup()

write_csv(cell_count, "Supplementary_celltype_count_and_proportion_by_sample.csv")

############################################
############################################

group_mean <- cell_count %>%
  group_by(Disease, CellType) %>%
  summarise(
    MeanProportion = mean(Proportion),
    SDProportion   = sd(Proportion),
    MeanPercent    = mean(Percent),
    SDPercent      = sd(Percent),
    .groups = "drop"
  )

write_csv(group_mean, "Supplementary_celltype_mean_proportion_by_group.csv")

############################################
############################################

p_stacked <- ggplot(cell_count, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.15
  ) +
  facet_grid(. ~ Disease, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = NULL,
    y = "Cell proportion",
    fill = "Cell type"
  ) +
  theme_sci(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    panel.spacing.x = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.justification = "center",
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.key.width = unit(0.36, "cm"),
    legend.key.height = unit(0.28, "cm"),
    legend.spacing.x = unit(0.28, "cm"),
    legend.spacing.y = unit(0.08, "cm"),
    
    legend.box.margin = margin(t = 2, r = 2, b = 4, l = 2),
    plot.margin = margin(t = 6, r = 8, b = 12, l = 8)
  ) +
  guides(
    fill = guide_legend(
      ncol = 7,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0
    )
  )

ggsave(
  filename = "Fig_celltype_composition_stacked_by_sample.pdf",
  plot = p_stacked,
  width = 10.2,
  height = 6.2,
  device = cairo_pdf
)

ggsave(
  filename = "Fig_celltype_composition_stacked_by_sample.png",
  plot = p_stacked,
  width = 10.2,
  height = 6.2,
  dpi = 600
)

############################################
############################################

p_heatmap <- ggplot(cell_count, aes(x = Sample, y = CellType, fill = Proportion)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(. ~ Disease, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low = "grey95",
    high = "#B2182B",
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Proportion"
  ) +
  theme_sci(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    panel.spacing.x = unit(0.8, "lines")
  )

ggsave(
  filename = "Fig_celltype_composition_heatmap_by_sample.pdf",
  plot = p_heatmap,
  width = 8.5,
  height = max(4.2, 0.28 * length(celltype_order) + 2),
  device = cairo_pdf
)

ggsave(
  filename = "Fig_celltype_composition_heatmap_by_sample.png",
  plot = p_heatmap,
  width = 8.5,
  height = max(4.2, 0.28 * length(celltype_order) + 2),
  dpi = 600
)

############################################
############################################

key_celltypes <- celltype_order[1:min(8, length(celltype_order))]

## key_celltypes <- c("AML_Blast", "AML_Progenitor", "MonoBlast", "T_cell", "NK_cell", "B_cell")

plot_box <- cell_count %>%
  filter(CellType %in% key_celltypes) %>%
  mutate(CellType = factor(CellType, levels = key_celltypes))

p_box <- ggplot(plot_box, aes(x = Disease, y = Proportion, color = Disease)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.45,
    fill = "white"
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.8,
    alpha = 0.85
  ) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
  scale_color_manual(values = disease_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Cell proportion",
    color = "Group"
  ) +
  theme_sci(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 9)
  )

ggsave(
  filename = "Fig_key_celltype_proportion_boxplot_by_group.pdf",
  plot = p_box,
  width = 8,
  height = 5.6,
  device = cairo_pdf
)

ggsave(
  filename = "Fig_key_celltype_proportion_boxplot_by_group.png",
  plot = p_box,
  width = 8,
  height = 5.6,
  dpi = 600
)

############################################
############################################

p_bubble <- ggplot(cell_count, aes(x = Sample, y = CellType)) +
  geom_point(
    aes(size = Proportion, fill = Proportion),
    shape = 21,
    color = "grey25",
    stroke = 0.2
  ) +
  facet_grid(. ~ Disease, scales = "free_x", space = "free_x") +
  scale_size_continuous(
    range = c(0.2, 6),
    labels = percent_format(accuracy = 1)
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "#B2182B",
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Proportion",
    fill = "Proportion"
  ) +
  theme_sci(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    panel.spacing.x = unit(0.8, "lines")
  )

ggsave(
  filename = "Supplementary_celltype_proportion_bubbleplot_by_sample.pdf",
  plot = p_bubble,
  width = 9,
  height = max(4.2, 0.28 * length(celltype_order) + 2),
  device = cairo_pdf
)

ggsave(
  filename = "Supplementary_celltype_proportion_bubbleplot_by_sample.png",
  plot = p_bubble,
  width = 9,
  height = max(4.2, 0.28 * length(celltype_order) + 2),
  dpi = 600
)

############################################
############################################

overall_prop <- meta %>%
  count(CellType, name = "CellNumber") %>%
  mutate(
    Proportion = CellNumber / sum(CellNumber),
    Percent = Proportion * 100
  ) %>%
  arrange(desc(Proportion)) %>%
  mutate(CellType = factor(CellType, levels = rev(CellType)))

write_csv(overall_prop, "Supplementary_overall_celltype_proportion.csv")

p_overall <- ggplot(overall_prop, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_col(
    width = 0.7,
    color = "white",
    linewidth = 0.2
  ) +
  coord_flip() +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Overall cell proportion",
    fill = "Cell type"
  ) +
  theme_sci(base_size = 10) +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "Supplementary_overall_celltype_proportion_barplot.pdf",
  plot = p_overall,
  width = 5.8,
  height = max(4.2, 0.28 * length(celltype_order) + 2),
  device = cairo_pdf
)

ggsave(
  filename = "Supplementary_overall_celltype_proportion_barplot.png",
  plot = p_overall,
  width = 5.8,
  height = max(4.2, 0.28 * length(celltype_order) + 2),
  dpi = 600
)

############################################
############################################

sample_summary <- meta %>%
  count(Sample, Disease, name = "TotalCell") %>%
  arrange(Disease, Sample)

write_csv(sample_summary, "Supplementary_sample_total_cell_number.csv")

############################################
############################################

message("绘图完成。主要输出文件：")
message("1. Fig_celltype_composition_stacked_by_sample.pdf")
message("2. Fig_celltype_composition_heatmap_by_sample.pdf")
message("3. Fig_key_celltype_proportion_boxplot_by_group.pdf")
message("4. Supplementary_celltype_proportion_bubbleplot_by_sample.pdf")
message("5. Supplementary_overall_celltype_proportion_barplot.pdf")
