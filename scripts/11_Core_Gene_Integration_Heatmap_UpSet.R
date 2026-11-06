############################################################
## AML strict 62-gene core screening + SCI circular heatmap
## + center UpSet
##
## One-step integrated pipeline
##
## PART A:
##   Strict bulk/scRNA/WGCNA/hdWGCNA core-gene screening
##   Threshold sensitivity analysis
##   Core tables, CMap inputs, QC and auxiliary figures
##   Strict UpSet gene-set export
##
## PART B:
##   Directly read the just-generated final strict core result
##   Verify exact 62-gene consistency
##   Draw directional evidence and network-strength tracks
##   Draw enlarged center UpSet with intersection counts
##   Highlight Final-core-related intersections
##
## IMPORTANT:
##   Run this script from top to bottom.
##   Do not run PART B alone unless PART A outputs already exist.
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
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexUpset)
})

############################################################
############################################################

BASE_DIR <- CORE_GENE_RESULTS_DIR
dir.create(BASE_DIR, recursive = TRUE, showWarnings = FALSE)
bulk_deg_file <- CORE_BULK_DEG_FILE
bulk_wgcna_gene_file <- CORE_BULK_WGCNA_GENE_FILE
sc_deg_file <- CORE_SCRNA_DEG_FILE
sc_module_gene_file <- CORE_SCRNA_MODULE_FILE
sc_hub_gene_file <- CORE_SCRNA_HUB_FILE
sc_expr_pct_file <- CORE_SCRNA_EXPRESSION_PCT_FILE

############################################################
############################################################

bulk_logFC_cutoff <- 0.50
bulk_fdr_cutoff <- 0.05

sc_logFC_cutoff <- 0.30
sc_fdr_cutoff <- 0.05

pct_cutoff <- 0.05

bulk_kME_cutoff <- 0.40
bulk_GS_cutoff <- 0.20
sc_kME_cutoff <- 0.15

############################################################
## bulk WGCNA AML-associated modules
############################################################

bulk_module_trait <- data.frame(
  bulk_module = c("black", "turquoise", "magenta", "blue", "brown"),
  bulk_module_cor_AML = c(-0.54, -0.47, 0.31, 0.52, 0.52),
  bulk_module_p_AML = c(5.8e-17, 5.3e-13, 5.9e-06, 1.2e-15, 8.1e-16),
  bulk_module_direction = c(
    "AML_negative",
    "AML_negative",
    "AML_positive",
    "AML_positive",
    "AML_positive"
  ),
  stringsAsFactors = FALSE
)

bulk_disease_modules <- bulk_module_trait$bulk_module

bulk_AML_positive_modules <- bulk_module_trait %>%
  filter(bulk_module_cor_AML > 0) %>%
  pull(bulk_module)

bulk_AML_negative_modules <- bulk_module_trait %>%
  filter(bulk_module_cor_AML < 0) %>%
  pull(bulk_module)

cat("\nSelected bulk AML-associated WGCNA modules:\n")
print(bulk_module_trait)

if (length(bulk_disease_modules) == 0) {
  stop("必须填写 bulk_disease_modules，即 bulk WGCNA 中与 AML status 显著相关的模块。")
}

############################################################
############################################################

check_file <- function(x) {
  if (!file.exists(x)) stop("找不到文件：", x)
}

invisible(lapply(
  c(
    bulk_deg_file,
    bulk_wgcna_gene_file,
    sc_deg_file,
    sc_module_gene_file,
    sc_hub_gene_file
  ),
  check_file
))

############################################################
############################################################

bulk_deg <- fread(bulk_deg_file, data.table = FALSE)

if ("gene_symbol" %in% colnames(bulk_deg)) {
  bulk_deg$gene <- bulk_deg$gene_symbol
} else if (!"gene" %in% colnames(bulk_deg)) {
  stop("bulk DEG 表中找不到 gene_symbol 或 gene 列")
}

bulk_deg2 <- bulk_deg %>%
  transmute(
    gene = as.character(gene),
    bulk_logFC = as.numeric(logFC),
    bulk_P.Value = as.numeric(P.Value),
    bulk_adj.P.Val = as.numeric(adj.P.Val),
    bulk_DEG_pass = bulk_adj.P.Val < bulk_fdr_cutoff &
      abs(bulk_logFC) >= bulk_logFC_cutoff,
    bulk_direction = case_when(
      bulk_DEG_pass & bulk_logFC > 0 ~ "bulk_up",
      bulk_DEG_pass & bulk_logFC < 0 ~ "bulk_down",
      TRUE ~ "bulk_not_sig"
    )
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene, .keep_all = TRUE)

cat("\nBulk DEG direction:\n")
print(table(bulk_deg2$bulk_direction))

############################################################
############################################################

bulk_wgcna <- fread(bulk_wgcna_gene_file, data.table = FALSE)

if ("gene_symbol" %in% colnames(bulk_wgcna)) {
  bulk_wgcna$gene <- bulk_wgcna$gene_symbol
} else if (!"gene" %in% colnames(bulk_wgcna)) {
  stop("bulk WGCNA 表中找不到 gene_symbol 或 gene 列")
}

bulk_wgcna2 <- bulk_wgcna %>%
  transmute(
    gene = as.character(gene),
    bulk_module = as.character(module),
    bulk_kME = if ("kME_own_module" %in% colnames(bulk_wgcna)) {
      as.numeric(kME_own_module)
    } else {
      NA_real_
    },
    bulk_GS_AML = if ("GS_AML_status" %in% colnames(bulk_wgcna)) {
      as.numeric(GS_AML_status)
    } else {
      NA_real_
    },
    bulk_hub_score = if ("hub_score" %in% colnames(bulk_wgcna)) {
      as.numeric(hub_score)
    } else {
      NA_real_
    }
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene, .keep_all = TRUE)
bulk_wgcna2 <- bulk_wgcna2 %>%
  left_join(bulk_module_trait, by = "bulk_module")

module_qc <- bulk_wgcna2 %>%
  filter(!is.na(bulk_module), bulk_module != "grey") %>%
  group_by(bulk_module) %>%
  summarise(
    n_genes = n(),
    median_abs_GS = median(abs(bulk_GS_AML), na.rm = TRUE),
    mean_abs_GS = mean(abs(bulk_GS_AML), na.rm = TRUE),
    median_abs_kME = median(abs(bulk_kME), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_abs_GS))

fwrite(
  module_qc,
  file.path(BASE_DIR, "00_bulk_WGCNA_module_QC_check.tsv"),
  sep = "\t"
)

cat("\nBulk disease modules used:\n")
print(bulk_disease_modules)

cat("\nBulk WGCNA module QC saved:\n")
cat(file.path(BASE_DIR, "00_bulk_WGCNA_module_QC_check.tsv"), "\n")

############################################################
############################################################

sc_deg <- fread(sc_deg_file, data.table = FALSE)

if (!"gene" %in% colnames(sc_deg)) {
  stop("scRNA DEG 表中找不到 gene 列")
}

sc_deg2 <- sc_deg %>%
  transmute(
    gene = as.character(gene),
    sc_logFC = as.numeric(logFC),
    sc_P.Value = as.numeric(P.Value),
    sc_adj.P.Val = as.numeric(adj.P.Val),
    sc_DEG_pass = sc_adj.P.Val < sc_fdr_cutoff &
      abs(sc_logFC) >= sc_logFC_cutoff,
    sc_direction = case_when(
      sc_DEG_pass & sc_logFC > 0 ~ "sc_up",
      sc_DEG_pass & sc_logFC < 0 ~ "sc_down",
      TRUE ~ "sc_not_sig"
    )
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene, .keep_all = TRUE)

cat("\nscRNA DEG direction:\n")
print(table(sc_deg2$sc_direction))

############################################################
############################################################

get_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

pct_malignant_col <- get_first_existing_col(
  sc_deg,
  c("pct_malignant", "pct_AML_malignant", "pct_aml_malignant", "pct.1", "pct_1")
)

pct_normal_col <- get_first_existing_col(
  sc_deg,
  c("pct_normal", "pct_NormalBM", "pct_normal_bm", "pct.2", "pct_2")
)

if (!is.na(pct_malignant_col) && !is.na(pct_normal_col)) {
  
  sc_pct <- sc_deg %>%
    transmute(
      gene = as.character(gene),
      pct_malignant = as.numeric(.data[[pct_malignant_col]]),
      pct_normal = as.numeric(.data[[pct_normal_col]])
    ) %>%
    filter(!is.na(gene), gene != "") %>%
    distinct(gene, .keep_all = TRUE)
  
  cat("\n使用 sc DEG 表中的表达比例列：\n")
  cat("pct_malignant:", pct_malignant_col, "\n")
  cat("pct_normal:", pct_normal_col, "\n")
  
} else if (nzchar(sc_expr_pct_file) && file.exists(sc_expr_pct_file)) {
  
  sc_pct <- fread(sc_expr_pct_file, data.table = FALSE)
  
  if (!all(c("gene", "pct_malignant", "pct_normal") %in% colnames(sc_pct))) {
    stop("sc_expr_pct_file 必须包含列：gene, pct_malignant, pct_normal")
  }
  
  sc_pct <- sc_pct %>%
    transmute(
      gene = as.character(gene),
      pct_malignant = as.numeric(pct_malignant),
      pct_normal = as.numeric(pct_normal)
    ) %>%
    filter(!is.na(gene), gene != "") %>%
    distinct(gene, .keep_all = TRUE)
  
  cat("\n使用外部 scRNA 表达比例文件：\n")
  cat(sc_expr_pct_file, "\n")
  
} else {
  
  stop(
    "\n没有找到 scRNA 表达比例信息。\n",
    "你需要满足以下任一条件：\n",
    "1. sc DEG 表中含有 pct.1 和 pct.2，且 pct.1=AML malignant，pct.2=Normal BM；\n",
    "2. sc DEG 表中含有 pct_malignant 和 pct_normal；\n",
    "3. 提供 sc_expr_pct_file，列名为 gene, pct_malignant, pct_normal。\n"
  )
}

if (max(sc_pct$pct_malignant, na.rm = TRUE) > 1.5 ||
    max(sc_pct$pct_normal, na.rm = TRUE) > 1.5) {
  sc_pct <- sc_pct %>%
    mutate(
      pct_malignant = pct_malignant / 100,
      pct_normal = pct_normal / 100
    )
}

############################################################
############################################################

sc_module_gene <- fread(sc_module_gene_file, data.table = FALSE)
sc_hub_gene <- fread(sc_hub_gene_file, data.table = FALSE)

if (!all(c("gene_name", "module") %in% colnames(sc_module_gene))) {
  stop("Tab01_moduleGene.csv 中需要 gene_name 和 module 列")
}

if (!all(c("gene_name", "module", "kME") %in% colnames(sc_hub_gene))) {
  stop("Tab02_hubGene.csv 中需要 gene_name、module 和 kME 列")
}

sc_module2 <- sc_module_gene %>%
  transmute(
    gene = as.character(gene_name),
    sc_hdWGCNA_module = as.character(module)
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(sc_hdWGCNA_module),
    sc_hdWGCNA_module != "grey"
  ) %>%
  distinct(gene, .keep_all = TRUE)

sc_hub2 <- sc_hub_gene %>%
  transmute(
    gene = as.character(gene_name),
    sc_hub_module = as.character(module),
    sc_hub_kME = as.numeric(kME)
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene, .keep_all = TRUE)

sc_own_kME <- vapply(seq_len(nrow(sc_module_gene)), function(i) {
  mod <- as.character(sc_module_gene$module[i])
  kme_col <- paste0("kME_", mod)
  
  if (!kme_col %in% colnames(sc_module_gene)) {
    return(NA_real_)
  }
  
  as.numeric(sc_module_gene[[kme_col]][i])
}, numeric(1))

sc_own_kME_df <- sc_module_gene %>%
  mutate(sc_own_kME = sc_own_kME) %>%
  transmute(
    gene = as.character(gene_name),
    sc_hdWGCNA_module = as.character(module),
    sc_own_kME = as.numeric(sc_own_kME)
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(sc_hdWGCNA_module),
    sc_hdWGCNA_module != "grey"
  ) %>%
  distinct(gene, sc_hdWGCNA_module, .keep_all = TRUE)

############################################################
############################################################

rank01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- median(x, na.rm = TRUE)
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    return(rep(0, length(x)))
  }
  as.numeric((rank(x, ties.method = "average") - 1) / (length(x) - 1))
}

core_all <- bulk_deg2 %>%
  inner_join(sc_deg2, by = "gene") %>%
  left_join(bulk_wgcna2, by = "gene") %>%
  left_join(sc_module2, by = "gene") %>%
  left_join(sc_hub2, by = "gene") %>%
  left_join(sc_own_kME_df, by = c("gene", "sc_hdWGCNA_module")) %>%
  left_join(sc_pct, by = "gene") %>%
  mutate(
    sc_kME_final = dplyr::coalesce(sc_own_kME, sc_hub_kME),
    
    direction_consistency = case_when(
      bulk_logFC > 0 & sc_logFC > 0 ~ "same_up",
      bulk_logFC < 0 & sc_logFC < 0 ~ "same_down",
      TRUE ~ "discordant"
    ),
    
    core_direction = case_when(
      bulk_DEG_pass &
        sc_DEG_pass &
        bulk_logFC > 0 &
        sc_logFC > 0 ~ "core_up",
      
      bulk_DEG_pass &
        sc_DEG_pass &
        bulk_logFC < 0 &
        sc_logFC < 0 ~ "core_down",
      
      TRUE ~ "not_core"
    ),
    
    in_bulk_AML_WGCNA_module = !is.na(bulk_module) &
      bulk_module %in% bulk_disease_modules,
    
    bulk_WGCNA_direction_pass = case_when(
      core_direction == "core_up" ~ bulk_module %in% bulk_AML_positive_modules,
      core_direction == "core_down" ~ bulk_module %in% bulk_AML_negative_modules,
      TRUE ~ FALSE
    ),
    
    in_sc_hdWGCNA = !is.na(sc_hdWGCNA_module),
    
    expression_prevalence_pass = case_when(
      core_direction == "core_up" ~ !is.na(pct_malignant) &
        pct_malignant >= pct_cutoff,
      
      core_direction == "core_down" ~ !is.na(pct_normal) &
        pct_normal >= pct_cutoff,
      
      TRUE ~ FALSE
    ),
    
    network_centrality_pass =
      !is.na(bulk_kME) &
      !is.na(bulk_GS_AML) &
      !is.na(sc_kME_final) &
      abs(bulk_kME) >= bulk_kME_cutoff &
      abs(bulk_GS_AML) >= bulk_GS_cutoff &
      abs(sc_kME_final) >= sc_kME_cutoff,
    
    final_core =
      core_direction %in% c("core_up", "core_down") &
      in_bulk_AML_WGCNA_module &
      bulk_WGCNA_direction_pass &
      in_sc_hdWGCNA &
      expression_prevalence_pass &
      network_centrality_pass,
    
    score_bulk_fdr = -log10(pmax(bulk_adj.P.Val, 1e-300)),
    score_sc_fdr = -log10(pmax(sc_adj.P.Val, 1e-300)),
    score_bulk_logFC = abs(bulk_logFC),
    score_sc_logFC = abs(sc_logFC),
    score_bulk_kME = abs(bulk_kME),
    score_bulk_GS = abs(bulk_GS_AML),
    score_sc_kME = abs(sc_kME_final)
  )

core_all$evidence_score <-
  rank01(core_all$score_bulk_fdr) +
  rank01(core_all$score_sc_fdr) +
  rank01(core_all$score_bulk_logFC) +
  rank01(core_all$score_sc_logFC) +
  rank01(core_all$score_bulk_kME) +
  rank01(core_all$score_bulk_GS) +
  rank01(core_all$score_sc_kME)

fwrite(
  core_all,
  file.path(BASE_DIR, "01_integrated_bulk_scRNA_WGCNA_gene_table_v2.tsv"),
  sep = "\t"
)

############################################################
## 7.1 Threshold sensitivity analysis
############################################################

make_core_by_threshold <- function(df,
                                   bulk_logFC_cut,
                                   sc_logFC_cut,
                                   pct_cut,
                                   bulk_kME_cut,
                                   bulk_GS_cut,
                                   sc_kME_cut) {
  
  tmp <- df %>%
    mutate(
      bulk_DEG_pass2 = bulk_adj.P.Val < bulk_fdr_cutoff &
        abs(bulk_logFC) >= bulk_logFC_cut,
      
      sc_DEG_pass2 = sc_adj.P.Val < sc_fdr_cutoff &
        abs(sc_logFC) >= sc_logFC_cut,
      
      core_direction2 = case_when(
        bulk_DEG_pass2 & sc_DEG_pass2 & bulk_logFC > 0 & sc_logFC > 0 ~ "core_up",
        bulk_DEG_pass2 & sc_DEG_pass2 & bulk_logFC < 0 & sc_logFC < 0 ~ "core_down",
        TRUE ~ "not_core"
      ),
      
      expression_prevalence_pass2 = case_when(
        core_direction2 == "core_up" ~ !is.na(pct_malignant) &
          pct_malignant >= pct_cut,
        core_direction2 == "core_down" ~ !is.na(pct_normal) &
          pct_normal >= pct_cut,
        TRUE ~ FALSE
      ),
      
      bulk_WGCNA_direction_pass2 = case_when(
        core_direction2 == "core_up" ~ bulk_module %in% bulk_AML_positive_modules,
        core_direction2 == "core_down" ~ bulk_module %in% bulk_AML_negative_modules,
        TRUE ~ FALSE
      ),
      
      network_centrality_pass2 =
        !is.na(bulk_kME) &
        !is.na(bulk_GS_AML) &
        !is.na(sc_kME_final) &
        abs(bulk_kME) >= bulk_kME_cut &
        abs(bulk_GS_AML) >= bulk_GS_cut &
        abs(sc_kME_final) >= sc_kME_cut,
      
      final_core2 =
        core_direction2 %in% c("core_up", "core_down") &
        in_bulk_AML_WGCNA_module &
        bulk_WGCNA_direction_pass2 &
        in_sc_hdWGCNA &
        expression_prevalence_pass2 &
        network_centrality_pass2
    )
  
  tmp %>% filter(final_core2)
}

threshold_grid <- expand.grid(
  bulk_logFC_cut = c(0.50, 0.75),
  sc_logFC_cut   = c(0.30, 0.40),
  pct_cut        = c(0.05, 0.10),
  bulk_kME_cut   = c(0.40, 0.45),
  bulk_GS_cut    = c(0.15, 0.20),
  sc_kME_cut     = c(0.10, 0.15, 0.20),
  stringsAsFactors = FALSE
)

threshold_sensitivity <- lapply(seq_len(nrow(threshold_grid)), function(i) {
  
  p <- threshold_grid[i, ]
  
  tmp_core <- make_core_by_threshold(
    df = core_all,
    bulk_logFC_cut = p$bulk_logFC_cut,
    sc_logFC_cut   = p$sc_logFC_cut,
    pct_cut        = p$pct_cut,
    bulk_kME_cut   = p$bulk_kME_cut,
    bulk_GS_cut    = p$bulk_GS_cut,
    sc_kME_cut     = p$sc_kME_cut
  )
  
  data.frame(
    p,
    n_total = nrow(tmp_core),
    n_up = sum(tmp_core$core_direction2 == "core_up"),
    n_down = sum(tmp_core$core_direction2 == "core_down")
  )
}) %>%
  bind_rows() %>%
  arrange(desc(n_total))

fwrite(
  threshold_sensitivity,
  file.path(BASE_DIR, "00_threshold_sensitivity_analysis.tsv"),
  sep = "\t"
)

print(threshold_sensitivity)


############################################################
############################################################

final_core <- core_all %>%
  filter(final_core) %>%
  arrange(core_direction, desc(evidence_score))

fwrite(
  final_core,
  file.path(BASE_DIR, "02_Final_core_signature_v2.tsv"),
  sep = "\t"
)

cat("\nFinal core genes v2:\n")
cat("Total:", nrow(final_core), "\n")
print(table(final_core$core_direction))

cat("\nFinal core by bulk module:\n")
print(table(final_core$bulk_module, final_core$core_direction))

cat("\nFinal core by sc hdWGCNA module:\n")
print(table(final_core$sc_hdWGCNA_module, final_core$core_direction))

CMAP_N <- 150

core_up <- final_core %>%
  filter(core_direction == "core_up") %>%
  arrange(desc(evidence_score)) %>%
  slice_head(n = CMAP_N) %>%
  pull(gene) %>%
  unique()

core_down <- final_core %>%
  filter(core_direction == "core_down") %>%
  arrange(desc(evidence_score)) %>%
  slice_head(n = CMAP_N) %>%
  pull(gene) %>%
  unique()

writeLines(
  core_up,
  file.path(BASE_DIR, "CMap_Final_core_v2_UP_top150.txt")
)

writeLines(
  core_down,
  file.path(BASE_DIR, "CMap_Final_core_v2_DOWN_top150.txt")
)

cat("\nCMap input genes v2:\n")
cat("UP:", length(core_up), "\n")
cat("DOWN:", length(core_down), "\n")

############################################################
## 8.2 Evidence matrix for top final core genes
############################################################

plot_core_top <- final_core %>%
  group_by(core_direction) %>%
  arrange(desc(evidence_score), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup() %>%
  mutate(
    gene_label = paste0(gene, " | ", bulk_module, " | ", sc_hdWGCNA_module),
    gene_label = factor(gene_label, levels = rev(unique(gene_label))),
    direction_label = case_when(
      core_direction == "core_up" ~ "UP in AML",
      core_direction == "core_down" ~ "DOWN in AML",
      TRUE ~ core_direction
    )
  )

evidence_long <- plot_core_top %>%
  transmute(
    gene_label,
    direction_label,
    `Bulk |log2FC|` = abs(bulk_logFC),
    `scRNA |log2FC|` = abs(sc_logFC),
    `Bulk |kME|` = abs(bulk_kME),
    `Bulk |GS|` = abs(bulk_GS_AML),
    `sc |kME|` = abs(sc_kME_final),
    `Expression pct` = ifelse(core_direction == "core_up", pct_malignant, pct_normal),
    `Evidence score` = evidence_score
  ) %>%
  pivot_longer(
    cols = -c(gene_label, direction_label),
    names_to = "evidence",
    values_to = "value"
  ) %>%
  group_by(evidence) %>%
  mutate(value_scaled = rank01(value)) %>%
  ungroup()

p_evidence_heatmap <- ggplot(
  evidence_long,
  aes(x = evidence, y = gene_label, fill = value_scaled)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_grid(
    direction_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#B2182B",
    name = "Relative\nstrength"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text.y = element_text(face = "bold", color = "#222222"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "#222222"),
    axis.text.y = element_text(size = 8.5, color = "#222222"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  filename = file.path(BASE_DIR, "Fig_Final_core_evidence_matrix.pdf"),
  plot = p_evidence_heatmap,
  width = 8.2,
  height = 8.8
)

ggsave(
  filename = file.path(BASE_DIR, "Fig_Final_core_evidence_matrix.png"),
  plot = p_evidence_heatmap,
  width = 8.2,
  height = 8.8,
  dpi = 600
)

############################################################
## 9. Filtering diagnostics
############################################################

diagnostic_df <- data.frame(
  step = c(
    "Integrated genes",
    "Bulk DEG pass",
    "scRNA DEG pass",
    "Bulk + scRNA concordant DEG",
    "AML-associated bulk WGCNA module",
    "Bulk WGCNA direction matched",
    "sc hdWGCNA module",
    "Expression prevalence pass",
    "Network centrality pass",
    "Final core"
  ),
  n = c(
    nrow(core_all),
    sum(core_all$bulk_DEG_pass, na.rm = TRUE),
    sum(core_all$sc_DEG_pass, na.rm = TRUE),
    sum(core_all$core_direction %in% c("core_up", "core_down"), na.rm = TRUE),
    sum(
      core_all$core_direction %in% c("core_up", "core_down") &
        core_all$in_bulk_AML_WGCNA_module,
      na.rm = TRUE
    ),
    sum(
      core_all$core_direction %in% c("core_up", "core_down") &
        core_all$in_bulk_AML_WGCNA_module &
        core_all$bulk_WGCNA_direction_pass,
      na.rm = TRUE
    ),
    sum(
      core_all$core_direction %in% c("core_up", "core_down") &
        core_all$in_bulk_AML_WGCNA_module &
        core_all$bulk_WGCNA_direction_pass &
        core_all$in_sc_hdWGCNA,
      na.rm = TRUE
    ),
    sum(
      core_all$core_direction %in% c("core_up", "core_down") &
        core_all$in_bulk_AML_WGCNA_module &
        core_all$bulk_WGCNA_direction_pass &
        core_all$in_sc_hdWGCNA &
        core_all$expression_prevalence_pass,
      na.rm = TRUE
    ),
    sum(
      core_all$core_direction %in% c("core_up", "core_down") &
        core_all$in_bulk_AML_WGCNA_module &
        core_all$bulk_WGCNA_direction_pass &
        core_all$in_sc_hdWGCNA &
        core_all$expression_prevalence_pass &
        core_all$network_centrality_pass,
      na.rm = TRUE
    ),
    nrow(final_core)
  )
)

print(diagnostic_df)

fwrite(
  diagnostic_df,
  file.path(BASE_DIR, "03_Final_core_filtering_diagnostics_v2.tsv"),
  sep = "\t"
)

############################################################
## 9.1 Publication-style filtering flow plot
############################################################

diagnostic_df$step <- factor(
  diagnostic_df$step,
  levels = diagnostic_df$step
)

p_filter_flow <- ggplot(diagnostic_df, aes(x = step, y = n)) +
  geom_col(width = 0.68, fill = "#4E79A7") +
  geom_text(
    aes(label = n),
    vjust = -0.35,
    size = 3.8,
    color = "#222222"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL,
    y = "Number of genes"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "#222222"
    ),
    axis.text.y = element_text(color = "#222222"),
    axis.title.y = element_text(face = "bold"),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  filename = file.path(BASE_DIR, "Fig_Final_core_filtering_flow.pdf"),
  plot = p_filter_flow,
  width = 8.5,
  height = 5.2
)

ggsave(
  filename = file.path(BASE_DIR, "Fig_Final_core_filtering_flow.png"),
  plot = p_filter_flow,
  width = 8.5,
  height = 5.2,
  dpi = 600
)

############################################################
############################################################

set_order <- c(
  "Bulk_DEG",
  "Bulk_AML_WGCNA",
  "scRNA_DEG",
  "sc_hdWGCNA",
  "Final_core"
)

upset_df <- core_all %>%
  transmute(
    gene,
    Bulk_DEG = bulk_DEG_pass,
    Bulk_AML_WGCNA = in_bulk_AML_WGCNA_module,
    scRNA_DEG = sc_DEG_pass,
    sc_hdWGCNA = in_sc_hdWGCNA,
    Final_core = final_core
  )

upset_df[set_order] <- lapply(upset_df[set_order], function(x) {
  x <- as.logical(x)
  x[is.na(x)] <- FALSE
  x
})

cat("\nUpSet set sizes:\n")
print(colSums(upset_df[set_order]))
############################################################
## 10.1 Export strict UpSet gene sets for circular heatmap
############################################################

strict_upset_groups <- lapply(set_order, function(ss) {
  unique(upset_df$gene[upset_df[[ss]]])
})

names(strict_upset_groups) <- set_order

strict_upset_groups$Final_core <- unique(final_core$gene)

strict_upset_rds <- file.path(BASE_DIR, "strict_upset_groups.rds")

saveRDS(strict_upset_groups, strict_upset_rds)

strict_upset_count <- data.table(
  set = names(strict_upset_groups),
  n_genes = sapply(strict_upset_groups, length)
)

fwrite(
  strict_upset_count,
  file.path(BASE_DIR, "strict_upset_groups_count.tsv"),
  sep = "\t"
)

fwrite(
  as.data.table(upset_df),
  file.path(BASE_DIR, "strict_upset_membership.tsv"),
  sep = "\t"
)

cat("\nExported strict UpSet groups:\n")
print(strict_upset_count)
cat("Saved RDS:\n", strict_upset_rds, "\n")


col_intersection <- "#5A9EC4"
col_set_size <- "#E9B3B8"
col_matrix <- "#3E8C8A"
col_core <- "#B2182B"
col_main_text <- "#222222"
col_grid <- "#EAEAEA"

p_upset <- upset(
  upset_df,
  intersect = set_order,
  name = "Gene set",
  width_ratio = 0.22,
  height_ratio = 0.60,
  n_intersections = 20,
  min_size = 1,
  sort_sets = FALSE,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  
  base_annotations = list(
    "Intersection size" = (
      intersection_size(
        counts = TRUE,
        text = list(
          size = 3.5,
          vjust = -0.25,
          color = col_main_text
        ),
        mapping = aes(fill = "Intersection")
      ) +
        scale_fill_manual(
          values = c("Intersection" = col_intersection),
          guide = "none"
        ) +
        ylab("Intersection size") +
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major.y = element_line(color = col_grid, linewidth = 0.35),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 11, face = "bold", color = col_main_text),
          axis.text.y = element_text(size = 9, color = col_main_text),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    )
  ),
  
  set_sizes = (
    upset_set_size(
      geom = geom_bar(
        fill = col_set_size,
        width = 0.68
      )
    ) +
      xlab("Set size") +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.x = element_line(color = col_grid, linewidth = 0.35),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10, face = "bold", color = col_main_text),
        axis.text.x = element_text(size = 9, color = col_main_text),
        axis.text.y = element_text(size = 10, color = col_main_text)
      )
  ),
  
  stripes = c("#FFFFFF", "#F8F8F8"),
  
  queries = list(
    upset_query(
      set = "Final_core",
      fill = col_core,
      color = col_core,
      only_components = c("Intersection size", "intersections_matrix")
    )
  ),
  
  themes = upset_modify_themes(
    list(
      "intersections_matrix" = theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 10, color = col_main_text),
        axis.title.x = element_text(size = 11, face = "bold", color = col_main_text),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ),
      "overall_sizes" = theme(
        axis.text.y = element_text(size = 10, color = col_main_text)
      )
    )
  )
)

ggsave(
  filename = file.path(BASE_DIR, "Fig_Final_core_UpSet.pdf"),
  plot = p_upset,
  width = 9.8,
  height = 6.0
)

ggsave(
  filename = file.path(BASE_DIR, "Fig_Final_core_UpSet.png"),
  plot = p_upset,
  width = 9.8,
  height = 6.0,
  dpi = 600
)

############################################################
## 11. Publication-style bulk logFC vs scRNA logFC scatter
############################################################

plot_df <- core_all %>%
  mutate(
    plot_group = case_when(
      final_core & core_direction == "core_up" ~ "Final core up",
      final_core & core_direction == "core_down" ~ "Final core down",
      bulk_DEG_pass & sc_DEG_pass & direction_consistency != "discordant" ~ "Concordant DEG",
      TRUE ~ "Other genes"
    ),
    plot_group = factor(
      plot_group,
      levels = c("Other genes", "Concordant DEG", "Final core down", "Final core up")
    )
  ) %>%
  arrange(plot_group)

label_df <- plot_df %>%
  filter(final_core) %>%
  group_by(core_direction) %>%
  arrange(desc(evidence_score), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

p_scatter <- ggplot(plot_df, aes(x = bulk_logFC, y = sc_logFC)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey55") +
  geom_point(
    data = subset(plot_df, plot_group == "Other genes"),
    aes(color = plot_group),
    alpha = 0.35,
    size = 0.9
  ) +
  geom_point(
    data = subset(plot_df, plot_group == "Concordant DEG"),
    aes(color = plot_group),
    alpha = 0.65,
    size = 1.2
  ) +
  geom_point(
    data = subset(plot_df, plot_group %in% c("Final core up", "Final core down")),
    aes(color = plot_group),
    alpha = 0.90,
    size = 2.0
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.1,
    color = "#222222",
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey60",
    segment.size = 0.25,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Other genes" = "grey88",
      "Concordant DEG" = "#BDBDBD",
      "Final core down" = "#2166AC",
      "Final core up" = "#B2182B"
    )
  ) +
  labs(
    x = "Bulk RNA log2FC: AML vs Healthy",
    y = "scRNA pseudo-bulk log2FC: AML malignant vs Normal BM",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold", color = "#222222"),
    axis.text = element_text(color = "#222222"),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  filename = file.path(BASE_DIR, "Fig_Final_core_bulk_sc_logFC_concordance_pub.pdf"),
  plot = p_scatter,
  width = 7.2,
  height = 6.0
)

ggsave(
  filename = file.path(BASE_DIR, "Fig_Final_core_bulk_sc_logFC_concordance_pub.png"),
  plot = p_scatter,
  width = 7.2,
  height = 6.0,
  dpi = 600
)

############################################################
############################################################

cat("\n==============================\n")
cat("Final core signature analysis completed.\n")
cat("Output directory:\n")
cat(BASE_DIR, "\n")
cat("==============================\n\n")

cat("Integrated genes:", nrow(core_all), "\n")
cat("Final core genes:", nrow(final_core), "\n")
cat("Final core UP:", length(core_up), "\n")
cat("Final core DOWN:", length(core_down), "\n\n")

cat("Main output files:\n")
cat("01_integrated_bulk_scRNA_WGCNA_gene_table_v2.tsv\n")
cat("02_Final_core_signature_v2.tsv\n")
cat("CMap_Final_core_v2_UP_top150.txt\n")
cat("CMap_Final_core_v2_DOWN_top150.txt\n")
cat("Fig_Final_core_UpSet.pdf/png\n")
cat("Fig_Final_core_bulk_sc_logFC_concordance_pub.pdf/png\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(ggplot2)
  library(cowplot)
})

############################################################
############################################################

FINAL_DIR <- BASE_DIR

CORE_FILE <- file.path(
  FINAL_DIR,
  "02_Final_core_signature_v2.tsv"
)

STRICT_UPSET_FILE <- file.path(
  FINAL_DIR,
  "strict_upset_groups.rds"
)

OUT_DIR <- file.path(
  dirname(FINAL_DIR),
  "Fig_main_circular_upset_core_SCI_final"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_PDF <- file.path(
  OUT_DIR,
  "Fig_main_circular_heatmap_with_upset_SCI_final.pdf"
)

OUT_PNG <- file.path(
  OUT_DIR,
  "Fig_main_circular_heatmap_with_upset_SCI_final.png"
)

SHOW_TITLE <- FALSE

FONT_FAMILY <- "sans"
options(bitmapType = "cairo")

############################################################
## Safe PNG device
##
## Your macOS/R environment may not have XQuartz/Cairo libraries.
## Therefore, the final 600-dpi PNG is exported with ragg instead of:
## png(..., type = "cairo")
##
## PDF export is still the primary SCI-ready vector output.
############################################################

export_png_safe <- function(
    filename,
    width = 14.0,
    height = 15.0,
    res = 600,
    draw_fun
) {
  
  if (!requireNamespace("ragg", quietly = TRUE)) {
  warning(
    "Package 'ragg' is unavailable. PNG export is skipped; ",
    "the PDF remains the primary vector output."
  )
  return(invisible(FALSE))
}
  
  if (!requireNamespace("ragg", quietly = TRUE)) {
    
    warning(
      "Package 'ragg' is still unavailable. ",
      "PNG export was skipped to avoid macOS X11/Cairo errors. ",
      "The PDF file has already been exported and should be used as the ",
      "primary SCI vector figure."
    )
    
    return(invisible(FALSE))
  }
  
  ok <- FALSE
  
  tryCatch(
    {
      ragg::agg_png(
        filename = filename,
        width = width,
        height = height,
        units = "in",
        res = res,
        background = "white"
      )
      
      ok <- TRUE
      
      draw_fun()
    },
    error = function(e) {
      
      warning(
        "PNG export failed: ",
        conditionMessage(e),
        "\nThe PDF export is unaffected."
      )
    },
    finally = {
      
      if (ok && grDevices::dev.cur() != 1) {
        grDevices::dev.off()
      }
      
      try(
        circos.clear(),
        silent = TRUE
      )
      
      invisible(gc())
    }
  )
  
  invisible(TRUE)
}

############################################################
############################################################

if (!file.exists(CORE_FILE)) {
  stop(
    "找不到最终严格 core 文件：\n",
    CORE_FILE,
    "\n请确认已运行严格版本 coregene 脚本。"
  )
}

if (!file.exists(STRICT_UPSET_FILE)) {
  stop(
    "找不到严格版 UpSet RDS：\n",
    STRICT_UPSET_FILE,
    "\n请先将严格版本 coregene 脚本中 strict_upset_groups.rds ",
    "统一导出到 FINAL_DIR，再重新运行本脚本。"
  )
}

############################################################
## 2. Helper functions
############################################################

zscore_cap <- function(x, cap = 2) {
  x <- as.numeric(x)
  z <- as.numeric(scale(x))
  z[is.na(z)] <- 0
  z <- pmax(pmin(z, cap), -cap)
  z
}

minmax01 <- function(x) {
  x <- as.numeric(x)
  
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  
  x[is.na(x)] <- median(x, na.rm = TRUE)
  
  rx <- range(x, na.rm = TRUE)
  
  if (!all(is.finite(rx)) || diff(rx) == 0) {
    return(rep(0.5, length(x)))
  }
  
  (x - rx[1]) / diff(rx)
}

get_hclust_order <- function(sub_mat) {
  if (is.null(sub_mat) || nrow(sub_mat) <= 1) {
    return(rownames(sub_mat))
  }
  
  hc <- hclust(
    dist(sub_mat),
    method = "complete"
  )
  
  rownames(sub_mat)[hc$order]
}

############################################################
############################################################

core_raw <- fread(
  CORE_FILE,
  data.table = TRUE
)

required_core_cols <- c(
  "gene",
  "core_direction",
  "bulk_logFC",
  "bulk_adj.P.Val",
  "sc_logFC",
  "sc_adj.P.Val",
  "bulk_kME",
  "bulk_GS_AML",
  "sc_kME_final"
)

missing_core_cols <- setdiff(
  required_core_cols,
  colnames(core_raw)
)

if (length(missing_core_cols) > 0) {
  stop(
    "02_Final_core_signature_v2.tsv 缺少必要列：",
    paste(missing_core_cols, collapse = ", ")
  )
}

core_dt <- core_raw[, .(
  gene = as.character(gene),
  
  direction = fifelse(
    core_direction == "core_up",
    "up",
    fifelse(
      core_direction == "core_down",
      "down",
      NA_character_
    )
  ),
  
  bulk_logFC = as.numeric(bulk_logFC),
  bulk_adjP = as.numeric(`bulk_adj.P.Val`),
  
  sc_logFC = as.numeric(sc_logFC),
  sc_adjP = as.numeric(`sc_adj.P.Val`),
  
  bulk_kME = as.numeric(bulk_kME),
  bulk_GS_AML = as.numeric(bulk_GS_AML),
  sc_kME_final = as.numeric(sc_kME_final)
)]

core_dt <- core_dt[
  !is.na(gene) &
    gene != "" &
    direction %in% c("up", "down")
]

if (anyDuplicated(core_dt$gene) > 0) {
  dup_gene <- unique(
    core_dt$gene[duplicated(core_dt$gene)]
  )
  
  stop(
    "最终 core 文件存在重复基因：",
    paste(dup_gene, collapse = ", ")
  )
}

if (nrow(core_dt) != 62) {
  stop(
    "最终严格 core 必须为 62 genes；当前读取到 ",
    nrow(core_dt),
    " genes。停止绘图，避免误用旧版数据。"
  )
}

if (anyNA(core_dt)) {
  warning(
    "core_dt 中存在 NA；请重点检查最终 core 文件。"
  )
}

bad_direction <- core_dt[
  (direction == "up" & (bulk_logFC <= 0 | sc_logFC <= 0)) |
    (direction == "down" & (bulk_logFC >= 0 | sc_logFC >= 0))
]

if (nrow(bad_direction) > 0) {
  stop(
    "发现 core_direction 与 bulk/scRNA logFC 方向不一致的基因：",
    paste(bad_direction$gene, collapse = ", ")
  )
}

cat("\nFinal strict core loaded:\n")
cat("Total:", nrow(core_dt), "\n")
print(table(core_dt$direction))

############################################################
##
## 4.1 Directional evidence:
##     Bulk logFC
##     scRNA logFC
##     signed -log10(FDR)
##
## 4.2 Network-strength evidence:
##     |Bulk kME|
##     |Bulk GS|
##     |sc kME final|
############################################################

core_dt[, bulk_signed_sig :=
          sign(bulk_logFC) *
          (-log10(pmax(bulk_adjP, 1e-300)))]

core_dt[, sc_signed_sig :=
          sign(sc_logFC) *
          (-log10(pmax(sc_adjP, 1e-300)))]

core_dt[, bulk_kME_abs := abs(bulk_kME)]
core_dt[, bulk_GS_abs := abs(bulk_GS_AML)]
core_dt[, sc_kME_abs := abs(sc_kME_final)]

############################################################
##
############################################################

dir_cols <- c(
  "sc_signed_sig",
  "bulk_signed_sig",
  "sc_logFC",
  "bulk_logFC"
)

dir_labels <- c(
  "sc signed -log10(FDR)",
  "Bulk signed -log10(FDR)",
  "sc logFC",
  "Bulk logFC"
)

net_cols <- c(
  "sc_kME_abs",
  "bulk_GS_abs",
  "bulk_kME_abs"
)

net_labels <- c(
  "sc |kME|",
  "Bulk |GS|",
  "Bulk |kME|"
)

mat_dir_raw <- as.matrix(
  core_dt[, ..dir_cols]
)

rownames(mat_dir_raw) <- core_dt$gene
colnames(mat_dir_raw) <- dir_labels

mat_net_raw <- as.matrix(
  core_dt[, ..net_cols]
)

rownames(mat_net_raw) <- core_dt$gene
colnames(mat_net_raw) <- net_labels

## Directional evidence:
mat_dir_z <- apply(
  mat_dir_raw,
  2,
  zscore_cap,
  cap = 2
)

mat_dir_z <- as.matrix(mat_dir_z)
rownames(mat_dir_z) <- rownames(mat_dir_raw)
colnames(mat_dir_z) <- colnames(mat_dir_raw)

## Network strength:
mat_net_01 <- apply(
  mat_net_raw,
  2,
  minmax01
)

mat_net_01 <- as.matrix(mat_net_01)
rownames(mat_net_01) <- rownames(mat_net_raw)
colnames(mat_net_01) <- colnames(mat_net_raw)

############################################################
##
############################################################

mat_net_cluster <- apply(
  mat_net_raw,
  2,
  zscore_cap,
  cap = 2
)

mat_net_cluster <- as.matrix(mat_net_cluster)
rownames(mat_net_cluster) <- rownames(mat_net_raw)
colnames(mat_net_cluster) <- colnames(mat_net_raw)

mat_cluster <- cbind(
  mat_dir_z,
  mat_net_cluster
)

core_dt[, direction := factor(
  direction,
  levels = c("down", "up")
)]

down_genes <- core_dt[
  direction == "down",
  gene
]

up_genes <- core_dt[
  direction == "up",
  gene
]

ord_down <- get_hclust_order(
  mat_cluster[
    down_genes,
    ,
    drop = FALSE
  ]
)

ord_up <- get_hclust_order(
  mat_cluster[
    up_genes,
    ,
    drop = FALSE
  ]
)

gene_order <- c(
  ord_down,
  ord_up
)

core_dt <- core_dt[
  match(gene_order, gene)
]

mat_dir_z <- mat_dir_z[
  gene_order,
  ,
  drop = FALSE
]

mat_net_01 <- mat_net_01[
  gene_order,
  ,
  drop = FALSE
]

core_dt$direction <- factor(
  as.character(core_dt$direction),
  levels = c("down", "up")
)

cat("\nSector size for plotting:\n")
cat(
  "AML-down:",
  sum(core_dt$direction == "down"),
  "\n"
)
cat(
  "AML-up:",
  sum(core_dt$direction == "up"),
  "\n"
)

############################################################
############################################################

gene_groups_raw <- readRDS(
  STRICT_UPSET_FILE
)

required_sets_raw <- c(
  "Bulk_DEG",
  "Bulk_AML_WGCNA",
  "scRNA_DEG",
  "sc_hdWGCNA",
  "Final_core"
)

missing_sets <- setdiff(
  required_sets_raw,
  names(gene_groups_raw)
)

if (length(missing_sets) > 0) {
  stop(
    "strict_upset_groups.rds 缺少集合：",
    paste(missing_sets, collapse = ", ")
  )
}

gene_groups <- list(
  Bulk_DEG = unique(
    as.character(gene_groups_raw$Bulk_DEG)
  ),
  
  Bulk_WGCNA = unique(
    as.character(gene_groups_raw$Bulk_AML_WGCNA)
  ),
  
  scRNA_DEG = unique(
    as.character(gene_groups_raw$scRNA_DEG)
  ),
  
  sc_hdWGCNA = unique(
    as.character(gene_groups_raw$sc_hdWGCNA)
  ),
  
  Final_core = unique(
    as.character(gene_groups_raw$Final_core)
  )
)

if (!setequal(
  gene_groups$Final_core,
  core_dt$gene
)) {
  stop(
    "strict_upset_groups.rds 中 Final_core 与 ",
    "02_Final_core_signature_v2.tsv 不完全一致。\n",
    "停止绘图，避免混用不同版本结果。"
  )
}

if (length(gene_groups$Final_core) != 62) {
  stop(
    "UpSet Final_core 不是 62 genes。当前为：",
    length(gene_groups$Final_core)
  )
}

cat("\nStrict UpSet groups loaded:\n")
print(
  data.table(
    set = names(gene_groups),
    n_genes = sapply(
      gene_groups,
      length
    )
  )
)

############################################################
## 8. SCI-style center UpSet
##
############################################################

make_sci_upset_plot <- function(
    gene_groups,
    set_order = c(
      "Final_core",
      "sc_hdWGCNA",
      "scRNA_DEG",
      "Bulk_WGCNA",
      "Bulk_DEG"
    ),
    set_labels = c(
      "Final_core" = "Final core",
      "sc_hdWGCNA" = "sc hdWGCNA",
      "scRNA_DEG" = "scRNA DEG",
      "Bulk_WGCNA" = "Bulk WGCNA",
      "Bulk_DEG" = "Bulk DEG"
    ),
    remove_singletons = TRUE,
    top_n = 12
) {
  
  gene_groups <- lapply(
    gene_groups,
    function(x) unique(as.character(x))
  )
  
  gene_groups <- gene_groups[
    set_order
  ]
  
  all_genes <- sort(
    unique(
      unlist(
        gene_groups,
        use.names = FALSE
      )
    )
  )
  
  mem <- sapply(
    set_order,
    function(ss) {
      all_genes %in% gene_groups[[ss]]
    }
  )
  
  mem <- as.data.table(mem)
  
  mem[, gene := all_genes]
  
  mem[, n_active := rowSums(
    .SD
  ),
  .SDcols = set_order]
  
  mem <- mem[
    n_active > 0
  ]
  
  if (remove_singletons) {
    mem <- mem[
      n_active >= 2
    ]
  }
  
  comb_dt <- mem[
    ,
    .N,
    by = set_order
  ]
  
  setnames(
    comb_dt,
    "N",
    "size"
  )
  
  comb_dt[, n_active := rowSums(
    .SD
  ),
  .SDcols = set_order]
  
  comb_dt[, final_core_hit :=
            as.logical(Final_core)]
  
  comb_dt <- comb_dt[
    order(
      -size,
      -n_active
    )
  ]
  
  core_comb <- comb_dt[
    final_core_hit == TRUE
  ]
  
  other_comb <- comb_dt[
    final_core_hit == FALSE
  ]
  
  if (!is.null(top_n)) {
    
    n_other_keep <- max(
      0,
      top_n - nrow(core_comb)
    )
    
    other_comb <- head(
      other_comb,
      n_other_keep
    )
    
    comb_dt <- rbind(
      core_comb,
      other_comb,
      use.names = TRUE
    )
    
    comb_dt <- comb_dt[
      order(
        -size,
        -n_active
      )
    ]
  }
  
  comb_dt[, x := seq_len(.N)]
  
  ##########################################################
  ## dot matrix data
  ##########################################################
  
  mat_dt <- melt(
    comb_dt,
    id.vars = c(
      "x",
      "size",
      "n_active",
      "final_core_hit"
    ),
    measure.vars = set_order,
    variable.name = "set",
    value.name = "active"
  )
  
  set_bottom_to_top <- rev(
    set_order
  )
  
  mat_dt[, y := match(
    set,
    set_bottom_to_top
  )]
  
  line_dt <- mat_dt[
    active == TRUE,
    .(
      ymin = min(y),
      ymax = max(y),
      final_core_hit = final_core_hit[1]
    ),
    by = x
  ]
  
  line_dt <- line_dt[
    ymax > ymin
  ]
  
  bg_dt <- data.table(
    y = seq_along(
      set_bottom_to_top
    ),
    set = set_bottom_to_top
  )
  
  ##########################################################
  ## set size data
  ##########################################################
  
  set_size_dt <- data.table(
    set = set_order,
    n = sapply(
      gene_groups[set_order],
      length
    )
  )
  
  set_size_dt[, y := match(
    set,
    set_bottom_to_top
  )]
  
  x_max <- nrow(
    comb_dt
  )
  
  col_intersection <- "#4F8FB3"
  col_core <- "#B2182B"
  col_set_size <- "#C98E96"
  col_matrix <- "#222222"
  col_inactive <- "#D2D2D2"
  
  ##########################################################
  ## upper intersection bars
  ##########################################################
  
  p_top <- ggplot(
    comb_dt,
    aes(
      x = x,
      y = size
    )
  ) +
    geom_col(
      aes(
        fill = final_core_hit
      ),
      width = 0.64
    ) +
    geom_text(
      aes(
        label = size
      ),
      vjust = -0.30,
      size = 2.15,
      family = FONT_FAMILY,
      color = "#222222"
    ) +
    scale_fill_manual(
      values = c(
        `FALSE` = col_intersection,
        `TRUE` = col_core
      ),
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(
        0.5,
        x_max + 0.5
      ),
      expand = c(
        0,
        0
      )
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0,
          0.15
        )
      )
    ) +
    labs(
      y = "Intersection\nsize",
      x = NULL
    ) +
    theme_minimal(
      base_family = FONT_FAMILY,
      base_size = 7.4
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_text(
        face = "bold",
        size = 7.0,
        color = "#222222"
      ),
      axis.text.y = element_text(
        size = 6.2,
        color = "#222222"
      ),
      axis.ticks.x = element_blank(),
      plot.margin = margin(
        0,
        0,
        0,
        0
      ),
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.background = element_rect(
        fill = "transparent",
        color = NA
      )
    )
  
  ##########################################################
  ## left set-size bars
  ##########################################################
  
  left_max <- max(
    set_size_dt$n,
    na.rm = TRUE
  )
  
  left_breaks <- pretty(
    c(
      0,
      left_max
    ),
    n = 3
  )
  
  left_breaks <- left_breaks[
    left_breaks >= 0 &
      left_breaks <= left_max
  ]
  
  p_left <- ggplot(
    set_size_dt
  ) +
    geom_rect(
      aes(
        xmin = 0,
        xmax = n,
        ymin = y - 0.29,
        ymax = y + 0.29
      ),
      fill = col_set_size,
      color = "white",
      linewidth = 0.10
    ) +
    scale_y_continuous(
      limits = c(
        0.5,
        length(
          set_bottom_to_top
        ) + 0.5
      ),
      expand = c(
        0,
        0
      )
    ) +
    scale_x_reverse(
      limits = c(
        left_max * 1.06,
        0
      ),
      breaks = left_breaks,
      expand = c(
        0,
        0
      )
    ) +
    labs(
      x = "Set size",
      y = NULL
    ) +
    theme_minimal(
      base_family = FONT_FAMILY,
      base_size = 7.4
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_text(
        face = "bold",
        size = 6.8,
        color = "#222222"
      ),
      axis.text.x = element_text(
        size = 5.8,
        color = "#222222"
      ),
      plot.margin = margin(
        0,
        1,
        0,
        0
      ),
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.background = element_rect(
        fill = "transparent",
        color = NA
      )
    )
  
  ##########################################################
  ## center matrix
  ##########################################################
  
  p_mat <- ggplot() +
    geom_rect(
      data = bg_dt[
        y %% 2 == 0
      ],
      aes(
        xmin = 0.5,
        xmax = x_max + 0.5,
        ymin = y - 0.5,
        ymax = y + 0.5
      ),
      fill = "#F4F4F4",
      color = NA
    ) +
    geom_point(
      data = mat_dt,
      aes(
        x = x,
        y = y
      ),
      color = col_inactive,
      size = 1.25
    ) +
    geom_segment(
      data = line_dt[
        final_core_hit == FALSE
      ],
      aes(
        x = x,
        xend = x,
        y = ymin,
        yend = ymax
      ),
      color = col_matrix,
      linewidth = 0.38
    ) +
    geom_segment(
      data = line_dt[
        final_core_hit == TRUE
      ],
      aes(
        x = x,
        xend = x,
        y = ymin,
        yend = ymax
      ),
      color = col_core,
      linewidth = 0.48
    ) +
    geom_point(
      data = mat_dt[
        active == TRUE &
          final_core_hit == FALSE
      ],
      aes(
        x = x,
        y = y
      ),
      color = col_matrix,
      size = 1.48
    ) +
    geom_point(
      data = mat_dt[
        active == TRUE &
          final_core_hit == TRUE
      ],
      aes(
        x = x,
        y = y
      ),
      color = col_core,
      size = 1.62
    ) +
    scale_x_continuous(
      limits = c(
        0.5,
        x_max + 0.5
      ),
      expand = c(
        0,
        0
      )
    ) +
    scale_y_continuous(
      breaks = seq_along(
        set_bottom_to_top
      ),
      labels = set_labels[
        set_bottom_to_top
      ],
      limits = c(
        0.5,
        length(
          set_bottom_to_top
        ) + 0.5
      ),
      expand = c(
        0,
        0
      )
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal(
      base_family = FONT_FAMILY,
      base_size = 7.4
    ) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(
        size = 6.7,
        color = "#222222"
      ),
      plot.margin = margin(
        0,
        0,
        0,
        0
      ),
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      panel.background = element_rect(
        fill = "transparent",
        color = NA
      )
    )
  
  ##########################################################
  ## combine
  ##########################################################
  
  blank <- ggplot() +
    theme_void()
  
  p_all <- cowplot::plot_grid(
    blank,
    p_top,
    p_left,
    p_mat,
    ncol = 2,
    rel_widths = c(
      0.25,
      0.75
    ),
    rel_heights = c(
      0.64,
      0.36
    ),
    align = "hv",
    axis = "tblr"
  )
  
  return(
    p_all
  )
}

p_center_upset <- make_sci_upset_plot(
  gene_groups = gene_groups,
  remove_singletons = TRUE,
  top_n = 12
)

############################################################
############################################################

dir_col_fun <- colorRamp2(
  c(
    -2,
    0,
    2
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

net_col_fun <- colorRamp2(
  c(
    0,
    0.5,
    1
  ),
  c(
    "#F7F7F7",
    "#C2A5CF",
    "#762A83"
  )
)

dir_fill <- c(
  "down" = "#A6CEE3",
  "up" = "#FDAE61"
)

############################################################
## 10. track label helper
############################################################

add_heatmap_track_labels <- function(
    track_index,
    labels,
    sector_name = "up",
    line_mm = 1.7,
    text_mm = 2.3,
    cex = 1.03
) {
  
  circos.track(
    track.index = track_index,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      if (
        CELL_META$sector.index ==
        sector_name
      ) {
        
        n <- length(
          labels
        )
        
        cell_height <- (
          CELL_META$cell.ylim[2] -
            CELL_META$cell.ylim[1]
        ) / n
        
        y_coords <- seq(
          CELL_META$cell.ylim[1] +
            cell_height / 2,
          CELL_META$cell.ylim[2] -
            cell_height / 2,
          length.out = n
        )
        
        for (
          i in seq_len(n)
        ) {
          
          circos.lines(
            c(
              CELL_META$cell.xlim[2],
              CELL_META$cell.xlim[2] +
                convert_x(
                  line_mm,
                  "mm"
                )
            ),
            c(
              y_coords[i],
              y_coords[i]
            ),
            col = "grey35",
            lwd = 0.75
          )
        }
        
        circos.text(
          rep(
            CELL_META$cell.xlim[2],
            n
          ) +
            convert_x(
              text_mm,
              "mm"
            ),
          y_coords,
          rev(
            labels
          ),
          cex = cex,
          adj = c(
            0,
            0.5
          ),
          facing = "inside",
          font = 1,
          col = "#222222"
        )
      }
    }
  )
}

############################################################
############################################################

draw_main_figure <- function() {
  
  try(
    circos.clear(),
    silent = TRUE
  )
  
  circos.par(
    start.degree = 70,
    gap.after = c(
      8,
      24
    ),
    track.margin = c(
      0.002,
      0.004
    ),
    cell.padding = c(
      0,
      0,
      0,
      0
    ),
    points.overflow.warning = FALSE
  )
  
  ##########################################################
  ##########################################################
  
  circos.heatmap.initialize(
    mat_dir_z,
    split = core_dt$direction,
    cluster = FALSE
  )
  
  ##########################################################
  ##########################################################
  
  circos.heatmap(
    mat_dir_z,
    split = core_dt$direction,
    cluster = FALSE,
    bg.border = "grey35",
    bg.lwd = 0.75,
    cell.border = "white",
    cell.lwd = 0.26,
    rownames.side = "outside",
    rownames.cex = 0.98,
    rownames.font = 1,
    rownames.col = "#222222",
    col = dir_col_fun,
    track.height = 0.165
  )
  
  dir_track_index <- get.current.track.index()
  
  add_heatmap_track_labels(
    track_index = dir_track_index,
    labels = colnames(
      mat_dir_z
    ),
    sector_name = "up",
    line_mm = 1.7,
    text_mm = 2.3,
    cex = 1.00
  )
  
  ##########################################################
  ##########################################################
  
  circos.heatmap(
    mat_net_01,
    split = core_dt$direction,
    cluster = FALSE,
    bg.border = "grey35",
    bg.lwd = 0.75,
    cell.border = "white",
    cell.lwd = 0.26,
    rownames.side = NULL,
    col = net_col_fun,
    track.height = 0.120
  )
  
  net_track_index <- get.current.track.index()
  
  add_heatmap_track_labels(
    track_index = net_track_index,
    labels = colnames(
      mat_net_01
    ),
    sector_name = "up",
    line_mm = 1.7,
    text_mm = 2.3,
    cex = 1.00
  )
  
  ##########################################################
  ## 11.3 AML-up / AML-down sector labels
  ##########################################################
  
  circos.track(
    ylim = c(
      0,
      1
    ),
    track.height = 0.045,
    bg.col = adjustcolor(
      dir_fill[
        levels(
          core_dt$direction
        )
      ],
      alpha.f = 0.35
    ),
    bg.border = "grey35",
    panel.fun = function(x, y) {
      
      lab <- ifelse(
        CELL_META$sector.index == "up",
        "AML-up",
        "AML-down"
      )
      
      circos.text(
        CELL_META$xcenter,
        0.52,
        lab,
        facing = "bending.inside",
        niceFacing = TRUE,
        cex = 1.08,
        font = 2,
        col = "#222222"
      )
    }
  )
  
  ##########################################################
  ##########################################################
  
  legend_dir <- Legend(
    title = "Directional\nz score",
    col_fun = dir_col_fun,
    at = c(
      -2,
      -1,
      0,
      1,
      2
    ),
    title_gp = gpar(
      fontsize = 10.5,
      fontface = "bold"
    ),
    labels_gp = gpar(
      fontsize = 8.5
    ),
    grid_width = unit(
      3.6,
      "mm"
    ),
    grid_height = unit(
      3.2,
      "mm"
    )
  )
  
  legend_net <- Legend(
    title = "Relative network\nstrength",
    col_fun = net_col_fun,
    at = c(
      0,
      0.5,
      1
    ),
    labels = c(
      "0",
      "0.5",
      "1"
    ),
    title_gp = gpar(
      fontsize = 10.5,
      fontface = "bold"
    ),
    labels_gp = gpar(
      fontsize = 8.5
    ),
    grid_width = unit(
      3.6,
      "mm"
    ),
    grid_height = unit(
      3.2,
      "mm"
    )
  )
  
  packed_legend <- packLegend(
    legend_dir,
    legend_net,
    direction = "vertical",
    gap = unit(
      2.5,
      "mm"
    )
  )
  
  draw(
    packed_legend,
    x = unit(
      0.955,
      "npc"
    ),
    y = unit(
      0.925,
      "npc"
    ),
    just = c(
      "right",
      "top"
    )
  )
  
  ##########################################################
  ##
  ## width 0.40 -> 0.48
  ## height 0.21 -> 0.265
  ##########################################################
  
  pushViewport(
    viewport(
      x = 0.478,
      y = 0.520,
      width = 0.45,
      height = 0.265,
      clip = "on"
    )
  )
  
  grid.draw(
    ggplotGrob(
      p_center_upset
    )
  )
  
  popViewport()
  
  ##########################################################
  ## 11.6 Final core label
  ##########################################################
  
  grid.text(
    paste0(
      "Final core genes\nn = ",
      length(
        gene_groups$Final_core
      )
    ),
    x = unit(
      0.50,
      "npc"
    ),
    y = unit(
      0.335,
      "npc"
    ),
    gp = gpar(
      fontsize = 17.0,
      fontface = "bold",
      fontfamily = FONT_FAMILY,
      lineheight = 0.92,
      col = "#222222"
    )
  )
  
  ##########################################################
  ## 11.7 Optional title
  ##########################################################
  
  if (SHOW_TITLE) {
    
    grid.text(
      "Cross-platform evidence landscape of the 62-gene AML core signature",
      x = unit(
        0.50,
        "npc"
      ),
      y = unit(
        0.965,
        "npc"
      ),
      gp = gpar(
        fontsize = 20,
        fontface = "bold",
        fontfamily = FONT_FAMILY,
        col = "#222222"
      )
    )
  }
  
  try(
    circos.clear(),
    silent = TRUE
  )
}

############################################################
## 12. Export PDF
############################################################

pdf(
  OUT_PDF,
  width = 14.0,
  height = 15.0,
  useDingbats = FALSE,
  family = FONT_FAMILY
)

draw_main_figure()

dev.off()

try(
  circos.clear(),
  silent = TRUE
)

invisible(gc())

############################################################
## 13. Export 600-dpi PNG
##
## Use ragg::agg_png instead of png(type = "cairo").
## This avoids the macOS/XQuartz errors:
##   libSM.6.dylib missing
##   libXrender.1.dylib missing
##   failed to load cairo DLL
############################################################

export_png_safe(
  filename = OUT_PNG,
  width = 14.0,
  height = 15.0,
  res = 600,
  draw_fun = draw_main_figure
)

############################################################
############################################################

plot_qc <- data.table(
  item = c(
    "Final core genes",
    "AML-up genes",
    "AML-down genes",
    "Bulk DEG set",
    "Bulk WGCNA set",
    "scRNA DEG set",
    "sc hdWGCNA set"
  ),
  n = c(
    nrow(core_dt),
    sum(
      core_dt$direction == "up"
    ),
    sum(
      core_dt$direction == "down"
    ),
    length(
      gene_groups$Bulk_DEG
    ),
    length(
      gene_groups$Bulk_WGCNA
    ),
    length(
      gene_groups$scRNA_DEG
    ),
    length(
      gene_groups$sc_hdWGCNA
    )
  )
)

fwrite(
  plot_qc,
  file.path(
    OUT_DIR,
    "Fig_main_circular_heatmap_plot_QC.tsv"
  ),
  sep = "\t"
)

cat("\n============================================\n")
cat("SCI final circular evidence figure completed.\n")
cat("PDF:\n", OUT_PDF, "\n")
cat("PNG 600 dpi:\n", OUT_PNG, "\n")
cat("QC:\n")
print(plot_qc)
cat("============================================\n")
