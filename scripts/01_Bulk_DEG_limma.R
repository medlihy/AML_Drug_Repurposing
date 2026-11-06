############################################################
## GSE114868 Affymetrix HTA 2.0 microarray DEG analysis
## Comparison: AML patients vs Healthy donors
## Method: GEOquery + limma
## Author: revised for high-quality DEG workflow
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
options(stringsAsFactors = FALSE)
set.seed(20260501)
# Paths are supplied by config/paths.R (or paths_example.R).
## =========================
## 0. User parameters
## =========================
GSE_ID <- "GSE114868"
OUTDIR <- BULK_DEG_RESULTS_DIR
FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5  # Absolute log2 fold-change threshold.
USE_ARRAY_WEIGHTS <- TRUE
USE_ROBUST_EBAYES <- TRUE
USE_TREND_EBAYES <- TRUE
KEEP_MAIN_CATEGORY <- TRUE
KEEP_LOCUS_TYPES <- NULL
TOP_N_HEATMAP <- 50

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

## =========================
## 1. Packages
## =========================

library("Biobase")
library("limma")
library("dplyr")
library("tibble")
library("ggplot2")
library("ggrepel")
library("pheatmap")
library("statmod")
library("data.table")

## =========================
## 2. Read local GEO files
## =========================

LOCAL_DIR <- GSE114868_DATA_DIR

SERIES_MATRIX_FILE <- file.path(LOCAL_DIR, "GSE114868_series_matrix.txt.gz")
GPL_FILE <- file.path(LOCAL_DIR, "GPL17586-45144.txt")

stopifnot(file.exists(SERIES_MATRIX_FILE))
stopifnot(file.exists(GPL_FILE))

eset <- GEOquery::getGEO(
  filename = SERIES_MATRIX_FILE,
  GSEMatrix = TRUE,
  getGPL = FALSE
)

if (is.list(eset)) {
  eset <- eset[[1]]
}

expr_raw <- Biobase::exprs(eset)
pd <- Biobase::pData(eset)

fd <- data.table::fread(
  GPL_FILE,
  data.table = FALSE,
  check.names = FALSE
)

if (!"ID" %in% colnames(fd)) {
  stop("Annotation file does not contain column named 'ID'. Please check GPL17586 annotation file.")
}

rownames(fd) <- fd$ID

pd <- pd[colnames(expr_raw), , drop = FALSE]
stopifnot(identical(rownames(pd), colnames(expr_raw)))

missing_probe <- setdiff(rownames(expr_raw), rownames(fd))
if (length(missing_probe) > 0) {
  stop("Some probes are missing from annotation: ", paste(head(missing_probe, 10), collapse = ", "))
}

fd <- fd[rownames(expr_raw), , drop = FALSE]
stopifnot(identical(rownames(fd), rownames(expr_raw)))

message("Expression matrix: ", nrow(expr_raw), " probes/transcript clusters x ", ncol(expr_raw), " samples")
message("Sample metadata: ", nrow(pd), " samples")
message("GPL annotation: ", nrow(fd), " probes/transcript clusters")
## =========================
## 3. Build phenotype table: AML vs Healthy
## =========================
pd$group <- NA_character_
pd$group[grepl("AML patient", pd$title, ignore.case = TRUE)] <- "AML"
pd$group[grepl("Healthy donor", pd$title, ignore.case = TRUE)] <- "Healthy"

if (any(is.na(pd$group))) {
  bad <- rownames(pd)[is.na(pd$group)]
  stop("Some samples could not be assigned to AML/Healthy by title: ", paste(bad, collapse = ", "))
}

pd$group <- factor(pd$group, levels = c("Healthy", "AML"))
print(table(pd$group))
write.table(pd, file = file.path(OUTDIR, "01_sample_metadata.tsv"), sep = "\t", quote = FALSE, col.names = NA)

## For GSE114868, expected: AML = 194, Healthy = 20.
if (!all(table(pd$group)[c("Healthy", "AML")] == c(20, 194))) {
  warning("Sample numbers are not 20 Healthy and 194 AML. Check pd$title and group parsing manually.")
}

## =========================
## 4. Expression preprocessing
## =========================
expr <- expr_raw
mode(expr) <- "numeric"

## 4.1 Remove probes with too many missing values
keep_na <- rowMeans(is.na(expr)) <= 0.20
expr <- expr[keep_na, , drop = FALSE]
fd <- fd[rownames(expr), , drop = FALSE]

## 4.2 Log2 transform only when needed
qx <- as.numeric(stats::quantile(expr, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
need_log2 <- (qx[5] > 100) || ((qx[6] - qx[1] > 50) && qx[2] > 0)
if (need_log2) {
  expr[expr <= 0] <- NA
  expr <- log2(expr)
}

## 4.3 Median impute a small number of remaining NA values, if present
if (anyNA(expr)) {
  message("Imputing remaining NA values by row median. Rows with >20% NA have already been removed.")
  row_median <- apply(expr, 1, stats::median, na.rm = TRUE)
  idx <- which(is.na(expr), arr.ind = TRUE)
  expr[idx] <- row_median[idx[, 1]]
}

## 4.4 QC before normalization
pdf(file.path(OUTDIR, "02_boxplot_before_normalization.pdf"), width = 12, height = 5)
boxplot(expr, outline = FALSE, las = 2, cex.axis = 0.45, main = "Before normalization", ylab = "Expression")
dev.off()

pdf(file.path(OUTDIR, "03_density_before_normalization.pdf"), width = 7, height = 5)
plotDensities(expr, group = pd$group, main = "Density before normalization")
dev.off()

## 4.5 Between-array normalization
## For GEO processed matrices this may be redundant, but it is kept for a consistent single-channel array workflow.
expr_norm <- limma::normalizeBetweenArrays(expr, method = "quantile")

pdf(file.path(OUTDIR, "04_boxplot_after_normalization.pdf"), width = 12, height = 5)
boxplot(expr_norm, outline = FALSE, las = 2, cex.axis = 0.45, main = "After quantile normalization", ylab = "log2 expression")
dev.off()

pdf(file.path(OUTDIR, "05_density_after_normalization.pdf"), width = 7, height = 5)
plotDensities(expr_norm, group = pd$group, main = "Density after normalization")
dev.off()

## =========================
## 5. Probe annotation and gene-level summarization
## =========================
find_col <- function(dat, choices) {
  nm <- names(dat)
  hit <- nm[tolower(nm) %in% tolower(choices)]
  if (length(hit) > 0) return(hit[1])
  hit2 <- nm[grepl(paste(choices, collapse = "|"), nm, ignore.case = TRUE)]
  if (length(hit2) > 0) return(hit2[1])
  return(NA_character_)
}

parse_affy_gene_assignment <- function(x, field_index = 2) {
  vapply(x, function(z) {
    if (is.na(z) || z == "" || z == "---") return(NA_character_)
    entries <- unlist(strsplit(z, " /// ", fixed = TRUE))
    vals <- vapply(entries, function(e) {
      fields <- unlist(strsplit(e, " // ", fixed = TRUE))
      if (length(fields) >= field_index) trimws(fields[field_index]) else NA_character_
    }, character(1))
    vals <- unique(vals[!is.na(vals) & vals != "" & vals != "---" & vals != "NULL"])
    if (length(vals) == 0) NA_character_ else paste(vals, collapse = "///")
  }, character(1))
}

gene_assign_col <- find_col(fd, c("gene_assignment", "Gene Assignment"))
category_col <- find_col(fd, c("category"))
locus_col <- find_col(fd, c("locus type", "locus_type"))

if (is.na(gene_assign_col)) {
  stop("Cannot find gene_assignment column in GPL annotation. Check fData(eset) column names.")
}

ann <- data.frame(
  probe_id = rownames(fd),
  gene_assignment = fd[[gene_assign_col]],
  category = if (!is.na(category_col)) fd[[category_col]] else NA_character_,
  locus_type = if (!is.na(locus_col)) fd[[locus_col]] else NA_character_,
  stringsAsFactors = FALSE
)

ann$gene_symbol_all <- parse_affy_gene_assignment(ann$gene_assignment, field_index = 2)
ann$entrez_id_all <- parse_affy_gene_assignment(ann$gene_assignment, field_index = 5)
ann$gene_symbol <- sub("///.*$", "", ann$gene_symbol_all)
ann$entrez_id <- sub("///.*$", "", ann$entrez_id_all)

## Optional HTA filtering
keep_ann <- !is.na(ann$gene_symbol) & ann$gene_symbol != ""
if (KEEP_MAIN_CATEGORY && !all(is.na(ann$category))) {
  keep_ann <- keep_ann & ann$category == "main"
}
if (!is.null(KEEP_LOCUS_TYPES) && !all(is.na(ann$locus_type))) {
  keep_ann <- keep_ann & ann$locus_type %in% KEEP_LOCUS_TYPES
}

expr_annot <- expr_norm[keep_ann, , drop = FALSE]
ann_annot <- ann[keep_ann, , drop = FALSE]

## Collapse multiple transcript clusters/probes to one gene symbol.
## Strategy: retain the probe/transcript cluster with the highest mean expression for each gene.
ann_annot$mean_expr <- rowMeans(expr_annot, na.rm = TRUE)
ord <- order(ann_annot$gene_symbol, -ann_annot$mean_expr)
ann_annot <- ann_annot[ord, , drop = FALSE]
expr_annot <- expr_annot[ord, , drop = FALSE]
keep_best <- !duplicated(ann_annot$gene_symbol)
ann_gene <- ann_annot[keep_best, , drop = FALSE]
expr_gene <- expr_annot[keep_best, , drop = FALSE]
rownames(expr_gene) <- ann_gene$gene_symbol

## Independent low-expression filter: remove bottom 20% by mean expression
mean_cut <- stats::quantile(rowMeans(expr_gene), probs = 0.20, na.rm = TRUE)
keep_mean <- rowMeans(expr_gene) > mean_cut
expr_gene <- expr_gene[keep_mean, , drop = FALSE]
ann_gene <- ann_gene[keep_mean, , drop = FALSE]

write.table(ann_gene, file = file.path(OUTDIR, "06_gene_annotation_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expr_gene, file = file.path(OUTDIR, "07_expression_gene_symbol_log2_norm.tsv"), sep = "\t", quote = FALSE, col.names = NA)

## =========================
## 6. QC after gene summarization
## =========================
## PCA
pca <- prcomp(t(expr_gene), scale. = TRUE)
pca_df <- data.frame(
  sample = rownames(pd),
  group = pd$group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
pc_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  labs(x = paste0("PC1 (", pc_var[1], "%)"), y = paste0("PC2 (", pc_var[2], "%)"), color = NULL,
       title = "PCA: GSE114868 AML vs Healthy") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTDIR, "08_PCA_AML_vs_Healthy.pdf"), p_pca, width = 6.5, height = 5)

## Sample correlation heatmap; 214 samples can be dense, so sample names are hidden
pdf(file.path(OUTDIR, "09_sample_correlation_heatmap.pdf"), width = 8, height = 7)
pheatmap::pheatmap(cor(expr_gene, method = "pearson"),
                   show_rownames = FALSE, show_colnames = FALSE,
                   annotation_col = data.frame(group = pd$group, row.names = rownames(pd)),
                   main = "Sample correlation")
dev.off()

## =========================
## 7. limma DEG analysis
## =========================
design <- model.matrix(~ 0 + pd$group)
colnames(design) <- levels(pd$group)
rownames(design) <- rownames(pd)

contrast_matrix <- limma::makeContrasts(AML_vs_Healthy = AML - Healthy, levels = design)

if (USE_ARRAY_WEIGHTS) {
  aw <- limma::arrayWeights(expr_gene, design = design)
  write.table(data.frame(sample = names(aw), group = pd$group, array_weight = aw),
              file = file.path(OUTDIR, "10_array_weights.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  fit <- limma::lmFit(expr_gene, design = design, weights = aw)
} else {
  fit <- limma::lmFit(expr_gene, design = design)
}

fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2, trend = USE_TREND_EBAYES, robust = USE_ROBUST_EBAYES)

res <- limma::topTable(fit2, coef = "AML_vs_Healthy", number = Inf, adjust.method = "BH", sort.by = "P")
res$gene_symbol <- rownames(res)
res <- dplyr::left_join(res, ann_gene, by = "gene_symbol")
res <- res[, c("gene_symbol", "probe_id", "entrez_id", "gene_symbol_all", "entrez_id_all",
               "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
               "category", "locus_type", "gene_assignment")]

res$regulation <- "Not_significant"
res$regulation[res$adj.P.Val < FDR_CUTOFF & res$logFC >= LOGFC_CUTOFF] <- "Up_in_AML"
res$regulation[res$adj.P.Val < FDR_CUTOFF & res$logFC <= -LOGFC_CUTOFF] <- "Down_in_AML"

write.table(res, file = file.path(OUTDIR, "11_limma_all_results_AML_vs_Healthy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

deg <- res[res$regulation != "Not_significant", ]
deg_up <- res[res$regulation == "Up_in_AML", ]
deg_down <- res[res$regulation == "Down_in_AML", ]

write.table(deg, file = file.path(OUTDIR, sprintf("12_DEG_FDR%.2f_logFC%.2f.tsv", FDR_CUTOFF, LOGFC_CUTOFF)), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_up, file = file.path(OUTDIR, "13_DEG_up_in_AML.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_down, file = file.path(OUTDIR, "14_DEG_down_in_AML.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_tab <- data.frame(
  comparison = "AML_vs_Healthy",
  n_genes_tested = nrow(res),
  FDR_cutoff = FDR_CUTOFF,
  logFC_cutoff = LOGFC_CUTOFF,
  DEG_total = nrow(deg),
  Up_in_AML = nrow(deg_up),
  Down_in_AML = nrow(deg_down)
)
write.table(summary_tab, file = file.path(OUTDIR, "15_DEG_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_tab)

## =========================
## 8. Visualization: volcano, MA, DEG heatmap
## =========================
plot_df <- res
plot_df$neg_log10_FDR <- -log10(pmax(plot_df$adj.P.Val, .Machine$double.xmin))
plot_df$label <- NA_character_
label_genes <- plot_df %>%
  dplyr::filter(regulation != "Not_significant") %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(gene_symbol)
plot_df$label[plot_df$gene_symbol %in% label_genes] <- plot_df$gene_symbol[plot_df$gene_symbol %in% label_genes]

p_volcano <- ggplot(plot_df, aes(x = logFC, y = neg_log10_FDR, color = regulation)) +
  geom_point(size = 1.2, alpha = 0.65) +
  geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    aes(label = label),
    max.overlaps = 30,
    size = 3,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Down_in_AML" = "#377EB8", "Not_significant" = "grey75", "Up_in_AML" = "#E41A1C")) +
  labs(title = "GSE114868: AML vs Healthy", x = "log2FC (AML / Healthy)", y = "-log10(FDR)", color = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(OUTDIR, "16_volcano_AML_vs_Healthy.pdf"), p_volcano, width = 7, height = 6)

p_ma <- ggplot(plot_df, aes(x = AveExpr, y = logFC, color = regulation)) +
  geom_point(size = 1.0, alpha = 0.6) +
  geom_hline(yintercept = c(-LOGFC_CUTOFF, 0, LOGFC_CUTOFF), linetype = c("dashed", "solid", "dashed"), linewidth = 0.3) +
  scale_color_manual(values = c("Down_in_AML" = "#377EB8", "Not_significant" = "grey75", "Up_in_AML" = "#E41A1C")) +
  labs(title = "MA plot", x = "Average expression", y = "log2FC (AML / Healthy)", color = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(OUTDIR, "17_MAplot_AML_vs_Healthy.pdf"), p_ma, width = 6.5, height = 5)

if (nrow(deg) >= 2) {
  top_genes <- deg %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = min(TOP_N_HEATMAP, nrow(deg))) %>%
    dplyr::pull(gene_symbol)
  hm <- expr_gene[top_genes, , drop = FALSE]
  hm_z <- t(scale(t(hm)))
  hm_z[hm_z > 2] <- 2
  hm_z[hm_z < -2] <- -2
  anno_col <- data.frame(group = pd$group, row.names = rownames(pd))
  pdf(file.path(OUTDIR, "18_top_DEG_heatmap.pdf"), width = 8, height = 10)
  pheatmap::pheatmap(hm_z,
                     annotation_col = anno_col,
                     show_colnames = FALSE,
                     fontsize_row = 8,
                     main = paste0("Top ", length(top_genes), " DEGs"))
  dev.off()
}

## =========================
## 9. Save R objects and session info
## =========================
saveRDS(list(expr_gene = expr_gene, metadata = pd, annotation = ann_gene, fit = fit2, results = res),
        file = file.path(OUTDIR, "19_GSE114868_limma_objects.rds"))
writeLines(capture.output(sessionInfo()), con = file.path(OUTDIR, "20_sessionInfo.txt"))

message("Analysis finished. Results saved to: ", normalizePath(OUTDIR))
