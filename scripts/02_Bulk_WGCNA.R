############################################################
## GSE114868 high-quality WGCNA workflow
## Input: outputs from GSE114868_limma_DEG_high_quality.R
## Network: signed WGCNA using normalized gene-level expression
## Integration: limma DEG annotation + module-trait association + module-DEG enrichment + hub genes
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
set.seed(20260502)

## =========================
## 0. User parameters
## =========================

LIMMA_OUTDIR <- BULK_DEG_RESULTS_DIR
OUTDIR <- BULK_WGCNA_RESULTS_DIR

EXPR_FILE <- file.path(LIMMA_OUTDIR, "07_expression_gene_symbol_log2_norm.tsv")
META_FILE <- file.path(LIMMA_OUTDIR, "01_sample_metadata.tsv")
DEG_ALL_FILE <- file.path(LIMMA_OUTDIR, "11_limma_all_results_AML_vs_Healthy.tsv")
DEG_FILE <- file.path(LIMMA_OUTDIR, "12_DEG_FDR0.05_logFC0.50.tsv")

ARRAY_WEIGHT_FILE_CANDIDATES <- c(
  file.path(LIMMA_OUTDIR, "10_array_weights_with_sample_id.tsv"),
  file.path(LIMMA_OUTDIR, "10_array_weights.tsv")
)

GENE_FILTER_MODE <- "top_mad"       # options: "top_mad", "mad_quantile", "all"
TOP_MAD_N <- 12000
MAD_QUANTILE <- 0.50
MIN_MAD <- 0

USE_ARRAY_WEIGHT_FILTER <- TRUE
ARRAY_WEIGHT_CUTOFF <- 0.30
REMOVE_SAMPLE_OUTLIERS_BY_ZK <- FALSE
SAMPLE_ZK_CUTOFF <- -2.5

## WGCNA parameters
N_THREADS <- 8
NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
COR_TYPE <- "bicor"
MAX_P_OUTLIERS <- 0.05
SCALE_FREE_R2_CUTOFF <- 0.85
POWERS <- c(1:20, seq(22, 30, by = 2))
MIN_MODULE_SIZE <- 60
DEEP_SPLIT <- 2
MERGE_CUT_HEIGHT <- 0.25
MAX_BLOCK_SIZE <- 12000
SAVE_TOM <- FALSE

## Hub / DEG integration
DEG_FDR_CUTOFF <- 0.05
DEG_LOGFC_CUTOFF <- 1
HUB_KME_CUTOFF <- 0.80
HUB_GS_CUTOFF <- 0.20
TOP_HUB_PER_MODULE <- 30
KEY_MODULE_FDR_CUTOFF <- 0.05

## Cytoscape export for key modules
RUN_CYTOSCAPE_EXPORT <- TRUE
KEY_MODULE_N_FOR_CYTOSCAPE <- 5
CYTO_TOP_GENES_PER_MODULE <- 150
CYTO_EDGE_THRESHOLD <- 0.08
CYTO_MAX_EDGES <- 15000

RUN_MODULE_ORA <- TRUE
ORA_P_CUTOFF <- 0.05
ORA_Q_CUTOFF <- 0.20

TRAITS_TO_PLOT_IN_MAIN_HEATMAP <- c("AML_status")
MODULE_TRAIT_HEATMAP_ROW_HEIGHT <- 0.42
MODULE_TRAIT_HEATMAP_TEXT_CEX <- 0.72
SOFT_POWER_LABEL_POWERS <- c(1:10, 12, 16, 20, 24, 28, 30)

## =========================
## 1. Packages and folders
## =========================

required_pkgs <- c("WGCNA", "limma", "data.table", "ggplot2", "pheatmap", "matrixStats")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not installed: ", pkg, ". Please install it before running this script.")
  }
}

suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
})

try(WGCNA::allowWGCNAThreads(nThreads = N_THREADS), silent = TRUE)

DIRS <- list(
  input = file.path(OUTDIR, "00_input"),
  qc = file.path(OUTDIR, "01_QC"),
  power = file.path(OUTDIR, "02_soft_power"),
  module = file.path(OUTDIR, "03_modules"),
  trait = file.path(OUTDIR, "04_module_trait"),
  hub = file.path(OUTDIR, "05_hub_genes"),
  cyto = file.path(OUTDIR, "06_cytoscape"),
  enrich = file.path(OUTDIR, "07_enrichment"),
  object = file.path(OUTDIR, "08_R_objects")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

log_file <- file.path(OUTDIR, "00_WGCNA_run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}
if (file.exists(log_file)) file.remove(log_file)

## =========================
## 2. Helper functions
## =========================

read_matrix_first_col <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  x <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  if (ncol(x) < 2) stop("Expression matrix must have gene column + sample columns: ", file)
  rn <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE])
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  mat
}

read_metadata <- function(file, sample_names) {
  if (!file.exists(file)) stop("File not found: ", file)
  pd <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  if ("sample_id" %in% colnames(pd)) {
    rownames(pd) <- as.character(pd$sample_id)
  } else {
    first_col <- colnames(pd)[1]
    if (all(as.character(pd[[first_col]]) %in% sample_names)) {
      pd$sample_id <- as.character(pd[[first_col]])
      rownames(pd) <- pd$sample_id
      pd[[first_col]] <- NULL
    } else if (all(sample_names %in% rownames(pd))) {
      pd$sample_id <- rownames(pd)
    } else {
      stop("Cannot infer sample IDs from metadata. Make sure 01_sample_metadata.tsv contains row names or sample_id.")
    }
  }
  missing <- setdiff(sample_names, rownames(pd))
  if (length(missing) > 0) stop("Metadata missing samples: ", paste(head(missing, 10), collapse = ", "))
  pd <- pd[sample_names, , drop = FALSE]
  if (!"group" %in% colnames(pd)) {
    if ("title" %in% colnames(pd)) {
      pd$group <- NA_character_
      pd$group[grepl("AML patient", pd$title, ignore.case = TRUE)] <- "AML"
      pd$group[grepl("Healthy donor", pd$title, ignore.case = TRUE)] <- "Healthy"
    } else {
      stop("Metadata must contain group column, or title column that includes AML patient / Healthy donor.")
    }
  }
  pd$group <- factor(pd$group, levels = c("Healthy", "AML"))
  if (any(is.na(pd$group))) stop("Some samples have NA group after parsing metadata.")
  pd
}

collapse_duplicate_genes <- function(expr) {
  if (!anyDuplicated(rownames(expr))) return(expr)
  madv <- matrixStats::rowMads(expr, na.rm = TRUE)
  meanv <- rowMeans(expr, na.rm = TRUE)
  ord <- order(rownames(expr), -madv, -meanv)
  expr2 <- expr[ord, , drop = FALSE]
  expr2 <- expr2[!duplicated(rownames(expr2)), , drop = FALSE]
  expr2
}

read_deg_annotation <- function(all_file, deg_file, fdr_cut = 0.05, logfc_cut = 1) {
  if (file.exists(all_file)) {
    deg <- data.table::fread(all_file, data.table = FALSE, check.names = FALSE)
  } else if (file.exists(deg_file)) {
    warning("Full limma result file not found. Falling back to DEG-only file: ", deg_file)
    deg <- data.table::fread(deg_file, data.table = FALSE, check.names = FALSE)
  } else {
    warning("No DEG annotation file found. WGCNA will run without DEG integration.")
    return(data.frame(gene_symbol = character(0)))
  }
  if (!"gene_symbol" %in% colnames(deg)) colnames(deg)[1] <- "gene_symbol"
  deg$gene_symbol <- as.character(deg$gene_symbol)
  deg <- deg[!is.na(deg$gene_symbol) & deg$gene_symbol != "", , drop = FALSE]
  deg <- deg[!duplicated(deg$gene_symbol), , drop = FALSE]
  if (!"regulation" %in% colnames(deg)) {
    deg$regulation <- "Not_significant"
    if (all(c("adj.P.Val", "logFC") %in% colnames(deg))) {
      deg$regulation[deg$adj.P.Val < fdr_cut & deg$logFC >= logfc_cut] <- "Up_in_AML"
      deg$regulation[deg$adj.P.Val < fdr_cut & deg$logFC <= -logfc_cut] <- "Down_in_AML"
    }
  }
  keep_cols <- intersect(c("gene_symbol", "probe_id", "entrez_id", "gene_symbol_all", "entrez_id_all",
                           "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
                           "category", "locus_type", "regulation"), colnames(deg))
  deg[, keep_cols, drop = FALSE]
}

read_array_weights <- function(files, sample_names) {
  file <- files[file.exists(files)][1]
  if (is.na(file)) return(NULL)
  aw <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  weight_col <- grep("array_weight|weight", colnames(aw), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(weight_col)) {
    warning("Array weight file found but no array_weight column: ", file)
    return(NULL)
  }
  sample_col <- grep("sample_id|sample", colnames(aw), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(sample_col)) {
    if (nrow(aw) == length(sample_names)) {
      aw$sample_id <- sample_names
    } else {
      warning("Cannot infer sample column from array weight file: ", file)
      return(NULL)
    }
  } else {
    aw$sample_id <- as.character(aw[[sample_col]])
    if (!all(aw$sample_id %in% sample_names)) {
      suppressWarnings(idx <- as.integer(aw$sample_id))
      if (all(!is.na(idx)) && all(idx >= 1 & idx <= length(sample_names))) {
        aw$sample_id <- sample_names[idx]
      }
    }
  }
  aw$array_weight <- as.numeric(aw[[weight_col]])
  aw <- aw[aw$sample_id %in% sample_names, c("sample_id", "array_weight"), drop = FALSE]
  aw <- aw[!duplicated(aw$sample_id), , drop = FALSE]
  rownames(aw) <- aw$sample_id
  aw
}

p_adjust_matrix <- function(pmat) {
  out <- apply(pmat, 2, p.adjust, method = "BH")
  out <- as.data.frame(out, check.names = FALSE)
  rownames(out) <- rownames(pmat)
  as.matrix(out)
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 1e-300, "<1e-300", formatC(x, format = "e", digits = 1)))
}

plot_module_trait_heatmap_readable <- function(cor_mat, p_mat, fdr_mat, traits, file, title) {
  traits <- intersect(traits, colnames(cor_mat))
  if (length(traits) == 0) return(invisible(NULL))
  cm <- cor_mat[, traits, drop = FALSE]
  pm <- p_mat[, traits, drop = FALSE]
  fm <- fdr_mat[, traits, drop = FALSE]
  textMatrix <- matrix(
    paste0(sprintf("%.2f", as.vector(cm)), "\nP=", fmt_p(as.vector(pm)), "\nFDR=", fmt_p(as.vector(fm))),
    nrow = nrow(cm), ncol = ncol(cm), dimnames = dimnames(cm)
  )
  row_labs <- sub("^ME", "", rownames(cm))
  x_labs <- gsub("_", " ", traits)
  fig_h <- max(6.2, 2.2 + MODULE_TRAIT_HEATMAP_ROW_HEIGHT * nrow(cm))
  fig_w <- max(5.2, 3.5 + 2.2 * ncol(cm))
  pdf(file, width = fig_w, height = fig_h)
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  par(mar = c(6, 9, 4, 3))
  WGCNA::labeledHeatmap(Matrix = cm,
                        xLabels = x_labs,
                        yLabels = row_labs,
                        ySymbols = row_labs,
                        colorLabels = FALSE,
                        colors = WGCNA::blueWhiteRed(50),
                        textMatrix = textMatrix,
                        setStdMargins = FALSE,
                        cex.text = MODULE_TRAIT_HEATMAP_TEXT_CEX,
                        cex.lab = 0.85,
                        zlim = c(-1, 1),
                        main = title)
  invisible(NULL)
}

calc_sample_connectivity <- function(datExpr, group) {
  sample_adj <- WGCNA::adjacency(t(datExpr), type = "signed", power = 2,
                                 corFnc = COR_TYPE,
                                 corOptions = list(use = "p", maxPOutliers = MAX_P_OUTLIERS))
  k <- colSums(sample_adj, na.rm = TRUE) - 1
  zk <- as.numeric(scale(k))
  data.frame(sample_id = rownames(datExpr), group = group[rownames(datExpr)], sample_connectivity = k, Z.k = zk,
             stringsAsFactors = FALSE)
}

## =========================
## 3. Read limma outputs and build WGCNA input
## =========================

log_msg("Reading limma expression matrix: ", EXPR_FILE)
expr0 <- read_matrix_first_col(EXPR_FILE)
expr0 <- collapse_duplicate_genes(expr0)

log_msg("Reading sample metadata: ", META_FILE)
pd0 <- read_metadata(META_FILE, colnames(expr0))

log_msg("Original group counts: ", paste(capture.output(print(table(pd0$group))), collapse = " "))
if (!all(table(pd0$group)[c("Healthy", "AML")] == c(20, 194))) {
  warning("Before WGCNA filtering, group counts are not exactly Healthy=20 and AML=194. Check metadata.")
}

## DEG annotation from limma
log_msg("Reading limma DEG annotation.")
deg_anno <- read_deg_annotation(DEG_ALL_FILE, DEG_FILE, DEG_FDR_CUTOFF, DEG_LOGFC_CUTOFF)
write.table(deg_anno, file = file.path(DIRS$input, "01_limma_annotation_loaded.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Remove samples with extremely low limma array weights if available
array_weights <- read_array_weights(ARRAY_WEIGHT_FILE_CANDIDATES, colnames(expr0))
if (!is.null(array_weights)) {
  array_weights$group <- as.character(pd0[array_weights$sample_id, "group"])
  write.table(array_weights, file = file.path(DIRS$input, "02_array_weights_matched.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  low_weight_samples <- array_weights$sample_id[array_weights$array_weight < ARRAY_WEIGHT_CUTOFF]
} else {
  low_weight_samples <- character(0)
}

samples_to_remove <- character(0)
if (USE_ARRAY_WEIGHT_FILTER && length(low_weight_samples) > 0) {
  samples_to_remove <- union(samples_to_remove, low_weight_samples)
}

if (length(samples_to_remove) > 0) {
  removed_df <- data.frame(sample_id = samples_to_remove,
                           group = as.character(pd0[samples_to_remove, "group"]),
                           reason = "array_weight_below_cutoff",
                           stringsAsFactors = FALSE)
  if (!is.null(array_weights)) {
    removed_df$array_weight <- array_weights[removed_df$sample_id, "array_weight"]
  }
  write.table(removed_df, file = file.path(DIRS$input, "03_samples_removed_before_WGCNA.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  expr0 <- expr0[, setdiff(colnames(expr0), samples_to_remove), drop = FALSE]
  pd0 <- pd0[colnames(expr0), , drop = FALSE]
}

log_msg("Group counts after pre-WGCNA sample filtering: ", paste(capture.output(print(table(pd0$group))), collapse = " "))

## Missing value filtering and imputation
keep_na <- rowMeans(is.na(expr0)) <= 0.20
expr0 <- expr0[keep_na, , drop = FALSE]
if (anyNA(expr0)) {
  row_median <- apply(expr0, 1, stats::median, na.rm = TRUE)
  idx <- which(is.na(expr0), arr.ind = TRUE)
  expr0[idx] <- row_median[idx[, 1]]
}
expr0 <- expr0[apply(expr0, 1, function(z) all(is.finite(z))), , drop = FALSE]

## Gene filtering for WGCNA
madv <- matrixStats::rowMads(expr0, na.rm = TRUE)
keep_var <- is.finite(madv) & madv > MIN_MAD
expr_var <- expr0[keep_var, , drop = FALSE]
madv <- madv[keep_var]

if (GENE_FILTER_MODE == "top_mad") {
  n_keep <- min(TOP_MAD_N, nrow(expr_var))
  keep_idx <- order(madv, decreasing = TRUE)[seq_len(n_keep)]
  expr_wgcna <- expr_var[keep_idx, , drop = FALSE]
} else if (GENE_FILTER_MODE == "mad_quantile") {
  cutoff <- stats::quantile(madv, probs = MAD_QUANTILE, na.rm = TRUE)
  expr_wgcna <- expr_var[madv >= cutoff, , drop = FALSE]
} else if (GENE_FILTER_MODE == "all") {
  expr_wgcna <- expr_var
} else {
  stop("Unknown GENE_FILTER_MODE: ", GENE_FILTER_MODE)
}

expr_wgcna <- expr_wgcna[order(rownames(expr_wgcna)), , drop = FALSE]
write.table(expr_wgcna, file = file.path(DIRS$input, "04_expression_used_for_WGCNA_gene_by_sample.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(pd0, file = file.path(DIRS$input, "05_metadata_used_for_WGCNA.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

log_msg("WGCNA input before goodSamplesGenes: ", nrow(expr_wgcna), " genes x ", ncol(expr_wgcna), " samples")

datExpr0 <- t(expr_wgcna)

## WGCNA goodSamplesGenes check
gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    write.table(data.frame(gene_symbol = colnames(datExpr0)[!gsg$goodGenes]),
                file = file.path(DIRS$qc, "01_removed_bad_genes_goodSamplesGenes.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (sum(!gsg$goodSamples) > 0) {
    write.table(data.frame(sample_id = rownames(datExpr0)[!gsg$goodSamples]),
                file = file.path(DIRS$qc, "02_removed_bad_samples_goodSamplesGenes.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  pd0 <- pd0[rownames(datExpr0), , drop = FALSE]
}

## Sample QC: dendrogram, correlation heatmap, sample connectivity
sampleTree <- hclust(dist(datExpr0), method = "average")
datTraits0 <- data.frame(
  AML_status = ifelse(pd0$group == "AML", 1, 0),
  Healthy_status = ifelse(pd0$group == "Healthy", 1, 0),
  row.names = rownames(pd0)
)
traitColors <- WGCNA::numbers2colors(datTraits0, signed = FALSE)

pdf(file.path(DIRS$qc, "03_sample_dendrogram_with_traits_before_Zk_filter.pdf"), width = 12, height = 7)
WGCNA::plotDendroAndColors(sampleTree, traitColors,
                           groupLabels = colnames(datTraits0),
                           main = "Sample dendrogram and clinical traits")
dev.off()

sample_conn <- calc_sample_connectivity(datExpr0, setNames(as.character(pd0$group), rownames(pd0)))
write.table(sample_conn, file = file.path(DIRS$qc, "04_sample_connectivity_Zk.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_zk <- ggplot(sample_conn, aes(x = group, y = Z.k, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.18, size = 2, alpha = 0.85) +
  geom_hline(yintercept = SAMPLE_ZK_CUTOFF, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(title = "Sample connectivity Z.k", x = NULL, y = "Z.k")
ggsave(file.path(DIRS$qc, "05_sample_connectivity_Zk.pdf"), p_zk, width = 5.5, height = 4.5)

zk_outliers <- sample_conn$sample_id[sample_conn$Z.k < SAMPLE_ZK_CUTOFF]
if (REMOVE_SAMPLE_OUTLIERS_BY_ZK && length(zk_outliers) > 0) {
  write.table(sample_conn[sample_conn$sample_id %in% zk_outliers, ],
              file = file.path(DIRS$qc, "06_samples_removed_by_Zk.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  datExpr0 <- datExpr0[setdiff(rownames(datExpr0), zk_outliers), , drop = FALSE]
  pd0 <- pd0[rownames(datExpr0), , drop = FALSE]
  datTraits0 <- datTraits0[rownames(datExpr0), , drop = FALSE]
}

## Sample correlation heatmap on final WGCNA input
sample_cor <- WGCNA::bicor(t(datExpr0), use = "p", maxPOutliers = MAX_P_OUTLIERS)
sample_cor[is.na(sample_cor)] <- 0
diag(sample_cor) <- 1
anno_col <- data.frame(group = pd0$group, row.names = rownames(pd0))
pdf(file.path(DIRS$qc, "07_sample_correlation_heatmap_final.pdf"), width = 10, height = 8)
pheatmap::pheatmap(sample_cor,
                   annotation_col = anno_col,
                   annotation_row = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main = "Sample correlation: final WGCNA input")
dev.off()

log_msg("Final WGCNA input: ", ncol(datExpr0), " genes x ", nrow(datExpr0), " samples")
final_group_counts <- as.data.frame(table(pd0$group), stringsAsFactors = FALSE)
colnames(final_group_counts) <- c("group", "n_samples")
write.table(final_group_counts, file = file.path(DIRS$qc, "08_final_group_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## =========================
## 4. Soft-thresholding power selection
## =========================

log_msg("Picking soft-thresholding power.")
sft <- WGCNA::pickSoftThreshold(datExpr0,
                                powerVector = POWERS,
                                networkType = NETWORK_TYPE,
                                corFnc = COR_TYPE,
                                corOptions = list(use = "p", maxPOutliers = MAX_P_OUTLIERS),
                                verbose = 5)

sft_tab <- as.data.frame(sft$fitIndices)
sft_tab$signed_R2 <- -sign(sft_tab$slope) * sft_tab$SFT.R.sq
write.table(sft_tab, file = file.path(DIRS$power, "01_soft_threshold_fit_indices.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.na(sft$powerEstimate)) {
  softPower <- sft$powerEstimate
} else {
  candidates <- sft_tab$Power[sft_tab$signed_R2 >= SCALE_FREE_R2_CUTOFF]
  if (length(candidates) > 0) {
    softPower <- min(candidates)
  } else {
    softPower <- sft_tab$Power[which.max(sft_tab$signed_R2)]
    warning("No power reached signed R2 cutoff. Using power with maximum signed R2: ", softPower)
  }
}
writeLines(as.character(softPower), con = file.path(DIRS$power, "02_selected_soft_power.txt"))
log_msg("Selected softPower = ", softPower)

## Readable soft-power diagnostics: use points/lines and label only selected powers to avoid overplotting.
soft_power_label_idx <- sft_tab$Power %in% unique(c(SOFT_POWER_LABEL_POWERS, softPower))
pdf(file.path(DIRS$power, "03_soft_power_diagnostics.pdf"), width = 11, height = 5.5)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
plot(sft_tab$Power, sft_tab$signed_R2,
     xlab = "Soft threshold power", ylab = "Signed scale-free topology fit R^2",
     type = "b", pch = 19, cex = 0.75, main = "Scale independence")
text(sft_tab$Power[soft_power_label_idx], sft_tab$signed_R2[soft_power_label_idx],
     labels = sft_tab$Power[soft_power_label_idx], cex = 0.75, pos = 3, offset = 0.35)
abline(h = SCALE_FREE_R2_CUTOFF, col = "red", lty = 2)
abline(v = softPower, col = "blue", lty = 2)
legend("bottomright", legend = paste0("Selected power = ", softPower), bty = "n", cex = 0.85)

plot(sft_tab$Power, sft_tab$mean.k.,
     xlab = "Soft threshold power", ylab = "Mean connectivity",
     type = "b", pch = 19, cex = 0.75, main = "Mean connectivity")
text(sft_tab$Power[soft_power_label_idx], sft_tab$mean.k.[soft_power_label_idx],
     labels = sft_tab$Power[soft_power_label_idx], cex = 0.75, pos = 3, offset = 0.35)
abline(v = softPower, col = "blue", lty = 2)
legend("topright", legend = paste0("Selected power = ", softPower), bty = "n", cex = 0.85)
par(oldpar)
dev.off()

## Also save two single-panel versions for manuscripts/supplementary figures.
pdf(file.path(DIRS$power, "03a_scale_independence_readable.pdf"), width = 6, height = 5)
plot(sft_tab$Power, sft_tab$signed_R2,
     xlab = "Soft threshold power", ylab = "Signed scale-free topology fit R^2",
     type = "b", pch = 19, cex = 0.8, main = "Scale independence")
text(sft_tab$Power[soft_power_label_idx], sft_tab$signed_R2[soft_power_label_idx],
     labels = sft_tab$Power[soft_power_label_idx], cex = 0.75, pos = 3, offset = 0.35)
abline(h = SCALE_FREE_R2_CUTOFF, col = "red", lty = 2)
abline(v = softPower, col = "blue", lty = 2)
dev.off()

pdf(file.path(DIRS$power, "03b_mean_connectivity_readable.pdf"), width = 6, height = 5)
plot(sft_tab$Power, sft_tab$mean.k.,
     xlab = "Soft threshold power", ylab = "Mean connectivity",
     type = "b", pch = 19, cex = 0.8, main = "Mean connectivity")
text(sft_tab$Power[soft_power_label_idx], sft_tab$mean.k.[soft_power_label_idx],
     labels = sft_tab$Power[soft_power_label_idx], cex = 0.75, pos = 3, offset = 0.35)
abline(v = softPower, col = "blue", lty = 2)
dev.off()

## =========================
## 5. Network construction and module detection
## =========================

log_msg("Running blockwiseModules.")
net <- WGCNA::blockwiseModules(datExpr0,
                               power = softPower,
                               maxBlockSize = MAX_BLOCK_SIZE,
                               networkType = NETWORK_TYPE,
                               TOMType = TOM_TYPE,
                               corType = COR_TYPE,
                               maxPOutliers = MAX_P_OUTLIERS,
                               minModuleSize = MIN_MODULE_SIZE,
                               deepSplit = DEEP_SPLIT,
                               pamRespectsDendro = FALSE,
                               mergeCutHeight = MERGE_CUT_HEIGHT,
                               numericLabels = FALSE,
                               saveTOMs = SAVE_TOM,
                               saveTOMFileBase = file.path(DIRS$object, "TOM"),
                               verbose = 3)

moduleColors <- net$colors
names(moduleColors) <- colnames(datExpr0)
MEs <- WGCNA::orderMEs(net$MEs)
moduleLabels <- match(moduleColors, WGCNA::standardColors(200))

module_size <- as.data.frame(table(moduleColors), stringsAsFactors = FALSE)
colnames(module_size) <- c("module", "n_genes")
module_size <- module_size[order(module_size$n_genes, decreasing = TRUE), ]
write.table(module_size, file = file.path(DIRS$module, "01_module_size.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Gene dendrogram and module colors. blockwiseModules may create multiple blocks.
pdf(file.path(DIRS$module, "02_gene_dendrogram_module_colors.pdf"), width = 12, height = 7)
for (b in seq_along(net$dendrograms)) {
  block_genes <- net$blockGenes[[b]]
  WGCNA::plotDendroAndColors(net$dendrograms[[b]],
                             moduleColors[block_genes],
                             "Module colors",
                             dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05,
                             main = paste("Gene dendrogram and module colors - block", b))
}
dev.off()

pdf(file.path(DIRS$module, "03_module_eigengene_clustering.pdf"), width = 7, height = 6)
ME_diss <- 1 - cor(MEs, use = "p")
ME_tree <- hclust(as.dist(ME_diss), method = "average")
plot(ME_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = MERGE_CUT_HEIGHT, col = "red")
dev.off()

saveRDS(net, file = file.path(DIRS$object, "01_blockwiseModules_net.rds"))
saveRDS(list(datExpr = datExpr0, metadata = pd0, datTraits = datTraits0,
             moduleColors = moduleColors, MEs = MEs, softPower = softPower),
        file = file.path(DIRS$object, "02_WGCNA_core_objects.rds"))

## =========================
## 6. Module-trait association
## =========================

nSamples <- nrow(datExpr0)
moduleTraitCor <- WGCNA::bicor(MEs, datTraits0, use = "p",
                               maxPOutliers = MAX_P_OUTLIERS,
                               robustY = FALSE)
moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitFDR <- p_adjust_matrix(moduleTraitPvalue)

write.table(moduleTraitCor, file = file.path(DIRS$trait, "01_module_trait_correlation.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(moduleTraitPvalue, file = file.path(DIRS$trait, "02_module_trait_pvalue.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(moduleTraitFDR, file = file.path(DIRS$trait, "03_module_trait_FDR.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

## Readable module-trait heatmaps.
## Main figure: only AML_status, because Healthy_status is mathematically the inverse of AML_status.
plot_module_trait_heatmap_readable(
  cor_mat = moduleTraitCor,
  p_mat = moduleTraitPvalue,
  fdr_mat = moduleTraitFDR,
  traits = TRAITS_TO_PLOT_IN_MAIN_HEATMAP,
  file = file.path(DIRS$trait, "04_module_trait_heatmap.pdf"),
  title = "Module-trait relationships: AML status"
)

## Supplementary/check figure: all traits.
plot_module_trait_heatmap_readable(
  cor_mat = moduleTraitCor,
  p_mat = moduleTraitPvalue,
  fdr_mat = moduleTraitFDR,
  traits = colnames(datTraits0),
  file = file.path(DIRS$trait, "04b_module_trait_heatmap_all_traits_readable.pdf"),
  title = "Module-trait relationships: all traits"
)

module_trait_long <- do.call(rbind, lapply(colnames(moduleTraitCor), function(tr) {
  data.frame(module_eigengene = rownames(moduleTraitCor),
             module = sub("^ME", "", rownames(moduleTraitCor)),
             trait = tr,
             correlation = moduleTraitCor[, tr],
             pvalue = moduleTraitPvalue[, tr],
             FDR = moduleTraitFDR[, tr],
             stringsAsFactors = FALSE)
}))
write.table(module_trait_long, file = file.path(DIRS$trait, "05_module_trait_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## A cleaner article-style plot than a dense text heatmap: AML-associated modules ranked by correlation.
aml_cor_plot <- module_trait_long[module_trait_long$trait == "AML_status" & module_trait_long$module != "grey", , drop = FALSE]
aml_cor_plot$label <- paste0("r=", sprintf("%.2f", aml_cor_plot$correlation), "; FDR=", fmt_p(aml_cor_plot$FDR))
aml_cor_plot$module <- factor(aml_cor_plot$module, levels = aml_cor_plot$module[order(aml_cor_plot$correlation)])
p_cor_bar <- ggplot(aml_cor_plot, aes(x = module, y = correlation, fill = correlation)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank()) +
  labs(title = "Module association with AML status", x = "Module", y = "Correlation with AML status")
ggsave(file.path(DIRS$trait, "05b_module_AML_correlation_barplot.pdf"),
       p_cor_bar, width = 6.2, height = max(4.5, 0.28 * nrow(aml_cor_plot)))

## Module eigengene group difference using limma; useful because the trait is AML vs Healthy.
ME_mat <- t(MEs)
design_ME <- model.matrix(~ 0 + group, data = pd0)
colnames(design_ME) <- levels(pd0$group)
contrast_ME <- limma::makeContrasts(AML_vs_Healthy = AML - Healthy, levels = design_ME)
fit_ME <- limma::lmFit(ME_mat, design_ME)
fit_ME <- limma::contrasts.fit(fit_ME, contrast_ME)
fit_ME <- limma::eBayes(fit_ME, robust = TRUE)
ME_group <- limma::topTable(fit_ME, coef = "AML_vs_Healthy", number = Inf, adjust.method = "BH", sort.by = "none")
ME_group$module_eigengene <- rownames(ME_group)
ME_group$module <- sub("^ME", "", ME_group$module_eigengene)
colnames(ME_group)[colnames(ME_group) == "logFC"] <- "ME_difference_AML_minus_Healthy"
write.table(ME_group, file = file.path(DIRS$trait, "06_module_eigengene_group_test_limma.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Eigengene boxplots for top modules associated with AML_status
aml_trait <- module_trait_long[module_trait_long$trait == "AML_status", , drop = FALSE]
aml_trait <- aml_trait[order(aml_trait$FDR, aml_trait$pvalue), ]
top_modules_for_boxplot <- head(aml_trait$module, min(12, nrow(aml_trait)))

pdf(file.path(DIRS$trait, "07_top_module_eigengene_boxplots.pdf"), width = 5.5, height = 4.5)
for (mod in top_modules_for_boxplot) {
  me_col <- paste0("ME", mod)
  plot_df <- data.frame(sample_id = rownames(MEs), group = pd0$group, eigengene = MEs[, me_col])
  p <- ggplot(plot_df, aes(x = group, y = eigengene, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.16, size = 1.8, alpha = 0.80) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(title = paste0("Module eigengene: ", mod), x = NULL, y = "Module eigengene")
  print(p)
}
dev.off()

## =========================
## 7. Gene-level MM/GS and DEG integration
## =========================

log_msg("Computing module membership and gene significance.")
modNames <- sub("^ME", "", colnames(MEs))

kME <- WGCNA::signedKME(datExpr0, MEs,
                        outputColumnName = "kME",
                        corFnc = COR_TYPE,
                        corOptions = paste0("use = 'p', maxPOutliers = ", MAX_P_OUTLIERS))
kME_p <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(kME), nSamples))
colnames(kME_p) <- paste0("p.", colnames(kME))

geneTraitSignificance <- as.data.frame(WGCNA::bicor(datExpr0, datTraits0,
                                                     use = "p",
                                                     maxPOutliers = MAX_P_OUTLIERS,
                                                     robustY = FALSE))
GSPvalue <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) <- paste0("GS.", colnames(datTraits0))
colnames(GSPvalue) <- paste0("p.GS.", colnames(datTraits0))

own_kME <- vapply(seq_along(moduleColors), function(i) {
  mod <- moduleColors[i]
  col <- paste0("kME", mod)
  if (col %in% colnames(kME)) kME[i, col] else NA_real_
}, numeric(1))
own_p_kME <- vapply(seq_along(moduleColors), function(i) {
  mod <- moduleColors[i]
  col <- paste0("p.kME", mod)
  if (col %in% colnames(kME_p)) kME_p[i, col] else NA_real_
}, numeric(1))

gene_info <- data.frame(
  gene_symbol = names(moduleColors),
  module = as.character(moduleColors),
  kME_own_module = own_kME,
  p_kME_own_module = own_p_kME,
  GS_AML_status = geneTraitSignificance$GS.AML_status,
  p_GS_AML_status = GSPvalue$p.GS.AML_status,
  GS_Healthy_status = geneTraitSignificance$GS.Healthy_status,
  p_GS_Healthy_status = GSPvalue$p.GS.Healthy_status,
  stringsAsFactors = FALSE
)

gene_info <- cbind(gene_info, kME, kME_p)

gene_info <- merge(gene_info, deg_anno, by = "gene_symbol", all.x = TRUE, sort = FALSE)
gene_info$regulation[is.na(gene_info$regulation)] <- "Not_significant"
gene_info$is_DEG <- gene_info$regulation %in% c("Up_in_AML", "Down_in_AML")
if (!"adj.P.Val" %in% colnames(gene_info)) gene_info$adj.P.Val <- NA_real_
if (!"logFC" %in% colnames(gene_info)) gene_info$logFC <- NA_real_

fdr_for_score <- ifelse(is.na(gene_info$adj.P.Val), 1, pmax(gene_info$adj.P.Val, .Machine$double.xmin))
gene_info$hub_score <- abs(gene_info$kME_own_module) * abs(gene_info$GS_AML_status) * (-log10(fdr_for_score))
gene_info$hub_score[!is.finite(gene_info$hub_score)] <- NA_real_

gene_info <- gene_info[order(gene_info$module, -gene_info$hub_score, -abs(gene_info$kME_own_module)), ]
write.table(gene_info, file = file.path(DIRS$hub, "01_gene_module_membership_GS_DEG.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Per-module gene lists with DEG information
module_gene_dir <- file.path(DIRS$module, "module_gene_lists")
dir.create(module_gene_dir, recursive = TRUE, showWarnings = FALSE)
for (mod in sort(unique(gene_info$module))) {
  x <- gene_info[gene_info$module == mod, , drop = FALSE]
  write.table(x, file = file.path(module_gene_dir, paste0("module_", mod, "_genes.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

## Hub gene tables
non_grey <- gene_info$module != "grey"
hub_candidates <- gene_info[non_grey &
                              abs(gene_info$kME_own_module) >= HUB_KME_CUTOFF &
                              abs(gene_info$GS_AML_status) >= HUB_GS_CUTOFF, , drop = FALSE]
hub_candidates <- hub_candidates[order(hub_candidates$module, -hub_candidates$hub_score), ]
write.table(hub_candidates, file = file.path(DIRS$hub, "02_hub_candidates_by_kME_GS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hub_deg <- hub_candidates[hub_candidates$is_DEG &
                            !is.na(hub_candidates$adj.P.Val) &
                            hub_candidates$adj.P.Val < DEG_FDR_CUTOFF &
                            !is.na(hub_candidates$logFC) &
                            abs(hub_candidates$logFC) >= DEG_LOGFC_CUTOFF, , drop = FALSE]
write.table(hub_deg, file = file.path(DIRS$hub, "03_high_confidence_DEG_hub_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

top_hubs <- do.call(rbind, lapply(split(gene_info[gene_info$module != "grey", ], gene_info$module[gene_info$module != "grey"]), function(x) {
  x <- x[order(-x$hub_score, -abs(x$kME_own_module)), , drop = FALSE]
  head(x, min(TOP_HUB_PER_MODULE, nrow(x)))
}))
write.table(top_hubs, file = file.path(DIRS$hub, "04_top_hubs_per_module.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## MM vs GS scatter plots for key modules
key_by_trait <- aml_trait$module[aml_trait$FDR < KEY_MODULE_FDR_CUTOFF]
if (length(key_by_trait) == 0) key_by_trait <- head(aml_trait$module, min(6, nrow(aml_trait)))
key_by_trait <- unique(key_by_trait[key_by_trait != "grey"])

pdf(file.path(DIRS$hub, "05_MM_vs_GS_key_modules.pdf"), width = 6, height = 5)
for (mod in key_by_trait) {
  x <- gene_info[gene_info$module == mod, , drop = FALSE]
  if (nrow(x) < 5) next
  p <- ggplot(x, aes(x = abs(kME_own_module), y = abs(GS_AML_status), color = regulation)) +
    geom_point(alpha = 0.7, size = 1.7) +
    geom_vline(xintercept = HUB_KME_CUTOFF, linetype = "dashed") +
    geom_hline(yintercept = HUB_GS_CUTOFF, linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(title = paste0("MM vs GS: ", mod, " module"),
         x = paste0("|kME| in ", mod, " module"),
         y = "|GS for AML status|", color = NULL)
  print(p)
}
dev.off()

## =========================
## 8. Module-DEG enrichment and module summary
## =========================

module_deg_enrich <- do.call(rbind, lapply(sort(unique(gene_info$module)), function(mod) {
  if (mod == "grey") return(NULL)
  in_mod <- gene_info$module == mod
  is_deg <- gene_info$is_DEG
  is_up <- gene_info$regulation == "Up_in_AML"
  is_down <- gene_info$regulation == "Down_in_AML"
  fisher_one <- function(target) {
    mat <- matrix(c(sum(in_mod & target), sum(in_mod & !target),
                    sum(!in_mod & target), sum(!in_mod & !target)),
                  nrow = 2, byrow = TRUE)
    ft <- fisher.test(mat, alternative = "greater")
    c(odds_ratio = unname(ft$estimate), pvalue = ft$p.value,
      n_in_module_target = sum(in_mod & target), n_target_total = sum(target))
  }
  all_deg <- fisher_one(is_deg)
  up_deg <- fisher_one(is_up)
  down_deg <- fisher_one(is_down)
  data.frame(module = mod,
             n_genes = sum(in_mod),
             DEG_in_module = sum(in_mod & is_deg),
             Up_in_module = sum(in_mod & is_up),
             Down_in_module = sum(in_mod & is_down),
             DEG_odds_ratio = all_deg["odds_ratio"], DEG_pvalue = all_deg["pvalue"],
             Up_odds_ratio = up_deg["odds_ratio"], Up_pvalue = up_deg["pvalue"],
             Down_odds_ratio = down_deg["odds_ratio"], Down_pvalue = down_deg["pvalue"],
             stringsAsFactors = FALSE)
}))
module_deg_enrich$DEG_FDR <- p.adjust(module_deg_enrich$DEG_pvalue, method = "BH")
module_deg_enrich$Up_FDR <- p.adjust(module_deg_enrich$Up_pvalue, method = "BH")
module_deg_enrich$Down_FDR <- p.adjust(module_deg_enrich$Down_pvalue, method = "BH")
module_deg_enrich <- module_deg_enrich[order(module_deg_enrich$DEG_FDR, -module_deg_enrich$DEG_odds_ratio), ]
write.table(module_deg_enrich, file = file.path(DIRS$trait, "08_module_DEG_enrichment_fisher.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_deg_bar <- ggplot(module_deg_enrich, aes(x = reorder(module, DEG_in_module), y = DEG_in_module, fill = module)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(title = "DEGs per WGCNA module", x = "Module", y = "Number of DEGs")
ggsave(file.path(DIRS$trait, "09_module_DEG_count_barplot.pdf"), p_deg_bar, width = 7, height = max(4, 0.22 * nrow(module_deg_enrich)))

## Build integrated module summary table
aml_cor <- module_trait_long[module_trait_long$trait == "AML_status", c("module", "correlation", "pvalue", "FDR")]
colnames(aml_cor) <- c("module", "cor_AML_status", "p_AML_status", "FDR_AML_status")
me_group_small <- ME_group[, c("module", "ME_difference_AML_minus_Healthy", "P.Value", "adj.P.Val")]
colnames(me_group_small) <- c("module", "ME_difference_AML_minus_Healthy", "ME_group_pvalue", "ME_group_FDR")
summary_df <- merge(module_size, aml_cor, by = "module", all.x = TRUE)
summary_df <- merge(summary_df, me_group_small, by = "module", all.x = TRUE)
summary_df <- merge(summary_df, module_deg_enrich, by = c("module", "n_genes"), all.x = TRUE)

top_hub_names <- tapply(top_hubs$gene_symbol, top_hubs$module, function(z) paste(head(z, 10), collapse = ";"))
summary_df$top_hub_genes <- top_hub_names[summary_df$module]
summary_df <- summary_df[order(summary_df$FDR_AML_status, summary_df$ME_group_FDR, summary_df$DEG_FDR), ]
write.table(summary_df, file = file.path(DIRS$trait, "10_integrated_module_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## =========================
## 9. Optional module ORA for key modules
## =========================

if (RUN_MODULE_ORA) {
  can_ora <- requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)
  if (!can_ora) {
    log_msg("Skipping ORA: clusterProfiler and/or org.Hs.eg.db not installed.")
  } else if (!"entrez_id" %in% colnames(gene_info)) {
    log_msg("Skipping ORA: entrez_id column not available in DEG annotation.")
  } else {
    log_msg("Running optional GO BP ORA for key modules.")
    gene_info$entrez_id_clean <- sub("///.*$", "", as.character(gene_info$entrez_id))
    universe <- unique(gene_info$entrez_id_clean[!is.na(gene_info$entrez_id_clean) & gene_info$entrez_id_clean != "" & gene_info$entrez_id_clean != "---"])
    ora_modules <- unique(c(key_by_trait, head(summary_df$module[summary_df$module != "grey"], 5)))
    ora_modules <- ora_modules[ora_modules != "grey"]
    ora_res <- list()
    for (mod in ora_modules) {
      entrez <- unique(gene_info$entrez_id_clean[gene_info$module == mod])
      entrez <- entrez[!is.na(entrez) & entrez != "" & entrez != "---"]
      if (length(entrez) < 10) next
      ego <- try(clusterProfiler::enrichGO(gene = entrez,
                                           universe = universe,
                                           OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                           keyType = "ENTREZID",
                                           ont = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = ORA_P_CUTOFF,
                                           qvalueCutoff = ORA_Q_CUTOFF,
                                           readable = TRUE), silent = TRUE)
      if (!inherits(ego, "try-error") && !is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        tab <- as.data.frame(ego)
        tab$module <- mod
        ora_res[[mod]] <- tab
      }
    }
    if (length(ora_res) > 0) {
      ora_tab <- do.call(rbind, ora_res)
      ora_tab <- ora_tab[, c("module", setdiff(colnames(ora_tab), "module")), drop = FALSE]
      write.table(ora_tab, file = file.path(DIRS$enrich, "01_key_module_GO_BP_ORA.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      log_msg("No significant GO BP ORA terms found for selected modules.")
    }
  }
}

## =========================
## 10. Cytoscape export for key modules
## =========================

if (RUN_CYTOSCAPE_EXPORT) {
  cyto_modules <- summary_df$module[summary_df$module != "grey"]
  cyto_modules <- unique(head(cyto_modules, min(KEY_MODULE_N_FOR_CYTOSCAPE, length(cyto_modules))))
  log_msg("Exporting Cytoscape files for modules: ", paste(cyto_modules, collapse = ", "))
  for (mod in cyto_modules) {
    x <- gene_info[gene_info$module == mod, , drop = FALSE]
    x <- x[order(-x$hub_score, -abs(x$kME_own_module)), ]
    genes <- head(x$gene_symbol, min(CYTO_TOP_GENES_PER_MODULE, nrow(x)))
    genes <- genes[genes %in% colnames(datExpr0)]
    if (length(genes) < 5) next
    TOM_mod <- WGCNA::TOMsimilarityFromExpr(datExpr0[, genes, drop = FALSE],
                                            power = softPower,
                                            networkType = NETWORK_TYPE,
                                            TOMType = TOM_TYPE,
                                            corType = COR_TYPE,
                                            maxPOutliers = MAX_P_OUTLIERS,
                                            verbose = 0)
    dimnames(TOM_mod) <- list(genes, genes)
    TOM_mod[lower.tri(TOM_mod, diag = TRUE)] <- NA
    edge_idx <- which(TOM_mod >= CYTO_EDGE_THRESHOLD, arr.ind = TRUE)
    if (nrow(edge_idx) > 0) {
      edge_df <- data.frame(from = rownames(TOM_mod)[edge_idx[, 1]],
                            to = colnames(TOM_mod)[edge_idx[, 2]],
                            weight = TOM_mod[edge_idx],
                            module = mod,
                            stringsAsFactors = FALSE)
      edge_df <- edge_df[order(-edge_df$weight), , drop = FALSE]
      if (nrow(edge_df) > CYTO_MAX_EDGES) edge_df <- head(edge_df, CYTO_MAX_EDGES)
    } else {
      edge_df <- data.frame(from = character(0), to = character(0), weight = numeric(0), module = character(0))
    }
    node_df <- x[x$gene_symbol %in% genes,
                 intersect(c("gene_symbol", "module", "kME_own_module", "GS_AML_status", "hub_score",
                             "logFC", "adj.P.Val", "regulation", "is_DEG", "entrez_id"), colnames(x)),
                 drop = FALSE]
    write.table(edge_df, file = file.path(DIRS$cyto, paste0("module_", mod, "_edges.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(node_df, file = file.path(DIRS$cyto, paste0("module_", mod, "_nodes.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

## =========================
## 11. Save final objects and session info
## =========================

saveRDS(list(datExpr = datExpr0,
             metadata = pd0,
             datTraits = datTraits0,
             moduleColors = moduleColors,
             MEs = MEs,
             softPower = softPower,
             module_summary = summary_df,
             gene_info = gene_info,
             module_trait = module_trait_long,
             module_DEG_enrichment = module_deg_enrich),
        file = file.path(DIRS$object, "03_final_WGCNA_integrated_results.rds"))

writeLines(capture.output(sessionInfo()), con = file.path(OUTDIR, "99_sessionInfo.txt"))

log_msg("WGCNA finished. Results saved to: ", normalizePath(OUTDIR))
