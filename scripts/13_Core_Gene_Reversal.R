############################################################
## AML 62-gene core directional drug reversal
## Final refactor based on work 5.6 critique
##
## Preserved formula:
## 1) non-zero rank normalization of |bulk logFC| and |scRNA logFC|
## 2) geometric mean integration
## 3) separate within-UP / within-DOWN weight normalization
## 4) equal 0.5 / 0.5 contribution of UP and DOWN axes
##
## Rebuilt components:
## A) raw-count hard QC: no round(), no negative-to-zero repair
## B) per cell-line x drug DESeq2 model with biological replicates
## C) optional batch-adjusted design
## D) DESeq2-shrunken log2FC used as drug effect
## E) explicit control-pairing checks
## F) gene-label alignment permutation test renamed q_align
## G) reduced manuscript-oriented figure set
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
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(1234)

############################################################
## 0. Packages
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggrepel)
  library(patchwork)
})

############################################################
## 1. Paths
############################################################

IN_CORE <- REVERSAL_CORE_SIGNATURE_FILE
IN_EXPR <- REVERSAL_EXPRESSION_FILE
IN_META <- REVERSAL_METADATA_FILE

OUT_DIR <- REVERSAL_RESULTS_DIR
TABLE_DIR <- file.path(OUT_DIR, "tables")
PDF_DIR <- file.path(OUT_DIR, "figures", "pdf")
PNG_DIR <- file.path(OUT_DIR, "figures", "png")

for (dd in c(OUT_DIR, TABLE_DIR, PDF_DIR, PNG_DIR)) {
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
}

############################################################
## 2. Analysis settings
############################################################

FOCUS_CELL_LINES <- c("OCI-AML3", "MV4-11")
DRUGS_OF_INTEREST <- c("albendazole", "mebendazole")

## Gene-label alignment permutation
N_PERM <- 10000
PERM_SEED <- 20260707


## DESeq2 settings
MIN_TOTAL_COUNT_PER_COMPARISON <- 10
MIN_FINITE_CORE_GENES_PER_COMPARISON <- 50
MIN_REPLICATES_PER_GROUP <- 3

## Preferred adaptive shrinkage estimator for ranking/visualization.
## The script stops rather than silently switching estimators.
LFC_SHRINK_TYPE <- "apeglm"

if (identical(LFC_SHRINK_TYPE, "apeglm") &&
    !requireNamespace("apeglm", quietly = TRUE)) {
  stop(
    "LFC_SHRINK_TYPE='apeglm' requires the Bioconductor package apeglm. ",
    "Install it explicitly (e.g. BiocManager::install('apeglm')) and rerun. ",
    "The script will not silently switch shrinkage estimators."
  )
}

## Same-batch control policy
## Your drug RNA-seq experiment was performed in the same batch
## (simultaneous treatment, library preparation, sequencing, and count generation).
## Therefore the main DESeq2 model is fixed as: design = ~ condition_model.
## No batch term is added.
##
## Optional matching columns below are used ONLY to select matched controls
## if such columns exist. They are NOT used as DESeq2 batch covariates.
PAIRING_COLUMN_CANDIDATES <- c(
  "comparison_id",
  "control_group"
)

## TRUE is appropriate here because you confirmed that contemporaneous
## same-batch controls apply within each cell line.
ALLOW_SHARED_CONTROLS_WITHOUT_PAIRING <- TRUE

############################################################
## 3. Colors
############################################################

direction_cols <- c("core_up" = "#B2182B", "core_down" = "#2166AC")
drug_cols <- c("albendazole" = "#1B9E77", "mebendazole" = "#7570B3")
component_cols <- c("R_UP" = "#B2182B", "R_DOWN" = "#2166AC", "R_TOTAL" = "#333333")
reversal_bg <- "#EAF4EA"
reinforce_bg <- "#F8EAEA"

############################################################
## 4. Helper functions
############################################################

message2 <- function(...) message("[AML-reversal] ", ...)

check_file <- function(x) {
  if (!file.exists(x)) stop("Cannot find file: ", x)
}

harmonize_gene <- function(x) {
  x <- toupper(stringr::str_trim(as.character(x)))
  x <- dplyr::case_when(
    x == "FYB" ~ "FYB1",
    x == "ATP5E" ~ "ATP5F1E",
    TRUE ~ x
  )
  x
}

standardize_cell_line <- function(x) {
  x0 <- as.character(x)
  key <- toupper(gsub("[^A-Z0-9]", "", x0))
  dplyr::case_when(
    key == "OCIAML3" ~ "OCI-AML3",
    key == "MV411" ~ "MV4-11",
    key == "HS5" ~ "HS-5",
    key == "ML2" ~ "ML-2",
    key == "SKM1" ~ "SKM-1",
    TRUE ~ x0
  )
}

standardize_condition <- function(condition, drug) {
  condition <- tolower(as.character(condition))
  drug <- tolower(as.character(drug))
  condition <- stringr::str_replace_all(condition, "-", "_")
  drug <- stringr::str_replace_all(drug, "-", "_")
  
  dplyr::case_when(
    condition %in% c("control", "ctrl", "vehicle", "dmso") ~ "control",
    drug %in% c("control", "ctrl", "vehicle", "dmso") ~ "control",
    condition %in% c("albendazole", "abz") ~ "albendazole",
    drug %in% c("albendazole", "abz") ~ "albendazole",
    condition %in% c("mebendazole", "mbz") ~ "mebendazole",
    drug %in% c("mebendazole", "mbz") ~ "mebendazole",
    TRUE ~ condition
  )
}

rank01 <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  n_ok <- sum(ok)
  if (n_ok == 0) return(out)
  r <- rank(x[ok], ties.method = "average")
  out[ok] <- (r - 0.5) / n_ok
  out
}

calc_weighted_reversal <- function(drug_logfc, direction, weight) {
  ok <- is.finite(drug_logfc) & is.finite(weight) &
    direction %in% c("core_up", "core_down")
  
  drug_logfc <- drug_logfc[ok]
  direction <- direction[ok]
  weight <- weight[ok]
  
  up_idx <- direction == "core_up"
  down_idx <- direction == "core_down"
  
  if (sum(up_idx) == 0 || sum(down_idx) == 0) {
    return(c(R_UP = NA_real_, R_DOWN = NA_real_, R_TOTAL = NA_real_))
  }
  
  w_up <- weight[up_idx] / sum(weight[up_idx])
  w_down <- weight[down_idx] / sum(weight[down_idx])
  
  R_UP <- sum(w_up * (-drug_logfc[up_idx]))
  R_DOWN <- sum(w_down * drug_logfc[down_idx])
  R_TOTAL <- 0.5 * R_UP + 0.5 * R_DOWN
  
  c(R_UP = R_UP, R_DOWN = R_DOWN, R_TOTAL = R_TOTAL)
}

theme_sci <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.key = element_blank(),
      plot.margin = margin(5, 6, 5, 5)
    )
}

save_gg <- function(p, filename, width = 6.5, height = 5.0) {
  ggsave(
    filename = file.path(PDF_DIR, paste0(filename, ".pdf")),
    plot = p, width = width, height = height,
    device = "pdf", useDingbats = FALSE
  )
  ggsave(
    filename = file.path(PNG_DIR, paste0(filename, ".png")),
    plot = p, width = width, height = height,
    dpi = 600, bg = "white"
  )
}

make_clean_plot <- function(plot) {
  clean_theme <- ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank()
  )
  
  if (inherits(plot, "patchwork")) {
    return(
      plot +
        patchwork::plot_annotation(title = NULL, subtitle = NULL, caption = NULL) &
        clean_theme
    )
  }
  
  plot + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) + clean_theme
}

## Fit DESeq2 to an already selected comparison matrix.
## condition_model must have levels control, treated.
fit_deseq_matrix <- function(
    count_sub,
    condition_model,
    batch_model = NULL,
    shrink_type = LFC_SHRINK_TYPE,
    min_total_count = MIN_TOTAL_COUNT_PER_COMPARISON
) {
  if (ncol(count_sub) != length(condition_model)) {
    stop("count_sub columns and condition_model length differ.")
  }
  
  condition_model <- factor(condition_model, levels = c("control", "treated"))
  if (any(is.na(condition_model))) stop("Invalid condition_model values.")
  
  coldata <- data.frame(
    row.names = colnames(count_sub),
    condition_model = condition_model
  )
  
  if (!is.null(batch_model) && length(batch_model) != ncol(count_sub)) {
    stop("batch_model length differs from count_sub columns.")
  }
  
  use_batch <- !is.null(batch_model) && length(unique(batch_model)) > 1
  
  if (use_batch) {
    coldata$batch_model <- factor(batch_model)
    design_formula <- ~ batch_model + condition_model
  } else {
    design_formula <- ~ condition_model
  }
  
  ## Fail early on confounded/non-full-rank designs instead of relying on
  ## a later opaque DESeq2 error.
  design_matrix <- stats::model.matrix(design_formula, data = coldata)
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop(
      "DESeq2 design matrix is not full rank. Treatment and batch/pairing ",
      "are confounded for this comparison. Review matched controls/batches."
    )
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_sub,
    colData = coldata,
    design = design_formula
  )
  
  keep <- rowSums(DESeq2::counts(dds)) >= min_total_count
  if (sum(keep) < 100) {
    stop("Too few genes remain after count prefiltering in a comparison.")
  }
  dds <- dds[keep, ]
  
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  coef_name <- grep(
    "condition_model_treated_vs_control",
    DESeq2::resultsNames(dds),
    value = TRUE,
    fixed = TRUE
  )
  
  if (length(coef_name) != 1) {
    stop(
      "Cannot uniquely identify treated-vs-control coefficient. resultsNames: ",
      paste(DESeq2::resultsNames(dds), collapse = ", ")
    )
  }
  
  res_raw <- DESeq2::results(dds, name = coef_name)
  res_shrunk <- DESeq2::lfcShrink(
    dds,
    coef = coef_name,
    type = shrink_type,
    quiet = TRUE
  )
  
  out <- tibble::tibble(
    gene = rownames(res_shrunk),
    baseMean = as.numeric(res_shrunk$baseMean),
    drug_induced_logFC = as.numeric(res_shrunk$log2FoldChange),
    shrunken_lfcSE = as.numeric(res_shrunk$lfcSE),
    unshrunk_log2FC = as.numeric(res_raw$log2FoldChange),
    deseq2_pvalue = as.numeric(res_raw$pvalue),
    deseq2_padj = as.numeric(res_raw$padj)
  )
  
  list(
    table = out,
    dds = dds,
    coefficient = coef_name,
    design = paste(deparse(design_formula), collapse = " "),
    use_batch = use_batch
  )
}

############################################################
## 5. Check and read input files
############################################################

invisible(lapply(c(IN_CORE, IN_EXPR, IN_META), check_file))

core_raw <- readr::read_tsv(IN_CORE, show_col_types = FALSE)
expr_raw <- readr::read_tsv(IN_EXPR, show_col_types = FALSE)
meta_raw <- readr::read_tsv(IN_META, show_col_types = FALSE)

############################################################
## 6. Core signature and disease-effect weights
############################################################

required_core_cols <- c("gene", "core_direction", "bulk_logFC", "sc_logFC")
missing_core_cols <- setdiff(required_core_cols, colnames(core_raw))
if (length(missing_core_cols) > 0) {
  stop("Missing columns in final core file: ", paste(missing_core_cols, collapse = ", "))
}

core <- core_raw %>%
  dplyr::transmute(
    gene = harmonize_gene(gene),
    direction = as.character(core_direction),
    bulk_logFC = as.numeric(bulk_logFC),
    sc_logFC = as.numeric(sc_logFC)
  ) %>%
  dplyr::filter(
    !is.na(gene), gene != "",
    direction %in% c("core_up", "core_down"),
    is.finite(bulk_logFC), is.finite(sc_logFC)
  )

if (anyDuplicated(core$gene) > 0) {
  dup_core <- core %>% dplyr::count(gene) %>% dplyr::filter(n > 1)
  readr::write_csv(dup_core, file.path(TABLE_DIR, "00_duplicate_core_genes.csv"))
  stop("Duplicated genes found in core signature. No arbitrary row was retained.")
}

if (nrow(core) != 62) {
  stop("Final core must contain exactly 62 unique genes; observed ", nrow(core), ".")
}

## Direction consistency QC
direction_qc <- core %>%
  dplyr::mutate(
    bulk_direction_ok = dplyr::case_when(
      direction == "core_up" ~ bulk_logFC > 0,
      direction == "core_down" ~ bulk_logFC < 0,
      TRUE ~ FALSE
    ),
    sc_direction_ok = dplyr::case_when(
      direction == "core_up" ~ sc_logFC > 0,
      direction == "core_down" ~ sc_logFC < 0,
      TRUE ~ FALSE
    )
  )

readr::write_csv(direction_qc, file.path(TABLE_DIR, "00_core_direction_QC.csv"))
if (any(!direction_qc$bulk_direction_ok) || any(!direction_qc$sc_direction_ok)) {
  stop("Core direction inconsistency detected. See 00_core_direction_QC.csv")
}

## Preserved non-zero rank-weight formula
core_weight <- core %>%
  dplyr::mutate(
    abs_bulk_logFC = abs(bulk_logFC),
    abs_sc_logFC = abs(sc_logFC),
    B_rank01 = rank01(abs_bulk_logFC),
    S_rank01 = rank01(abs_sc_logFC),
    disease_effect_magnitude = sqrt(B_rank01 * S_rank01)
  ) %>%
  dplyr::group_by(direction) %>%
  dplyr::mutate(
    disease_weight = disease_effect_magnitude / sum(disease_effect_magnitude)
  ) %>%
  dplyr::ungroup()

weight_qc <- core_weight %>%
  dplyr::group_by(direction) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    sum_weight = sum(disease_weight),
    min_weight = min(disease_weight),
    max_weight = max(disease_weight),
    .groups = "drop"
  )

readr::write_csv(weight_qc, file.path(TABLE_DIR, "01_weight_QC.csv"))
if (any(abs(weight_qc$sum_weight - 1) > 1e-10)) stop("UP/DOWN weight normalization failed.")
if (any(!is.finite(core_weight$disease_weight) | core_weight$disease_weight <= 0)) {
  stop("At least one disease weight is non-finite or <= 0.")
}

core_weight <- core_weight %>% dplyr::arrange(direction, dplyr::desc(disease_weight))
readr::write_csv(core_weight, file.path(TABLE_DIR, "02_cross_bulk_scRNA_disease_effect_weights.csv"))

############################################################
## 7. Metadata cleaning and explicit control-pairing policy
############################################################

required_meta_cols <- c("sample", "cell_line", "condition", "drug")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_raw))
if (length(missing_meta_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}

meta <- meta_raw %>%
  dplyr::mutate(
    sample = as.character(sample),
    cell_line = standardize_cell_line(cell_line),
    condition_std = standardize_condition(condition, drug)
  )

if (anyDuplicated(meta$sample) > 0) stop("Duplicated sample names found in metadata.")
if (any(is.na(meta$sample) | stringr::str_trim(meta$sample) == "")) stop("Blank sample names found in metadata.")

pairing_col <- intersect(PAIRING_COLUMN_CANDIDATES, colnames(meta))
pairing_col <- if (length(pairing_col) > 0) pairing_col[[1]] else NA_character_

control_policy <- tibble::tibble(
  setting = c(
    "detected_pairing_column",
    "allow_shared_controls_without_pairing",
    "minimum_replicates_per_group"
  ),
  value = c(
    ifelse(is.na(pairing_col), "NONE", pairing_col),
    as.character(ALLOW_SHARED_CONTROLS_WITHOUT_PAIRING),
    as.character(MIN_REPLICATES_PER_GROUP)
  )
)
readr::write_csv(control_policy, file.path(TABLE_DIR, "00_control_pairing_policy.csv"))

if (is.na(pairing_col) && !ALLOW_SHARED_CONTROLS_WITHOUT_PAIRING) {
  stop(
    "No explicit control-matching column was found in sample_metadata.txt. ",
    "The script will not silently pool all controls unless ALLOW_SHARED_CONTROLS_WITHOUT_PAIRING is TRUE. Add one of: ",
    paste(PAIRING_COLUMN_CANDIDATES, collapse = ", "),
    "; or keep ALLOW_SHARED_CONTROLS_WITHOUT_PAIRING <- TRUE only if the same ",
    "contemporaneous same-batch DMSO controls genuinely apply within each cell line."
  )
}

############################################################
## 8. Expression matrix: raw-count hard QC
##
## Critical rule:
## - no round()
## - no negative-to-zero repair
## - non-negative integer counts are required
############################################################

if (!"gene" %in% colnames(expr_raw)) {
  stop("expr_matrix.txt must contain a column named 'gene'.")
}

expr_raw <- expr_raw %>%
  dplyr::mutate(
    .input_row = dplyr::row_number(),
    gene = harmonize_gene(gene)
  )

invalid_gene_rows <- expr_raw %>%
  dplyr::filter(is.na(gene) | gene == "" | gene == "NA")
if (nrow(invalid_gene_rows) > 0) {
  readr::write_csv(
    invalid_gene_rows,
    file.path(TABLE_DIR, "00_invalid_gene_name_rows_in_expression.csv")
  )
  warning("Removed ", nrow(invalid_gene_rows), " rows with invalid gene names.")
}

expr_raw <- expr_raw %>%
  dplyr::filter(!is.na(gene), gene != "", gene != "NA") %>%
  dplyr::select(-.input_row) %>%
  dplyr::mutate(
    dplyr::across(-gene, ~ suppressWarnings(as.numeric(.x)))
  )

expr_numeric_mat <- as.matrix(dplyr::select(expr_raw, -gene))
storage.mode(expr_numeric_mat) <- "numeric"

if (any(!is.finite(expr_numeric_mat))) {
  stop(
    "Expression matrix contains ", sum(!is.finite(expr_numeric_mat)),
    " NA/NaN/Inf values after numeric conversion."
  )
}

if (any(expr_numeric_mat < 0)) {
  bad_n <- sum(expr_numeric_mat < 0)
  stop(
    "Expression matrix contains ", bad_n,
    " negative values. DESeq2 raw counts must not be repaired by setting them to zero."
  )
}

integer_deviation <- abs(expr_numeric_mat - round(expr_numeric_mat))
if (any(integer_deviation != 0)) {
  bad_n <- sum(integer_deviation != 0)
  bad_max <- max(integer_deviation)
  stop(
    "Expression matrix is not a non-negative integer count matrix: ",
    bad_n, " values are non-integer; maximum deviation from integer = ",
    signif(bad_max, 5), ". Do not input TPM/FPKM/log-expression values to DESeq2."
  )
}

## Aggregate duplicated gene symbols by sum only after raw-count QC.
expr_agg <- expr_raw %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), sum), .groups = "drop")

############################################################
## 9. Match samples and build exact integer count matrix
############################################################

expr_samples <- setdiff(colnames(expr_agg), "gene")
missing_in_meta <- setdiff(expr_samples, meta$sample)
missing_in_expr <- setdiff(meta$sample, expr_samples)

if (length(missing_in_meta) > 0) {
  stop("Expression samples missing in metadata: ", paste(missing_in_meta, collapse = ", "))
}
if (length(missing_in_expr) > 0) {
  stop("Metadata samples missing in expression: ", paste(missing_in_expr, collapse = ", "))
}

sample_match_idx <- base::match(expr_samples, meta$sample)
if (anyNA(sample_match_idx)) stop("Failed to reorder metadata to expression columns.")
meta <- tibble::as_tibble(meta[sample_match_idx, , drop = FALSE])
stopifnot(identical(as.character(meta$sample), as.character(expr_samples)))

if (anyNA(expr_agg$gene) || any(expr_agg$gene == "")) stop("Invalid gene names remain after aggregation.")
if (anyDuplicated(expr_agg$gene) > 0) stop("Duplicated genes remain after aggregation.")

count_mat <- as.matrix(dplyr::select(expr_agg, -gene))
rownames(count_mat) <- expr_agg$gene
storage.mode(count_mat) <- "integer"

if (any(count_mat < 0)) stop("Internal QC failed: negative counts after exact integer conversion.")
if (any(colSums(count_mat) <= 0)) stop("At least one sample has zero library size.")

############################################################
## 10. Build matched comparison metadata
############################################################

build_comparison_meta <- function(meta, cell_line_i, drug_i, pairing_col = NA_character_) {
  treated <- meta %>%
    dplyr::filter(cell_line == cell_line_i, condition_std == drug_i)
  controls <- meta %>%
    dplyr::filter(cell_line == cell_line_i, condition_std == "control")
  
  if (nrow(treated) == 0) stop("No treated samples for ", cell_line_i, " / ", drug_i)
  if (nrow(controls) == 0) stop("No control samples for ", cell_line_i)
  
  if (!is.na(pairing_col)) {
    treated_pair <- as.character(treated[[pairing_col]])
    if (any(is.na(treated_pair) | stringr::str_trim(treated_pair) == "")) {
      stop("Missing ", pairing_col, " values in treated samples for ", cell_line_i, " / ", drug_i)
    }
    
    treated_pair_levels <- unique(treated_pair)
    controls <- controls %>%
      dplyr::filter(as.character(.data[[pairing_col]]) %in% treated_pair_levels)
    
    if (nrow(controls) == 0) {
      stop("No matched controls after applying ", pairing_col, " for ", cell_line_i, " / ", drug_i)
    }
    
    control_pair_levels <- unique(as.character(controls[[pairing_col]]))
    missing_control_pairs <- setdiff(treated_pair_levels, control_pair_levels)
    if (length(missing_control_pairs) > 0) {
      stop(
        "Treated pairing levels without matched controls for ", cell_line_i, " / ", drug_i,
        ": ", paste(missing_control_pairs, collapse = ", ")
      )
    }
  }
  
  comp <- dplyr::bind_rows(controls, treated) %>%
    dplyr::mutate(
      comparison_cell_line = cell_line_i,
      comparison_drug = drug_i,
      condition_model = ifelse(condition_std == "control", "control", "treated")
    )
  
  n_ctrl <- sum(comp$condition_model == "control")
  n_treat <- sum(comp$condition_model == "treated")
  
  if (n_ctrl < MIN_REPLICATES_PER_GROUP || n_treat < MIN_REPLICATES_PER_GROUP) {
    stop(
      "Insufficient biological replicates for ", cell_line_i, " / ", drug_i,
      ": control n=", n_ctrl, ", treated n=", n_treat,
      ". Minimum per group = ", MIN_REPLICATES_PER_GROUP
    )
  }
  
  comp
}

comparison_specs <- tidyr::expand_grid(
  cell_line = FOCUS_CELL_LINES,
  drug = DRUGS_OF_INTEREST
)

comparison_meta_list <- vector("list", nrow(comparison_specs))
comparison_qc_list <- vector("list", nrow(comparison_specs))

for (i in seq_len(nrow(comparison_specs))) {
  cell_i <- comparison_specs$cell_line[[i]]
  drug_i <- comparison_specs$drug[[i]]
  
  comp_i <- build_comparison_meta(meta, cell_i, drug_i, pairing_col)
  comparison_meta_list[[i]] <- comp_i
  
  n_batch <- if (!is.na(pairing_col)) {
    dplyr::n_distinct(comp_i[[pairing_col]])
  } else {
    1L
  }
  
  comparison_qc_list[[i]] <- tibble::tibble(
    cell_line = cell_i,
    drug = drug_i,
    n_control = sum(comp_i$condition_model == "control"),
    n_treated = sum(comp_i$condition_model == "treated"),
    pairing_column = ifelse(is.na(pairing_col), "NONE", pairing_col),
    n_pairing_levels = n_batch,
    pooled_controls_without_pairing = is.na(pairing_col)
  )
}

names(comparison_meta_list) <- paste(comparison_specs$cell_line, comparison_specs$drug, sep = "__")
comparison_qc <- dplyr::bind_rows(comparison_qc_list)
readr::write_csv(comparison_qc, file.path(TABLE_DIR, "03_comparison_and_control_pairing_QC.csv"))

## Exact sample-level audit trail for every drug-vs-control comparison.
## This prevents silent control pooling and makes the fitted sample sets reviewable.
matched_sample_manifest <- purrr::imap_dfr(
  comparison_meta_list,
  function(comp_i, comparison_name_i) {
    out <- comp_i %>%
      dplyr::transmute(
        comparison = comparison_name_i,
        cell_line = comparison_cell_line,
        drug = comparison_drug,
        sample = sample,
        role = condition_model,
        original_condition = condition_std
      )
    
    if (!is.na(pairing_col)) {
      out$pairing_column <- pairing_col
      out$pairing_value <- as.character(comp_i[[pairing_col]])
    } else {
      out$pairing_column <- "NONE"
      out$pairing_value <- NA_character_
    }
    
    out
  }
)

readr::write_csv(
  matched_sample_manifest,
  file.path(TABLE_DIR, "03_matched_samples_by_comparison.csv")
)

############################################################
## 11. Per-comparison DESeq2 models and shrunken LFC
############################################################

lfc_result_list <- vector("list", length(comparison_meta_list))
model_qc_list <- vector("list", length(comparison_meta_list))

for (k in seq_along(comparison_meta_list)) {
  comp_i <- comparison_meta_list[[k]]
  cell_i <- comp_i$comparison_cell_line[[1]]
  drug_i <- comp_i$comparison_drug[[1]]
  
  sample_ids <- comp_i$sample
  count_sub <- count_mat[, sample_ids, drop = FALSE]
  
  ## Same-batch experiment: no batch covariate is included.
  ## The fitted model is design = ~ condition_model.
  batch_vec <- NULL
  
  fit_i <- fit_deseq_matrix(
    count_sub = count_sub,
    condition_model = comp_i$condition_model,
    batch_model = batch_vec
  )
  
  lfc_result_list[[k]] <- fit_i$table %>%
    dplyr::mutate(
      cell_line = cell_i,
      drug = drug_i,
      shrink_type = LFC_SHRINK_TYPE,
      model_design = fit_i$design
    )
  
  model_qc_list[[k]] <- tibble::tibble(
    cell_line = cell_i,
    drug = drug_i,
    n_samples = nrow(comp_i),
    n_control = sum(comp_i$condition_model == "control"),
    n_treated = sum(comp_i$condition_model == "treated"),
    use_batch = fit_i$use_batch,
    design = fit_i$design,
    coefficient = fit_i$coefficient,
    shrink_type = LFC_SHRINK_TYPE,
    n_modeled_genes = nrow(fit_i$table)
  )
  
  message2(
    "DESeq2 completed: ", cell_i, " / ", drug_i,
    "; design ", fit_i$design,
    "; shrunken LFC type = ", LFC_SHRINK_TYPE
  )
}

all_deseq_lfc <- dplyr::bind_rows(lfc_result_list)
model_qc <- dplyr::bind_rows(model_qc_list)

readr::write_csv(all_deseq_lfc, file.path(TABLE_DIR, "04_all_genes_DESeq2_shrunken_LFC.csv"))
readr::write_csv(model_qc, file.path(TABLE_DIR, "04_DESeq2_model_QC.csv"))

############################################################
## 12. Merge with 62-gene core and comparison-specific weights
############################################################

drug_logfc_pre <- all_deseq_lfc %>%
  dplyr::inner_join(
    core_weight %>%
      dplyr::select(
        gene, direction, bulk_logFC, sc_logFC,
        B_rank01, S_rank01, disease_effect_magnitude, disease_weight
      ),
    by = "gene"
  )

core_detection_qc <- drug_logfc_pre %>%
  dplyr::group_by(cell_line, drug) %>%
  dplyr::summarise(
    n_core_rows = dplyr::n(),
    n_finite_shrunken_lfc = sum(is.finite(drug_induced_logFC)),
    n_up_finite = sum(direction == "core_up" & is.finite(drug_induced_logFC)),
    n_down_finite = sum(direction == "core_down" & is.finite(drug_induced_logFC)),
    .groups = "drop"
  )

readr::write_csv(core_detection_qc, file.path(TABLE_DIR, "05_core_gene_detection_by_comparison_QC.csv"))

if (any(core_detection_qc$n_finite_shrunken_lfc < MIN_FINITE_CORE_GENES_PER_COMPARISON)) {
  stop(
    "At least one comparison has fewer than ", MIN_FINITE_CORE_GENES_PER_COMPARISON,
    " core genes with finite shrunken LFC. See QC table."
  )
}
if (any(core_detection_qc$n_up_finite == 0 | core_detection_qc$n_down_finite == 0)) {
  stop("At least one comparison lacks finite genes on one disease-direction axis.")
}

## Comparison-specific renormalization within direction after excluding NA LFC.
drug_logfc <- drug_logfc_pre %>%
  dplyr::filter(is.finite(drug_induced_logFC)) %>%
  dplyr::group_by(cell_line, drug, direction) %>%
  dplyr::mutate(
    analysis_weight = disease_weight / sum(disease_weight)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    directional_response = dplyr::case_when(
      direction == "core_up" ~ -drug_induced_logFC,
      direction == "core_down" ~ drug_induced_logFC,
      TRUE ~ NA_real_
    ),
    reversed = directional_response > 0,
    weighted_contribution = 0.5 * analysis_weight * directional_response
  )

readr::write_csv(
  drug_logfc,
  file.path(TABLE_DIR, "06_gene_level_DESeq2_shrunken_weighted_reversal.csv")
)

############################################################
## 13. Weighted reversal scores
############################################################

score_table <- drug_logfc %>%
  dplyr::group_by(cell_line, drug) %>%
  dplyr::group_modify(~ {
    s <- calc_weighted_reversal(
      drug_logfc = .x$drug_induced_logFC,
      direction = .x$direction,
      weight = .x$analysis_weight
    )
    
    up_idx <- .x$direction == "core_up"
    down_idx <- .x$direction == "core_down"
    
    tibble::tibble(
      n_genes = nrow(.x),
      n_up = sum(up_idx),
      n_down = sum(down_idx),
      R_UP = unname(s["R_UP"]),
      R_DOWN = unname(s["R_DOWN"]),
      R_TOTAL = unname(s["R_TOTAL"]),
      reversed_weight_UP = sum(.x$analysis_weight[up_idx] * (.x$drug_induced_logFC[up_idx] < 0)),
      reversed_weight_DOWN = sum(.x$analysis_weight[down_idx] * (.x$drug_induced_logFC[down_idx] > 0)),
      reversed_weight_TOTAL =
        0.5 * sum(.x$analysis_weight[up_idx] * (.x$drug_induced_logFC[up_idx] < 0)) +
        0.5 * sum(.x$analysis_weight[down_idx] * (.x$drug_induced_logFC[down_idx] > 0))
    )
  }) %>%
  dplyr::ungroup()

## Mathematical consistency
contribution_check <- drug_logfc %>%
  dplyr::group_by(cell_line, drug) %>%
  dplyr::summarise(contribution_sum = sum(weighted_contribution), .groups = "drop")

score_table <- score_table %>%
  dplyr::left_join(contribution_check, by = c("cell_line", "drug")) %>%
  dplyr::mutate(absolute_check_difference = abs(R_TOTAL - contribution_sum))

if (any(score_table$absolute_check_difference > 1e-10)) {
  stop("Internal score consistency check failed.")
}

############################################################
## 14. Gene-label alignment permutation test
##
## Interpretation:
## Tests whether observed drug effects align with the predefined
## disease directions/weights more strongly than after random
## reassignment of drug-effect values across the 62-gene labels.
##
## This is NOT DESeq2 differential-expression FDR.
############################################################

set.seed(PERM_SEED)

split_data <- split(
  drug_logfc,
  interaction(drug_logfc$cell_line, drug_logfc$drug, drop = TRUE)
)

perm_result_list <- vector("list", length(split_data))
perm_distribution_list <- vector("list", length(split_data))
names(perm_distribution_list) <- names(split_data)

for (k in seq_along(split_data)) {
  dat <- split_data[[k]]
  
  observed <- calc_weighted_reversal(
    dat$drug_induced_logFC,
    dat$direction,
    dat$analysis_weight
  )["R_TOTAL"]
  
  perm_scores <- replicate(N_PERM, {
    permuted_lfc <- sample(dat$drug_induced_logFC, size = nrow(dat), replace = FALSE)
    calc_weighted_reversal(
      permuted_lfc,
      dat$direction,
      dat$analysis_weight
    )["R_TOTAL"]
  })
  
  perm_scores <- as.numeric(perm_scores)
  alignment_p <- (1 + sum(perm_scores >= observed, na.rm = TRUE)) / (N_PERM + 1)
  null_mean <- mean(perm_scores, na.rm = TRUE)
  null_sd <- sd(perm_scores, na.rm = TRUE)
  null_z <- ifelse(is.finite(null_sd) && null_sd > 0, (observed - null_mean) / null_sd, NA_real_)
  
  perm_result_list[[k]] <- tibble::tibble(
    cell_line = dat$cell_line[[1]],
    drug = dat$drug[[1]],
    R_TOTAL_observed = as.numeric(observed),
    alignment_p = alignment_p,
    alignment_null_mean = null_mean,
    alignment_null_sd = null_sd,
    alignment_null_z = null_z,
    n_permutations = N_PERM
  )
  
  perm_distribution_list[[k]] <- perm_scores
}

alignment_results <- dplyr::bind_rows(perm_result_list) %>%
  dplyr::mutate(
    q_align = p.adjust(alignment_p, method = "BH")
  )

score_table <- score_table %>%
  dplyr::left_join(alignment_results, by = c("cell_line", "drug")) %>%
  dplyr::arrange(cell_line, drug)

readr::write_csv(alignment_results, file.path(TABLE_DIR, "07_gene_label_alignment_permutation_test.csv"))
saveRDS(perm_distribution_list, file.path(TABLE_DIR, "07_alignment_permutation_null_distributions.rds"))

############################################################
## 15. Export final score table
############################################################

readr::write_csv(
  score_table,
  file.path(TABLE_DIR, "09_weighted_directional_reversal_scores_final.csv")
)

############################################################
## 16. Console results
############################################################

cat("\n============================================================\n")
cat("DESeq2-shrunken weighted directional reversal results\n")
cat("============================================================\n")
print(
  score_table %>%
    dplyr::select(
      cell_line, drug, R_UP, R_DOWN, R_TOTAL,
      reversed_weight_TOTAL,
      alignment_p, q_align, alignment_null_z
    )
)

############################################################
## 17. Main Figure A: 2x2 shrunken-LFC reversal scatter
############################################################

N_LABEL_GENES <- 5

weighted_scatter_df <- drug_logfc %>%
  dplyr::mutate(
    signed_disease_effect = dplyr::case_when(
      direction == "core_up" ~ disease_effect_magnitude,
      direction == "core_down" ~ -disease_effect_magnitude,
      TRUE ~ NA_real_
    ),
    drug_label = dplyr::recode(
      drug,
      albendazole = "Albendazole",
      mebendazole = "Mebendazole"
    ),
    cell_line = factor(cell_line, levels = FOCUS_CELL_LINES),
    drug_label = factor(drug_label, levels = c("Albendazole", "Mebendazole"))
  )

scatter_label_df <- weighted_scatter_df %>%
  dplyr::group_by(cell_line, drug_label) %>%
  dplyr::slice_max(order_by = abs(weighted_contribution), n = N_LABEL_GENES, with_ties = FALSE) %>%
  dplyr::ungroup()

scatter_score_anno <- score_table %>%
  dplyr::mutate(
    drug_label = dplyr::recode(
      drug,
      albendazole = "Albendazole",
      mebendazole = "Mebendazole"
    ),
    cell_line = factor(cell_line, levels = FOCUS_CELL_LINES),
    drug_label = factor(drug_label, levels = c("Albendazole", "Mebendazole")),
    score_label = paste0(
      "R[total] = ", sprintf("%.3f", R_TOTAL),
      "\nq[align] = ", format.pval(q_align, digits = 2, eps = 1e-4)
    )
  )

x_lim <- max(abs(weighted_scatter_df$signed_disease_effect), na.rm = TRUE) * 1.10
y_lim <- max(abs(weighted_scatter_df$drug_induced_logFC), na.rm = TRUE) * 1.12
weight_breaks <- c(0.03, 0.06, 0.09, 0.12)
weight_breaks <- weight_breaks[
  weight_breaks <= max(weighted_scatter_df$analysis_weight, na.rm = TRUE)
]
p_scatter <- ggplot(
  weighted_scatter_df,
  aes(x = signed_disease_effect, y = drug_induced_logFC)
) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = reversal_bg, alpha = 0.45) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = reversal_bg, alpha = 0.45) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = reinforce_bg, alpha = 0.30) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = reinforce_bg, alpha = 0.30) +
  geom_hline(yintercept = 0, linewidth = 0.42, linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = 0, linewidth = 0.42, linetype = "dashed", color = "grey45") +
  geom_point(
    aes(size = analysis_weight, fill = direction),
    shape = 21, color = "white", stroke = 0.3, alpha = 0.88
  ) +
  ggrepel::geom_text_repel(
    data = scatter_label_df,
    aes(label = gene),
    size = 2.35, color = "black",
    box.padding = 0.32, point.padding = 0.18,
    min.segment.length = 0, segment.color = "grey55",
    segment.linewidth = 0.25, max.overlaps = Inf,
    seed = 1234, show.legend = FALSE
  ) +
  geom_label(
    data = scatter_score_anno,
    aes(x = Inf, y = Inf, label = score_label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.10,
    size = 2.55, lineheight = 1.02,
    color = "black", fill = scales::alpha("white", 0.82),
    label.size = 0,
    label.padding = grid::unit(0.12, "lines")
  ) +
  facet_grid(cell_line ~ drug_label) +
  scale_fill_manual(
    values = direction_cols,
    breaks = c("core_up", "core_down"),
    labels = c("Core-up", "Core-down"),
    name = "Core direction"
  ) +
  scale_size_area(
    max_size = 6.2,
    breaks = weight_breaks,
    limits = c(0, max(weighted_scatter_df$analysis_weight, na.rm = TRUE)),
    name = "Disease weight"
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        size = 4,
        shape = 21,
        color = "white",
        alpha = 0.90
      )
    ),
    size = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        fill = "grey55",
        color = "grey25",
        alpha = 0.95,
        stroke = 0.25
      )
    )
  ) +
  coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(-y_lim, y_lim), clip = "off") +
  labs(
    x = "Signed cross-omics disease effect\n(Core-down <- 0 -> Core-up)",
    y = "DESeq2-shrunken log2FC (treated vs control)",
    title = "Disease-effect-weighted alignment of drug responses with the AML core program"
  ) +
  theme_sci(base_size = 8.8) +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 7.3),
    legend.key.height = grid::unit(0.55, "cm"),
    legend.key.width = grid::unit(0.55, "cm"),
    legend.spacing.y = grid::unit(0.15, "cm"),
    panel.spacing = grid::unit(0.65, "lines")
  )

save_gg(p_scatter, "Fig_main_A_DESeq2_shrunken_2x2_scatter", width = 7.1, height = 5.15)
save_gg(make_clean_plot(p_scatter), "Fig_main_A_DESeq2_shrunken_2x2_scatter_clean", width = 7.1, height = 5.15)

############################################################
## 17A. Main Figure A individual clean panels
############################################################

panel_list <- weighted_scatter_df %>%
  dplyr::distinct(cell_line, drug_label) %>%
  dplyr::arrange(cell_line, drug_label)

for (ii in seq_len(nrow(panel_list))) {
  cell_i <- as.character(panel_list$cell_line[[ii]])
  drug_i <- as.character(panel_list$drug_label[[ii]])
  
  dat_i <- weighted_scatter_df %>%
    dplyr::filter(
      as.character(cell_line) == cell_i,
      as.character(drug_label) == drug_i
    )
  
  lab_i <- scatter_label_df %>%
    dplyr::filter(
      as.character(cell_line) == cell_i,
      as.character(drug_label) == drug_i
    )
  
  anno_i <- scatter_score_anno %>%
    dplyr::filter(
      as.character(cell_line) == cell_i,
      as.character(drug_label) == drug_i
    )
  
  p_i <- ggplot(
    dat_i,
    aes(x = signed_disease_effect, y = drug_induced_logFC)
  ) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = reversal_bg, alpha = 0.45) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = reversal_bg, alpha = 0.45) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = reinforce_bg, alpha = 0.30) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = reinforce_bg, alpha = 0.30) +
    geom_hline(yintercept = 0, linewidth = 0.42, linetype = "dashed", color = "grey45") +
    geom_vline(xintercept = 0, linewidth = 0.42, linetype = "dashed", color = "grey45") +
    geom_point(
      aes(size = analysis_weight, fill = direction),
      shape = 21, color = "white", stroke = 0.3, alpha = 0.88
    ) +
    ggrepel::geom_text_repel(
      data = lab_i,
      aes(label = gene),
      size = 2.6, color = "black",
      box.padding = 0.32, point.padding = 0.18,
      min.segment.length = 0, segment.color = "grey55",
      segment.linewidth = 0.25, max.overlaps = Inf,
      seed = 1234, show.legend = FALSE
    ) +
    geom_label(
      data = anno_i,
      aes(x = Inf, y = Inf, label = score_label),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = 1.10,
      size = 2.8, lineheight = 1.02,
      color = "black", fill = scales::alpha("white", 0.82),
      label.size = 0,
      label.padding = grid::unit(0.12, "lines")
    ) +
    scale_fill_manual(
      values = direction_cols,
      breaks = c("core_up", "core_down"),
      labels = c("Core-up", "Core-down"),
      name = "Core direction"
    ) +
    scale_size_area(
      max_size = 7.0,
      breaks = weight_breaks,
      limits = c(0, max(weighted_scatter_df$analysis_weight, na.rm = TRUE)),
      name = "Disease weight"
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          size = 4,
          shape = 21,
          color = "white",
          alpha = 0.90
        )
      ),
      size = guide_legend(
        order = 2,
        override.aes = list(
          shape = 21,
          fill = "grey55",
          color = "grey25",
          alpha = 0.95,
          stroke = 0.25
        )
      )
    ) +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(-y_lim, y_lim), clip = "off") +
    labs(
      x = "Signed cross-omics disease effect\n(Core-down <- 0 -> Core-up)",
      y = "DESeq2-shrunken log2FC (treated vs control)",
      title = paste(cell_i, drug_i, sep = " | ")
    ) +
    theme_sci(base_size = 10) +
    theme(
      aspect.ratio = 1,
      legend.position = "right",
      legend.title = element_text(size = 10.5, face = "bold"),
      legend.text = element_text(size = 9.5),
      legend.key.height = grid::unit(0.55, "cm"),
      legend.key.width = grid::unit(0.55, "cm"),
      legend.spacing.y = grid::unit(0.15, "cm")
    )
  
  safe_cell <- gsub("[^A-Za-z0-9]", "", cell_i)
  safe_drug <- gsub("[^A-Za-z0-9]", "", drug_i)
  
  save_gg(
    make_clean_plot(p_i),
    paste0("Fig_main_A_DESeq2_shrunken_scatter_", safe_cell, "_", safe_drug, "_clean"),
    width = 5.4, height = 5.4
  )
}

############################################################
## 18. Main Figure B: R_TOTAL summary
############################################################

score_plot_df <- score_table %>%
  dplyr::mutate(
    condition_label = dplyr::case_when(
      cell_line == "OCI-AML3" & drug == "albendazole" ~ "OCI-AML3 | ABZ",
      cell_line == "OCI-AML3" & drug == "mebendazole" ~ "OCI-AML3 | MBZ",
      cell_line == "MV4-11" & drug == "albendazole" ~ "MV4-11 | ABZ",
      cell_line == "MV4-11" & drug == "mebendazole" ~ "MV4-11 | MBZ",
      TRUE ~ paste(cell_line, drug, sep = " | ")
    ),
    condition_label = factor(
      condition_label,
      levels = c("OCI-AML3 | ABZ", "OCI-AML3 | MBZ", "MV4-11 | ABZ", "MV4-11 | MBZ")
    )
  )

p_score <- ggplot(score_plot_df, aes(x = condition_label, y = R_TOTAL, color = drug)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.42, color = "grey55") +
  geom_point(size = 3.0) +
  scale_color_manual(values = drug_cols, guide = "none") +
  labs(
    x = NULL,
    y = expression(R[total]),
    title = "Weighted reversal score"
  ) +
  theme_sci(base_size = 8.5) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1))

save_gg(p_score, "Fig_main_B_R_TOTAL", width = 3.5, height = 3.15)
save_gg(make_clean_plot(p_score), "Fig_main_B_R_TOTAL_clean", width = 3.5, height = 3.15)

############################################################
## 19. Main Figure C: compact R_UP/R_DOWN/R_TOTAL heatmap
############################################################

architecture_df <- score_table %>%
  dplyr::select(cell_line, drug, R_UP, R_DOWN, R_TOTAL) %>%
  tidyr::pivot_longer(
    cols = c(R_UP, R_DOWN, R_TOTAL),
    names_to = "component", values_to = "score"
  ) %>%
  dplyr::mutate(
    component = factor(component, levels = c("R_UP", "R_DOWN", "R_TOTAL")),
    condition_label = dplyr::case_when(
      cell_line == "OCI-AML3" & drug == "albendazole" ~ "OCI-AML3 | ABZ",
      cell_line == "OCI-AML3" & drug == "mebendazole" ~ "OCI-AML3 | MBZ",
      cell_line == "MV4-11" & drug == "albendazole" ~ "MV4-11 | ABZ",
      cell_line == "MV4-11" & drug == "mebendazole" ~ "MV4-11 | MBZ",
      TRUE ~ paste(cell_line, drug, sep = " | ")
    ),
    condition_label = factor(
      condition_label,
      levels = c("OCI-AML3 | ABZ", "OCI-AML3 | MBZ", "MV4-11 | ABZ", "MV4-11 | MBZ")
    )
  )

architecture_limit <- max(abs(architecture_df$score), na.rm = TRUE)

p_architecture <- ggplot(
  architecture_df,
  aes(x = component, y = condition_label, fill = score)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.3f", score)), size = 2.6, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "grey98", high = "#B2182B",
    midpoint = 0, limits = c(-architecture_limit, architecture_limit),
    name = "Score"
  ) +
  labs(x = NULL, y = NULL, title = "Directional reversal architecture") +
  theme_sci(base_size = 8.5) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 7.2),
    legend.position = "right"
  )

save_gg(p_architecture, "Fig_main_C_reversal_architecture_heatmap", width = 4.0, height = 4.0)
save_gg(make_clean_plot(p_architecture), "Fig_main_C_reversal_architecture_heatmap_clean", width = 4.0, height = 4.0)

############################################################
## 20. Combined manuscript-oriented main figure
############################################################

Fig_main <- patchwork::wrap_plots(
  p_scatter,
  patchwork::wrap_plots(p_score, p_architecture, ncol = 2),
  ncol = 1,
  heights = c(1.70, 1.0)
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "DESeq2-shrunken, disease-effect-weighted reversal of the AML core transcriptional program",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
      plot.tag = element_text(face = "bold", size = 11)
    )
  )

save_gg(Fig_main, "Figure_main_AML_reversal_final", width = 7.1, height = 8.15)

## Explicit clean combined figure: remove titles while retaining A-C panel tags.
Fig_main_clean <- patchwork::wrap_plots(
  make_clean_plot(p_scatter),
  patchwork::wrap_plots(make_clean_plot(p_score), make_clean_plot(p_architecture), ncol = 2),
  ncol = 1,
  heights = c(1.70, 1.0)
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    theme = theme(plot.tag = element_text(face = "bold", size = 11))
  )

save_gg(Fig_main_clean, "Figure_main_AML_reversal_final_clean", width = 7.1, height = 8.15)

############################################################
## 21. Supplementary Figure S1: full 62-gene directional-response heatmap
##
## Shows directional_response, NOT weighted_contribution.
## Thus the color encodes drug effect aligned to disease direction;
## predefined disease weights are not mixed into the color value.
############################################################

gene_order_heat <- core_weight %>%
  dplyr::mutate(direction = factor(direction, levels = c("core_up", "core_down"))) %>%
  dplyr::arrange(direction, dplyr::desc(disease_weight)) %>%
  dplyr::pull(gene)

response_heatmap_df <- drug_logfc %>%
  dplyr::mutate(
    condition_label = dplyr::case_when(
      cell_line == "OCI-AML3" & drug == "albendazole" ~ "OCI-AML3 | ABZ",
      cell_line == "OCI-AML3" & drug == "mebendazole" ~ "OCI-AML3 | MBZ",
      cell_line == "MV4-11" & drug == "albendazole" ~ "MV4-11 | ABZ",
      cell_line == "MV4-11" & drug == "mebendazole" ~ "MV4-11 | MBZ",
      TRUE ~ paste(cell_line, drug, sep = " | ")
    ),
    condition_label = factor(
      condition_label,
      levels = c("OCI-AML3 | ABZ", "OCI-AML3 | MBZ", "MV4-11 | ABZ", "MV4-11 | MBZ")
    ),
    direction = factor(direction, levels = c("core_up", "core_down")),
    gene = factor(gene, levels = rev(gene_order_heat))
  )

response_limit <- max(abs(response_heatmap_df$directional_response), na.rm = TRUE)

p_full_heat <- ggplot(
  response_heatmap_df,
  aes(x = condition_label, y = gene, fill = directional_response)
) +
  geom_tile(color = "white", linewidth = 0.16) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "grey98", high = "#B2182B",
    midpoint = 0, limits = c(-response_limit, response_limit),
    name = "Directional\nresponse"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Full-core directional response to drug treatment",
    subtitle = "Positive values indicate movement opposite to the predefined AML disease direction"
  ) +
  theme_sci(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5.5),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold")
  )

save_gg(p_full_heat, "FigS1_full_62gene_directional_response_heatmap", width = 6.5, height = 10.0)
save_gg(make_clean_plot(p_full_heat), "FigS1_full_62gene_directional_response_heatmap_clean", width = 6.5, height = 10.0)

############################################################
## 22. Supplementary Figure S2: cross-omics weight construction
############################################################

top_weight_labels <- core_weight %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(order_by = disease_weight, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup()

p_weight_map <- ggplot(core_weight, aes(x = B_rank01, y = S_rank01)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.42, color = "grey65") +
  geom_point(aes(size = disease_weight, color = direction), alpha = 0.82) +
  ggrepel::geom_text_repel(
    data = top_weight_labels,
    aes(label = gene),
    size = 2.6, show.legend = FALSE,
    box.padding = 0.3, point.padding = 0.2,
    min.segment.length = 0, max.overlaps = Inf, seed = 1234
  ) +
  scale_color_manual(values = direction_cols, name = "Core direction") +
  scale_size_area(max_size = 6.2, name = "Disease weight") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(
    x = "Rank-normalized |bulk RNA logFC|",
    y = "Rank-normalized |scRNA pseudo-bulk logFC|",
    title = "Cross-omics construction of disease-effect weights"
  ) +
  theme_sci(base_size = 9)

save_gg(p_weight_map, "FigS2_cross_omics_weight_construction", width = 5.6, height = 4.8)

############################################################
## 23. Supplementary Figure S3: top weighted contributions
############################################################

top_contrib <- drug_logfc %>%
  dplyr::group_by(cell_line, drug) %>%
  dplyr::slice_max(order_by = abs(weighted_contribution), n = 10, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    panel_gene = paste(cell_line, drug, gene, sep = "___")
  ) %>%
  dplyr::arrange(cell_line, drug, weighted_contribution) %>%
  dplyr::mutate(panel_gene = factor(panel_gene, levels = unique(panel_gene)))

p_contrib <- ggplot(
  top_contrib,
  aes(x = weighted_contribution, y = panel_gene, fill = direction)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.42, color = "grey55") +
  geom_col(width = 0.72) +
  facet_grid(cell_line ~ drug, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = direction_cols, name = "Core direction") +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  labs(
    x = "Weighted contribution to R_TOTAL",
    y = NULL,
    title = "Major gene-level contributors to weighted reversal"
  ) +
  theme_sci(base_size = 8.5) +
  theme(axis.text.y = element_text(size = 6.5))

save_gg(p_contrib, "FigS3_top_weighted_gene_contributions", width = 7.1, height = 6.2)

############################################################
## 24. Supplementary Figure S4: alignment permutation nulls
############################################################

perm_null_df <- purrr::imap_dfr(
  perm_distribution_list,
  function(scores_i, split_name_i) {
    dat_i <- split_data[[split_name_i]]
    tibble::tibble(
      cell_line = dat_i$cell_line[[1]],
      drug = dat_i$drug[[1]],
      perm_score = as.numeric(scores_i)
    )
  }
)

perm_observed_df <- score_table %>%
  dplyr::mutate(
    label = paste0(
      "R = ", sprintf("%.3f", R_TOTAL),
      "\nq_align = ", format.pval(q_align, digits = 2, eps = 1e-4)
    )
  )

p_perm_null <- ggplot(perm_null_df, aes(x = perm_score)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "grey82", color = "white", linewidth = 0.14) +
  geom_density(linewidth = 0.65, color = "grey25") +
  geom_vline(
    data = perm_observed_df,
    aes(xintercept = R_TOTAL),
    color = "#B2182B", linewidth = 0.85, inherit.aes = FALSE
  ) +
  geom_text(
    data = perm_observed_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.06, vjust = 1.20, size = 2.6,
    inherit.aes = FALSE
  ) +
  facet_grid(cell_line ~ drug, scales = "free") +
  labs(
    x = "Permuted alignment score",
    y = "Density",
    title = "Gene-label permutation null distributions"
  ) +
  theme_sci(base_size = 8.5)

save_gg(p_perm_null, "FigS4_gene_label_alignment_permutation_nulls", width = 6.8, height = 5.2)

############################################################
## 25. Export plot data and methods note
############################################################

readr::write_csv(weighted_scatter_df, file.path(TABLE_DIR, "10_main_scatter_plot_data.csv"))
readr::write_csv(scatter_label_df, file.path(TABLE_DIR, "11_main_scatter_labeled_genes.csv"))
readr::write_csv(architecture_df, file.path(TABLE_DIR, "12_architecture_heatmap_data.csv"))
readr::write_csv(response_heatmap_df, file.path(TABLE_DIR, "13_full_directional_response_heatmap_data.csv"))

methods_note <- c(
  "ANALYTICAL BOUNDARIES",
  "",
  "1. Input expression must be a non-negative integer raw-count matrix.",
  "   The script stops on negative or non-integer values and never rounds or truncates them.",
  "",
  "2. Each cell-line x drug comparison is modeled separately with DESeq2.",
  "   Same-batch design: ~ condition_model.",
  "   Pairing metadata, if present, is used for control matching but not as a DESeq2 batch covariate.",
  "",
  "3. Controls are not silently pooled.",
  "   By default an explicit comparison_id/control_group/experiment_batch/batch column is required.",
  "   Shared controls without such metadata require an explicit user declaration in the settings.",
  "",
  paste0("4. drug_induced_logFC is DESeq2-shrunken treated-vs-control log2FC using type='", LFC_SHRINK_TYPE, "'."),
  "   DESeq2 differential-expression p-values/padj are retained as gene-level model outputs but do not define the reversal statistic.",
  "",
  "5. Disease weights preserve the original formula:",
  "   B_i = non-zero rank01(|bulk logFC_i|); S_i = non-zero rank01(|scRNA logFC_i|);",
  "   M_i = sqrt(B_i * S_i); weights are normalized separately within core_up and core_down.",
  "",
  "6. R_UP = sum_i w_i * (-drug_LFC_i) among core_up genes.",
  "   R_DOWN = sum_i w_i * drug_LFC_i among core_down genes.",
  "   R_TOTAL = 0.5*R_UP + 0.5*R_DOWN.",
  "",
  "7. q_align is BH-adjusted one-sided gene-label permutation q-value.",
  "   It tests alignment of the 62-gene drug effects with predefined disease directions/weights.",
  "   It is not a drug-treatment differential-expression FDR."
)
writeLines(methods_note, file.path(OUT_DIR, "README_analysis_definitions.txt"))

############################################################
## 26. Session information
############################################################

sink(file.path(OUT_DIR, "sessionInfo.txt"))
sessionInfo()
sink()

############################################################
## 27. Final message
############################################################

message2("Analysis completed.")
message2("Output directory: ", OUT_DIR)
