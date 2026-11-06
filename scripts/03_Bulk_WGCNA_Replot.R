############################################################
## GSE114868 WGCNA replot only
## Purpose:
##   Re-draw WGCNA module-trait heatmap in the classic WGCNA
##   labeledHeatmap() style.
##
## Main changes compared with previous readable version:
##   1. Keep the original WGCNA heatmap style.
##   2. Do NOT add an extra left module-color block.
##   3. Do NOT display FDR in cells.
##   4. Display only correlation coefficient and P value:
##          r
##        (P value)
##   5. Increase PDF size, margins and tune cex.text to avoid label overlap.
##   6. Remove all fdr_mat dependency to avoid object-not-found errors.
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

## =========================
## 0. User settings
## =========================

OUTDIR <- BULK_WGCNA_RESULTS_DIR

TRAIT_DIR <- file.path(OUTDIR, "04_module_trait")
POWER_DIR <- file.path(OUTDIR, "02_soft_power")

PLOT_AML_ONLY <- TRUE
PLOT_ALL_TRAITS <- TRUE
PLOT_AML_BARPLOT <- TRUE
PLOT_SOFT_POWER <- TRUE

AML_ONLY_PDF_WIDTH  <- 6.5
AML_ONLY_PDF_HEIGHT <- 8.0
AML_ONLY_TEXT_CEX   <- 0.60

ALL_TRAITS_PDF_WIDTH  <- 7.2
ALL_TRAITS_PDF_HEIGHT <- 8.2
ALL_TRAITS_TEXT_CEX   <- 0.55

## =========================
## 1. Packages
## =========================

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("Package not installed: WGCNA")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package not installed: ggplot2")
}

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
})

## =========================
## 2. Helper functions
## =========================

read_mat <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  x <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    row.names = 1,
    quote = "",
    comment.char = ""
  )
  as.matrix(x)
}

## Classic WGCNA-like P value format.
## Keep it short to prevent overlapping.
format_p_wgcna <- function(p) {
  out <- rep(NA_character_, length(p))
  out[is.na(p)] <- "NA"

  idx <- !is.na(p)
  out[idx] <- ifelse(
    p[idx] < 1e-300,
    "<1e-300",
    format(signif(p[idx], 2), scientific = TRUE, trim = TRUE)
  )

  out
}

plot_module_trait_heatmap_wgcna_style <- function(
    cor_mat,
    p_mat,
    traits = colnames(cor_mat),
    file,
    title = "Module-trait relationships(GEO)",
    width = 7,
    height = 8,
    cex_text = 0.55,
    cex_lab = 0.85,
    bottom_margin = 7,
    left_margin = 11,
    top_margin = 3,
    right_margin = 3
) {
  traits <- intersect(traits, colnames(cor_mat))
  if (length(traits) == 0) {
    stop("No requested traits found in cor_mat.")
  }

  moduleTraitCor <- cor_mat[, traits, drop = FALSE]
  moduleTraitPvalue <- p_mat[, traits, drop = FALSE]

  if (!identical(dim(moduleTraitCor), dim(moduleTraitPvalue))) {
    stop("cor_mat and p_mat dimensions do not match after trait selection.")
  }

  textMatrix <- paste0(
    signif(moduleTraitCor, 2),
    "\n(",
    format_p_wgcna(moduleTraitPvalue),
    ")"
  )
  dim(textMatrix) <- dim(moduleTraitCor)

  pdf(file = file, width = width, height = height)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  par(mar = c(bottom_margin, left_margin, top_margin, right_margin))

  WGCNA::labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = traits,
    yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor),
    colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = cex_text,
    cex.lab = cex_lab,
    zlim = c(-1, 1),
    main = title
  )

  invisible(file)
}

## =========================
## 3. Read WGCNA module-trait results
## =========================

cor_mat <- read_mat(file.path(TRAIT_DIR, "01_module_trait_correlation.tsv"))
p_mat   <- read_mat(file.path(TRAIT_DIR, "02_module_trait_pvalue.tsv"))

if (!identical(rownames(cor_mat), rownames(p_mat))) {
  stop("Row names of correlation and p-value matrices are not identical.")
}
if (!identical(colnames(cor_mat), colnames(p_mat))) {
  stop("Column names of correlation and p-value matrices are not identical.")
}

message("Module-trait matrix loaded: ", nrow(cor_mat), " modules x ", ncol(cor_mat), " traits")
message("Traits found: ", paste(colnames(cor_mat), collapse = ", "))

## =========================
## 4. Replot module-trait heatmap
## =========================

if (PLOT_AML_ONLY && "AML_status" %in% colnames(cor_mat)) {
  plot_module_trait_heatmap_wgcna_style(
    cor_mat = cor_mat,
    p_mat = p_mat,
    traits = c("AML_status"),
    file = file.path(TRAIT_DIR, "04c_module_trait_heatmap_AML_only_WGCNAstyle_noFDR.pdf"),
    title = "Module-trait relationships(GEO)",
    width = AML_ONLY_PDF_WIDTH,
    height = AML_ONLY_PDF_HEIGHT,
    cex_text = AML_ONLY_TEXT_CEX,
    cex_lab = 0.88,
    bottom_margin = 7,
    left_margin = 11,
    top_margin = 3,
    right_margin = 3
  )
  message("Saved AML-only WGCNA-style heatmap.")
}

if (PLOT_ALL_TRAITS) {
  plot_module_trait_heatmap_wgcna_style(
    cor_mat = cor_mat,
    p_mat = p_mat,
    traits = colnames(cor_mat),
    file = file.path(TRAIT_DIR, "04d_module_trait_heatmap_all_traits_WGCNAstyle_noFDR.pdf"),
    title = "Module-trait relationships(GEO)",
    width = ALL_TRAITS_PDF_WIDTH,
    height = ALL_TRAITS_PDF_HEIGHT,
    cex_text = ALL_TRAITS_TEXT_CEX,
    cex_lab = 0.86,
    bottom_margin = 7,
    left_margin = 11,
    top_margin = 3,
    right_margin = 3
  )
  message("Saved all-traits WGCNA-style heatmap.")
}

## =========================
## 5. Optional: AML correlation barplot without FDR
## =========================

if (PLOT_AML_BARPLOT && "AML_status" %in% colnames(cor_mat)) {
  module_trait_long <- do.call(rbind, lapply(colnames(cor_mat), function(tr) {
    data.frame(
      module_eigengene = rownames(cor_mat),
      module = sub("^ME", "", rownames(cor_mat)),
      trait = tr,
      correlation = as.numeric(cor_mat[, tr]),
      pvalue = as.numeric(p_mat[, tr]),
      stringsAsFactors = FALSE
    )
  }))

  write.table(
    module_trait_long,
    file = file.path(TRAIT_DIR, "04e_module_trait_long_correlation_pvalue_noFDR.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  aml_cor_plot <- module_trait_long[
    module_trait_long$trait == "AML_status" & module_trait_long$module != "grey",
    ,
    drop = FALSE
  ]
  aml_cor_plot <- aml_cor_plot[order(aml_cor_plot$correlation), , drop = FALSE]
  aml_cor_plot$module <- factor(aml_cor_plot$module, levels = aml_cor_plot$module)

  p_cor_bar <- ggplot(aml_cor_plot, aes(x = module, y = correlation, fill = correlation)) +
    geom_col(width = 0.75, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Module association with AML status",
      x = "Module",
      y = "Correlation with AML status"
    )

  ggsave(
    filename = file.path(TRAIT_DIR, "05c_module_AML_correlation_barplot_noFDR.pdf"),
    plot = p_cor_bar,
    width = 6.2,
    height = max(4.5, 0.30 * nrow(aml_cor_plot))
  )

  message("Saved AML correlation barplot without FDR.")
}

## =========================
## 6. Optional: readable soft-power plots
## =========================

if (PLOT_SOFT_POWER) {
  sft_file <- file.path(POWER_DIR, "01_soft_threshold_fit_indices.tsv")
  power_file <- file.path(POWER_DIR, "02_selected_soft_power.txt")

  if (file.exists(sft_file) && file.exists(power_file)) {
    sft_tab <- read.table(sft_file, header = TRUE, sep = "\t", check.names = FALSE)
    softPower <- as.numeric(readLines(power_file, warn = FALSE)[1])

    if (!"Power" %in% colnames(sft_tab)) {
      stop("Soft-threshold table does not contain column: Power")
    }

    if (!"signed_R2" %in% colnames(sft_tab)) {
      if (all(c("slope", "SFT.R.sq") %in% colnames(sft_tab))) {
        sft_tab$signed_R2 <- -sign(sft_tab$slope) * sft_tab$SFT.R.sq
      } else {
        stop("Cannot find signed_R2, or slope + SFT.R.sq columns in soft-threshold table.")
      }
    }

    mean_k_col <- grep("^mean\\.k\\.?$|mean.k", colnames(sft_tab), value = TRUE)[1]
    if (is.na(mean_k_col)) {
      stop("Cannot find mean connectivity column in soft-threshold table.")
    }

    label_powers <- unique(c(1:10, 12, 16, 20, 24, 28, 30, softPower))
    label_idx <- sft_tab$Power %in% label_powers

    pdf(file.path(POWER_DIR, "03_soft_power_diagnostics_readable.pdf"), width = 11, height = 5.5)
    oldpar <- par(no.readonly = TRUE)
    on.exit({
      par(oldpar)
      dev.off()
    }, add = TRUE)

    par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

    plot(
      sft_tab$Power,
      sft_tab$signed_R2,
      xlab = "Soft threshold power",
      ylab = "Signed scale-free topology fit R^2",
      type = "b",
      pch = 19,
      cex = 0.75,
      main = "Scale independence"
    )
    text(
      sft_tab$Power[label_idx],
      sft_tab$signed_R2[label_idx],
      labels = sft_tab$Power[label_idx],
      cex = 0.75,
      pos = 3,
      offset = 0.35
    )
    abline(h = 0.85, col = "red", lty = 2)
    abline(v = softPower, col = "blue", lty = 2)
    legend("bottomright", legend = paste0("Selected power = ", softPower), bty = "n", cex = 0.85)

    plot(
      sft_tab$Power,
      sft_tab[[mean_k_col]],
      xlab = "Soft threshold power",
      ylab = "Mean connectivity",
      type = "b",
      pch = 19,
      cex = 0.75,
      main = "Mean connectivity"
    )
    text(
      sft_tab$Power[label_idx],
      sft_tab[[mean_k_col]][label_idx],
      labels = sft_tab$Power[label_idx],
      cex = 0.75,
      pos = 3,
      offset = 0.35
    )
    abline(v = softPower, col = "blue", lty = 2)
    legend("topright", legend = paste0("Selected power = ", softPower), bty = "n", cex = 0.85)

    message("Saved readable soft-power diagnostics.")
  } else {
    message("Soft-power input files not found. Skipped soft-power replot.")
  }
}

message("Done. Figures saved in: ", normalizePath(OUTDIR, mustWork = FALSE))
