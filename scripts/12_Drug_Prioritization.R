############################################################
## AML 62-gene perturbational drug screening
## Final revision based on work 5.6 critique
##
## Scientific boundaries enforced in code:
## 1) cMAP and L1000FWD are treated as two L1000-based
##    perturbational resources, NOT independent databases.
## 2) Albendazole / mebendazole are selected AFTER scoring as
##    experimentally feasible representatives of the mechanism axis.
## 3) cMAP and L1000FWD are scored separately first.
##    Cross-resource integration is calculated ONLY for same-drug
##    compounds present in both resources.
## 4) No clinical/AML relevance prior and no mechanism prior enters
##    the primary evidence scores.
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
set.seed(20260710)

############################################################
## 0. Packages
############################################################

pkgs <- c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "purrr",
  "ggplot2", "ggrepel", "forcats", "patchwork", "scales", "tibble"
)

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(pkgs, library, character.only = TRUE))

############################################################
## 1. Input / output
############################################################

CMAP_FILE <- DRUG_CMAP_FILE
L1000_FILE <- DRUG_L1000_FILE

OUT_DIR <- DRUG_PRIORITIZATION_DIR
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")

for (dd in c(OUT_DIR, FIG_DIR, TAB_DIR)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

############################################################
## 2. Analysis settings
############################################################

## cMAP export is a pre-filtered negative-hit table.
CMAP_INPUT_MODE <- "prefiltered_negative_hits"
CMAP_SIGNATURE_ID_COL <- "id"

CMAP_PERT_TYPE        <- "trt_cp"
CMAP_REQUIRE_QC_PASS  <- 1
CMAP_RAW_CS_CUT       <- 0
CMAP_NORM_CS_CUT      <- 0
CMAP_FDR_NLOG10_CUT   <- 1

## L1000FWD input must be a full/bidirectional similarity export
## if negative_fraction is used.
L1000_INPUT_SCOPE <- "full_or_bidirectional_similarity_export"
## Set TRUE only after confirming that L1000_FILE is the downloadable
## full similarity-result CSV, not merely the displayed top-N lists.
L1000_FULL_EXPORT_CONFIRMED_BY_USER <- TRUE
REQUIRE_L1000_BIDIRECTIONAL_SCORES <- TRUE
STOP_ON_L1000_DUPLICATE_SCORE_CONFLICT <- TRUE
L1000_NEG_SIM_CUT <- 0
MIN_L1000_SIGNATURE_N <- 3

## cMAP sensitivity score only
W_CMAP_SENS_BEST_NORM <- 0.45
W_CMAP_SENS_MEAN_NORM <- 0.25
W_CMAP_SENS_FDR       <- 0.15
W_CMAP_SENS_HIQ       <- 0.10
W_CMAP_SENS_TAS       <- 0.05

## Mechanism-axis display
MIN_AXIS_DRUGS_FOR_PLOT <- 3
TOP_N_AXES_FOR_PLOT <- 15

## Experimental validation compounds are defined post-scoring.
SELECTED_DRUGS <- c("albendazole", "mebendazole")

############################################################
## 3. Publication figure settings
############################################################

FIG_DPI <- 600
FIG_FONT <- "Helvetica"
USE_CAIRO_DEVICE <- FALSE
FIG_AXIS_TEXT <- 7.5
FIG_AXIS_TITLE <- 8.5
FIG_LEGEND_TEXT <- 7
FIG_LEGEND_TITLE <- 7.5

FIG_COLORS <- c(
  selected = "#B2182B",
  benzimidazole = "#EF8A62",
  microtubule = "#2166AC",
  neutral = "grey78",
  grid = "grey86",
  axis = "black"
)

.figure_export_log <- tibble::tibble(
  figure = character(),
  format = character(),
  width_in = numeric(),
  height_in = numeric(),
  dpi = numeric(),
  pixel_width = numeric(),
  pixel_height = numeric(),
  path = character()
)

############################################################
## 4. Curated structural-class table
##
## Do NOT infer class membership from a "bendazole" regex.
## This a-priori table is used only AFTER resource scoring.
## Add/remove compounds only with documented structural evidence.
############################################################

BENZIMIDAZOLE_CLASS_TABLE <- tibble::tribble(
  ~drug_std,       ~structural_class,                 ~curation_note,
  "albendazole",   "Benzimidazole compound class",   "a-priori curated",
  "mebendazole",   "Benzimidazole compound class",   "a-priori curated",
  "flubendazole",  "Benzimidazole compound class",   "a-priori curated",
  "parbendazole",  "Benzimidazole compound class",   "a-priori curated",
  "fenbendazole",  "Benzimidazole compound class",   "a-priori curated",
  "oxfendazole",   "Benzimidazole compound class",   "a-priori curated",
  "thiabendazole", "Benzimidazole compound class",   "a-priori curated",
  "ricobendazole", "Benzimidazole compound class",   "a-priori curated",
  "nocodazole",    "Benzimidazole compound class",   "a-priori curated"
)

readr::write_csv(
  BENZIMIDAZOLE_CLASS_TABLE,
  file.path(TAB_DIR, "00_curated_benzimidazole_class_table.csv")
)

############################################################
## 5. Helper functions
############################################################

clean_drug_name <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_to_lower() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("_", "-") %>%
    stringr::str_replace_all("\\s+", "-")
}

na666_chr <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x %in% c("", "-666", "NA", "NaN")] <- NA_character_
  x
}

to_num_clean <- function(x) {
  suppressWarnings(as.numeric(na666_chr(x)))
}

first_nonmissing_chr <- function(x) {
  x <- na666_chr(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

first_finite_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x[[1]]
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else median(x)
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

collapse_unique <- function(x, n = 12) {
  x <- na666_chr(x)
  x <- unique(x[!is.na(x)])
  if (length(x) == 0) return(NA_character_)
  paste(head(x, n), collapse = "; ")
}

rank01 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) == 0) return(out)
  xx <- x[ok]
  if (length(unique(xx)) == 1) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (rank(xx, ties.method = "average") - 1) / (length(xx) - 1)
  }
  out
}

weighted_available_mean <- function(df, cols, weights) {
  stopifnot(length(cols) == length(weights))
  mat <- as.matrix(df[, cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  ok <- is.finite(mat)
  mat0 <- mat
  mat0[!ok] <- 0
  w <- as.numeric(weights)
  num <- rowSums(sweep(mat0, 2, w, "*"), na.rm = TRUE)
  den <- rowSums(sweep(ok * 1, 2, w, "*"), na.rm = TRUE)
  ifelse(den > 0, num / den, NA_real_)
}

equal_available_mean <- function(df, cols) {
  mat <- as.matrix(df[, cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  ok <- is.finite(mat)
  mat0 <- mat
  mat0[!ok] <- 0
  num <- rowSums(mat0, na.rm = TRUE)
  den <- rowSums(ok)
  ifelse(den > 0, num / den, NA_real_)
}

primary_axis <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "Unannotated"]
  if (length(x) == 0) return("Unannotated")
  names(sort(table(x), decreasing = TRUE))[1]
}

annotate_cmap_mechanism <- function(moa_text, target_text) {
  moa_text <- ifelse(is.na(moa_text), "", as.character(moa_text))
  target_text <- ifelse(is.na(target_text), "", as.character(target_text))
  all_text <- paste(moa_text, target_text, sep = " | ") %>% stringr::str_to_lower()
  
  dplyr::case_when(
    stringr::str_detect(all_text, "tubulin|microtubule|\\btuba[0-9a-z]*\\b|\\btubb[0-9a-z]*\\b|spindle") ~ "Tubulin / microtubule axis",
    stringr::str_detect(all_text, "\\bhdac\\b|histone deacetylase") ~ "HDAC inhibitor",
    stringr::str_detect(all_text, "topoisomerase") ~ "Topoisomerase inhibitor",
    stringr::str_detect(all_text, "proteasome") ~ "Proteasome inhibitor",
    stringr::str_detect(all_text, "\\bcdk[0-9a-z]*\\b|cyclin-dependent kinase") ~ "CDK inhibitor",
    stringr::str_detect(all_text, "\\bmek\\b|\\bmap2k[0-9a-z]*\\b") ~ "MEK inhibitor",
    stringr::str_detect(all_text, "\\braf\\b|\\bbraf\\b|\\bcraf\\b") ~ "RAF inhibitor",
    stringr::str_detect(all_text, "\\begfr\\b") ~ "EGFR inhibitor",
    stringr::str_detect(all_text, "\\bjak[0-9a-z]*\\b") ~ "JAK inhibitor",
    stringr::str_detect(all_text, "\\bpi3k\\b|\\bpik3[0-9a-z]*\\b") ~ "PI3K inhibitor",
    stringr::str_detect(all_text, "\\bmtor\\b") ~ "mTOR inhibitor",
    stringr::str_detect(all_text, "\\bparp[0-9a-z]*\\b") ~ "PARP inhibitor",
    stringr::str_detect(all_text, "\\bhmgcr\\b") ~ "HMGCR inhibitor",
    stringr::str_detect(all_text, "retinoid|\\brar[abg]?\\b|\\brxr[abg]?\\b") ~ "Retinoid receptor axis",
    stringr::str_detect(all_text, "glucocorticoid") ~ "Glucocorticoid receptor axis",
    stringr::str_detect(all_text, "phosphodiesterase|\\bpde[0-9a-z]*\\b") ~ "Phosphodiesterase inhibitor",
    stringr::str_detect(all_text, "calcium channel") ~ "Calcium channel blocker",
    stringr::str_detect(all_text, "serotonin receptor|5-ht") ~ "Serotonin receptor axis",
    stringr::str_detect(all_text, "dopamine receptor") ~ "Dopamine receptor axis",
    stringr::str_detect(all_text, "angiotensin receptor") ~ "Angiotensin receptor axis",
    stringr::str_detect(all_text, "nf.?kb|\\bnfkb") ~ "NF-kB axis",
    stringr::str_detect(all_text, "na.?/?.?k.?atpase|sodium potassium atpase") ~ "Na/K-ATPase axis",
    moa_text != "" ~ "Other annotated mechanism",
    TRUE ~ "Unannotated"
  )
}

annotate_cmap_mechanisms_multilabel <- function(moa_text, target_text) {
  moa_text <- ifelse(is.na(moa_text), "", as.character(moa_text))
  target_text <- ifelse(is.na(target_text), "", as.character(target_text))
  all_text <- paste(moa_text, target_text, sep = " | ") %>% stringr::str_to_lower()
  axes <- character(0)
  add_axis <- function(condition, label) if (isTRUE(condition)) axes <<- c(axes, label)
  
  add_axis(stringr::str_detect(all_text, "tubulin|microtubule|\\btuba[0-9a-z]*\\b|\\btubb[0-9a-z]*\\b|spindle"), "Tubulin / microtubule axis")
  add_axis(stringr::str_detect(all_text, "\\bhdac\\b|histone deacetylase"), "HDAC inhibitor")
  add_axis(stringr::str_detect(all_text, "topoisomerase"), "Topoisomerase inhibitor")
  add_axis(stringr::str_detect(all_text, "proteasome"), "Proteasome inhibitor")
  add_axis(stringr::str_detect(all_text, "\\bcdk[0-9a-z]*\\b|cyclin-dependent kinase"), "CDK inhibitor")
  add_axis(stringr::str_detect(all_text, "\\bmek\\b|\\bmap2k[0-9a-z]*\\b"), "MEK inhibitor")
  add_axis(stringr::str_detect(all_text, "\\braf\\b|\\bbraf\\b|\\bcraf\\b"), "RAF inhibitor")
  add_axis(stringr::str_detect(all_text, "\\begfr\\b"), "EGFR inhibitor")
  add_axis(stringr::str_detect(all_text, "\\bjak[0-9a-z]*\\b"), "JAK inhibitor")
  add_axis(stringr::str_detect(all_text, "\\bpi3k\\b|\\bpik3[0-9a-z]*\\b"), "PI3K inhibitor")
  add_axis(stringr::str_detect(all_text, "\\bmtor\\b"), "mTOR inhibitor")
  add_axis(stringr::str_detect(all_text, "\\bparp[0-9a-z]*\\b"), "PARP inhibitor")
  add_axis(stringr::str_detect(all_text, "\\bhmgcr\\b"), "HMGCR inhibitor")
  add_axis(stringr::str_detect(all_text, "retinoid|\\brar[abg]?\\b|\\brxr[abg]?\\b"), "Retinoid receptor axis")
  add_axis(stringr::str_detect(all_text, "glucocorticoid"), "Glucocorticoid receptor axis")
  add_axis(stringr::str_detect(all_text, "phosphodiesterase|\\bpde[0-9a-z]*\\b"), "Phosphodiesterase inhibitor")
  add_axis(stringr::str_detect(all_text, "calcium channel"), "Calcium channel blocker")
  add_axis(stringr::str_detect(all_text, "serotonin receptor|5-ht"), "Serotonin receptor axis")
  add_axis(stringr::str_detect(all_text, "dopamine receptor"), "Dopamine receptor axis")
  add_axis(stringr::str_detect(all_text, "angiotensin receptor"), "Angiotensin receptor axis")
  add_axis(stringr::str_detect(all_text, "nf.?kb|\\bnfkb"), "NF-kB axis")
  add_axis(stringr::str_detect(all_text, "na.?/?.?k.?atpase|sodium potassium atpase"), "Na/K-ATPase axis")
  
  axes <- unique(axes)
  if (length(axes) == 0) axes <- if (moa_text != "") "Other annotated mechanism" else "Unannotated"
  axes
}

close_device_safely <- function(dev_id) {
  devs <- grDevices::dev.list()
  if (!is.null(devs) && dev_id %in% devs) try(grDevices::dev.off(which = dev_id), silent = TRUE)
}

save_plot_both <- function(plot, stem, width, height, dpi = FIG_DPI) {
  pdf_file <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(stem, ".png"))
  export_dpi <- dpi
  
  pdf_dev <- NA_integer_
  tryCatch({
    if (isTRUE(USE_CAIRO_DEVICE) && isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf(
        pdf_file, width = width, height = height,
        onefile = FALSE, family = FIG_FONT, fallback_resolution = dpi
      )
    } else {
      grDevices::pdf(
        pdf_file, width = width, height = height,
        onefile = FALSE, useDingbats = FALSE, family = FIG_FONT
      )
    }
    pdf_dev <- grDevices::dev.cur()
    print(plot)
    grDevices::dev.off(which = pdf_dev)
    pdf_dev <- NA_integer_
  }, error = function(e) {
    if (is.finite(pdf_dev)) close_device_safely(pdf_dev)
    stop("PDF export failed for '", stem, "': ", conditionMessage(e), call. = FALSE)
  })
  
  png_dev <- NA_integer_
  tryCatch({
    png_args <- list(
      filename = png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white"
    )
    if (isTRUE(USE_CAIRO_DEVICE) && isTRUE(capabilities("cairo"))) png_args$type <- "cairo-png"
    do.call(grDevices::png, png_args)
    png_dev <- grDevices::dev.cur()
    print(plot)
    grDevices::dev.off(which = png_dev)
    png_dev <- NA_integer_
  }, error = function(e) {
    if (is.finite(png_dev)) close_device_safely(png_dev)
    stop("PNG export failed for '", stem, "': ", conditionMessage(e), call. = FALSE)
  })
  .figure_export_log <<- dplyr::bind_rows(
    .figure_export_log,
    tibble::tibble(
      figure = stem,
      format = c("pdf", "png"),
      width_in = width,
      height_in = height,
      dpi = c(NA_real_, export_dpi),
      pixel_width = c(NA_real_, round(width * export_dpi)),
      pixel_height = c(NA_real_, round(height * export_dpi)),
      path = c(pdf_file, png_file)
    )
  )

  invisible(c(pdf = pdf_file, png = png_file))
}

## Clean submission version: remove titles/subtitles/captions but KEEP panel tags.
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

theme_sci <- function(base_size = 8.5, base_family = FIG_FONT) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = element_text(family = base_family, color = "black"),
      axis.title = element_text(face = "plain", color = "black", size = FIG_AXIS_TITLE),
      axis.text = element_text(color = "black", size = FIG_AXIS_TEXT),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.30, color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = 9.2),
      plot.subtitle = element_text(color = "grey30", hjust = 0, size = 7.2),
      legend.title = element_text(face = "bold", size = FIG_LEGEND_TITLE),
      legend.text = element_text(size = FIG_LEGEND_TEXT),
      legend.key = element_blank(),
      legend.spacing.x = grid::unit(0.06, "in"),
      legend.spacing.y = grid::unit(0.02, "in"),
      plot.margin = margin(5, 6, 5, 5)
    )
}

############################################################
## 6. Read inputs and check columns
############################################################

if (!file.exists(CMAP_FILE)) stop("Cannot find CMAP_FILE: ", CMAP_FILE)
if (!file.exists(L1000_FILE)) stop("Cannot find L1000_FILE: ", L1000_FILE)

cmap_raw <- readxl::read_excel(CMAP_FILE)
l1000_raw <- readr::read_csv(L1000_FILE, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))

required_cmap_cols <- c(
  CMAP_SIGNATURE_ID_COL, "pert_iname", "cell_iname", "pert_type", "qc_pass",
  "is_hiq", "moa", "target_name", "raw_cs", "fdr_q_nlog10", "norm_cs"
)
required_l1000_cols <- c("sig_id", "drug", "pert_id", "similarity_scores")

missing_cmap <- setdiff(required_cmap_cols, colnames(cmap_raw))
missing_l1000 <- setdiff(required_l1000_cols, colnames(l1000_raw))
if (length(missing_cmap) > 0) stop("Missing cMAP columns: ", paste(missing_cmap, collapse = ", "))
if (length(missing_l1000) > 0) stop("Missing L1000FWD columns: ", paste(missing_l1000, collapse = ", "))

message("cMAP raw dimensions: ", nrow(cmap_raw), " rows x ", ncol(cmap_raw), " columns")
message("L1000FWD raw dimensions: ", nrow(l1000_raw), " rows x ", ncol(l1000_raw), " columns")

############################################################
## 7. cMAP cleaning, -666 handling, and signature-ID deduplication
############################################################

cmap0_pre <- cmap_raw %>%
  dplyr::mutate(
    cmap_signature_id = na666_chr(.data[[CMAP_SIGNATURE_ID_COL]]),
    pert_iname = na666_chr(pert_iname),
    drug_std = clean_drug_name(pert_iname),
    cell_iname = na666_chr(cell_iname),
    pert_type = na666_chr(pert_type),
    raw_cs = to_num_clean(raw_cs),
    norm_cs = to_num_clean(norm_cs),
    fdr_q_nlog10 = to_num_clean(fdr_q_nlog10),
    qc_pass = to_num_clean(qc_pass),
    is_hiq = to_num_clean(is_hiq),
    tas = if ("tas" %in% colnames(.)) to_num_clean(tas) else NA_real_,
    moa_clean = na666_chr(moa),
    target_clean = na666_chr(target_name)
  )

if (any(is.na(cmap0_pre$cmap_signature_id))) {
  stop("cMAP contains missing/blank signature IDs after -666 cleaning; safe deduplication is impossible.")
}
if (any(is.na(cmap0_pre$drug_std) | cmap0_pre$drug_std == "")) {
  stop("cMAP contains missing standardized drug names.")
}

cmap_duplicate_signature_qc <- cmap0_pre %>%
  dplyr::group_by(cmap_signature_id) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_drugs = dplyr::n_distinct(drug_std, na.rm = TRUE),
    n_cells = dplyr::n_distinct(cell_iname, na.rm = TRUE),
    n_raw_cs = dplyr::n_distinct(raw_cs[is.finite(raw_cs)]),
    n_norm_cs = dplyr::n_distinct(norm_cs[is.finite(norm_cs)]),
    n_fdr = dplyr::n_distinct(fdr_q_nlog10[is.finite(fdr_q_nlog10)]),
    n_qc_pass = dplyr::n_distinct(qc_pass[is.finite(qc_pass)]),
    n_hiq = dplyr::n_distinct(is_hiq[is.finite(is_hiq)]),
    n_tas = dplyr::n_distinct(tas[is.finite(tas)]),
    n_pert_type = dplyr::n_distinct(pert_type, na.rm = TRUE),
    n_moa = dplyr::n_distinct(moa_clean, na.rm = TRUE),
    n_target = dplyr::n_distinct(target_clean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_rows > 1)

readr::write_csv(
  cmap_duplicate_signature_qc,
  file.path(TAB_DIR, "00_cMAP_duplicate_signature_QC.csv")
)

cmap_duplicate_conflicts <- cmap_duplicate_signature_qc %>%
  dplyr::filter(
    n_drugs > 1 | n_cells > 1 | n_raw_cs > 1 | n_norm_cs > 1 |
      n_fdr > 1 | n_qc_pass > 1 | n_hiq > 1 | n_tas > 1 |
      n_pert_type > 1 | n_moa > 1 | n_target > 1
  )

if (nrow(cmap_duplicate_conflicts) > 0) {
  readr::write_csv(
    cmap_duplicate_conflicts,
    file.path(TAB_DIR, "00_cMAP_duplicate_signature_CONFLICTS.csv")
  )
  stop(
    "Conflicting repeated cMAP signature IDs detected. No arbitrary row was retained. ",
    "See 00_cMAP_duplicate_signature_CONFLICTS.csv"
  )
}

## Repeated signature IDs are collapsed only after conflict checks.
## IMPORTANT: do not keep an arbitrary first row. When duplicate rows differ
## only by missingness, take the unique available value deterministically.
## Any genuinely conflicting non-missing values have already triggered stop().
cmap0 <- cmap0_pre %>%
  dplyr::group_by(cmap_signature_id) %>%
  dplyr::summarise(
    pert_iname = first_nonmissing_chr(pert_iname),
    drug_std = first_nonmissing_chr(drug_std),
    cell_iname = first_nonmissing_chr(cell_iname),
    pert_type = first_nonmissing_chr(pert_type),
    raw_cs = first_finite_num(raw_cs),
    norm_cs = first_finite_num(norm_cs),
    fdr_q_nlog10 = first_finite_num(fdr_q_nlog10),
    qc_pass = first_finite_num(qc_pass),
    is_hiq = first_finite_num(is_hiq),
    tas = first_finite_num(tas),
    moa_clean = first_nonmissing_chr(moa_clean),
    target_clean = first_nonmissing_chr(target_clean),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    passes_upstream_hit_definition =
      pert_type == CMAP_PERT_TYPE &
      qc_pass == CMAP_REQUIRE_QC_PASS &
      raw_cs < CMAP_RAW_CS_CUT &
      norm_cs < CMAP_NORM_CS_CUT &
      fdr_q_nlog10 >= CMAP_FDR_NLOG10_CUT,
    passes_upstream_hit_definition = dplyr::coalesce(passes_upstream_hit_definition, FALSE)
  )

cmap_qc <- tibble::tibble(
  metric = c(
    "Input rows",
    "Unique signature IDs after safe deduplication",
    "Exact duplicate rows collapsed",
    "Unique drugs",
    "Rows with cell_iname missing after -666 cleaning",
    "Rows with TAS missing after -666 cleaning",
    "Rows with HIQ missing after -666 cleaning",
    "Rows satisfying full upstream hit definition"
  ),
  value = c(
    nrow(cmap0_pre),
    nrow(cmap0),
    nrow(cmap0_pre) - nrow(cmap0),
    dplyr::n_distinct(cmap0$drug_std),
    sum(is.na(cmap0$cell_iname)),
    sum(is.na(cmap0$tas)),
    sum(is.na(cmap0$is_hiq)),
    sum(cmap0$passes_upstream_hit_definition)
  )
)
readr::write_csv(cmap_qc, file.path(TAB_DIR, "00_cMAP_input_QC.csv"))

if (CMAP_INPUT_MODE != "prefiltered_negative_hits") {
  stop("This script requires CMAP_INPUT_MODE = 'prefiltered_negative_hits'.")
}

cmap_hit <- cmap0 %>% dplyr::filter(passes_upstream_hit_definition)
if (nrow(cmap_hit) == 0) stop("No cMAP rows remain after declared upstream hit-definition check.")

############################################################
## 8. cMAP mechanism annotation: moa/target only
############################################################

cmap_hit <- cmap_hit %>%
  dplyr::mutate(
    cmap_mechanism_axis_row = annotate_cmap_mechanism(moa_clean, target_clean),
    cmap_mechanism_axes_row = purrr::map2(moa_clean, target_clean, annotate_cmap_mechanisms_multilabel)
  )

cmap_mechanism_drug_axis_long <- cmap_hit %>%
  dplyr::select(drug_std, cmap_mechanism_axes_row) %>%
  tidyr::unnest_longer(cmap_mechanism_axes_row, values_to = "cmap_mechanism_axis") %>%
  dplyr::distinct(drug_std, cmap_mechanism_axis)

cmap_mechanism_axes_all <- cmap_mechanism_drug_axis_long %>%
  dplyr::group_by(drug_std) %>%
  dplyr::summarise(
    cmap_mechanism_axes_all = collapse_unique(cmap_mechanism_axis, n = 20),
    .groups = "drop"
  )

cmap_mechanism_drug <- cmap_hit %>%
  dplyr::group_by(drug_std) %>%
  dplyr::summarise(
    cmap_mechanism_axis = primary_axis(cmap_mechanism_axis_row),
    moa_all = collapse_unique(moa_clean, n = 12),
    target_all = collapse_unique(target_clean, n = 12),
    .groups = "drop"
  ) %>%
  dplyr::left_join(cmap_mechanism_axes_all, by = "drug_std")

############################################################
## 9. cMAP drug-level score
############################################################

cmap_drug <- cmap_hit %>%
  dplyr::group_by(drug_std) %>%
  dplyr::summarise(
    drug_cmap_name = first_nonmissing_chr(pert_iname),
    n_cmap_hits = dplyr::n(),
    n_cmap_cells = dplyr::n_distinct(cell_iname, na.rm = TRUE),
    best_norm_cs = safe_min(norm_cs),
    mean_norm_cs = safe_mean(norm_cs),
    median_norm_cs = safe_median(norm_cs),
    best_raw_cs = safe_min(raw_cs),
    max_fdr_q_nlog10 = safe_max(fdr_q_nlog10),
    mean_fdr_q_nlog10 = safe_mean(fdr_q_nlog10),
    median_fdr_q_nlog10 = safe_median(fdr_q_nlog10),
    hiq_fraction = if (all(is.na(is_hiq))) NA_real_ else mean(is_hiq == 1, na.rm = TRUE),
    mean_tas = safe_mean(tas),
    .groups = "drop"
  ) %>%
  dplyr::left_join(cmap_mechanism_drug, by = "drug_std") %>%
  dplyr::mutate(
    cmap_rank_median_norm = rank01(-median_norm_cs),
    cmap_rank_median_fdr = rank01(median_fdr_q_nlog10),
    cmap_rank_hiq = rank01(hiq_fraction),
    cmap_rank_tas = rank01(mean_tas),
    cmap_rank_mean_norm = rank01(-mean_norm_cs),
    cmap_sens_rank_best_norm = rank01(-best_norm_cs),
    cmap_sens_rank_mean_norm = rank01(-mean_norm_cs),
    cmap_sens_rank_max_fdr = rank01(max_fdr_q_nlog10),
    cmap_sens_hiq = pmin(pmax(hiq_fraction, 0), 1),
    cmap_sens_rank_tas = rank01(mean_tas)
  )

cmap_drug$cMAP_hit_score <- equal_available_mean(
  cmap_drug,
  cols = c("cmap_rank_median_norm", "cmap_rank_median_fdr", "cmap_rank_hiq", "cmap_rank_tas")
)

cmap_primary_component_matrix <- as.matrix(
  cmap_drug[, c("cmap_rank_median_norm", "cmap_rank_median_fdr", "cmap_rank_hiq", "cmap_rank_tas")]
)
storage.mode(cmap_primary_component_matrix) <- "numeric"
cmap_drug$n_cmap_primary_components <- rowSums(is.finite(cmap_primary_component_matrix))

cmap_drug$cMAP_hit_score_weighted_sensitivity <- weighted_available_mean(
  cmap_drug,
  cols = c(
    "cmap_sens_rank_best_norm", "cmap_sens_rank_mean_norm",
    "cmap_sens_rank_max_fdr", "cmap_sens_hiq", "cmap_sens_rank_tas"
  ),
  weights = c(
    W_CMAP_SENS_BEST_NORM, W_CMAP_SENS_MEAN_NORM,
    W_CMAP_SENS_FDR, W_CMAP_SENS_HIQ, W_CMAP_SENS_TAS
  )
)

cmap_drug <- cmap_drug %>%
  dplyr::mutate(
    cMAP_rank = rank(-cMAP_hit_score, ties.method = "min", na.last = "keep"),
    cMAP_rank_weighted_sensitivity = rank(-cMAP_hit_score_weighted_sensitivity, ties.method = "min", na.last = "keep"),
    cMAP_rank_shift_sensitivity_minus_primary = cMAP_rank_weighted_sensitivity - cMAP_rank
  ) %>%
  dplyr::arrange(cMAP_rank)

readr::write_csv(cmap_drug, file.path(TAB_DIR, "01_cMAP_prefiltered_hit_drug_evidence.csv"))

############################################################
## 10. L1000FWD input-scope QC and safe signature deduplication
############################################################

l1000_pre <- l1000_raw %>%
  dplyr::mutate(
    sig_id = na666_chr(sig_id),
    drug = na666_chr(drug),
    pert_id = na666_chr(pert_id),
    drug_std = clean_drug_name(drug),
    similarity_scores = to_num_clean(similarity_scores)
  ) %>%
  dplyr::filter(
    !is.na(sig_id), sig_id != "",
    !is.na(drug_std), drug_std != "",
    is.finite(similarity_scores)
  )

if (nrow(l1000_pre) == 0) {
  stop("No finite/nonempty L1000FWD similarity rows remain after cleaning.")
}

if (!isTRUE(L1000_FULL_EXPORT_CONFIRMED_BY_USER)) {
  stop(
    "L1000 full-result export has not been explicitly confirmed. ",
    "Verify that L1000_FILE is the downloadable full similarity-result CSV ",
    "rather than a displayed top-N list, then set ",
    "L1000_FULL_EXPORT_CONFIRMED_BY_USER <- TRUE."
  )
}

l1000_score_sign_qc <- tibble::tibble(
  metric = c(
    "L1000 input scope declaration",
    "Full export explicitly confirmed by user",
    "Finite similarity rows",
    "Negative similarity rows",
    "Positive similarity rows",
    "Zero similarity rows",
    "Minimum similarity",
    "Maximum similarity"
  ),
  value = c(
    L1000_INPUT_SCOPE,
    as.character(L1000_FULL_EXPORT_CONFIRMED_BY_USER),
    as.character(nrow(l1000_pre)),
    as.character(sum(l1000_pre$similarity_scores < 0)),
    as.character(sum(l1000_pre$similarity_scores > 0)),
    as.character(sum(l1000_pre$similarity_scores == 0)),
    as.character(min(l1000_pre$similarity_scores)),
    as.character(max(l1000_pre$similarity_scores))
  )
)
readr::write_csv(l1000_score_sign_qc, file.path(TAB_DIR, "00_L1000FWD_score_sign_QC.csv"))

if (REQUIRE_L1000_BIDIRECTIONAL_SCORES) {
  if (!any(l1000_pre$similarity_scores < 0) || !any(l1000_pre$similarity_scores > 0)) {
    stop(
      "L1000FWD input does not contain both negative and positive similarity scores. ",
      "negative_fraction is not interpretable from a one-sided/top-hit-only export."
    )
  }
}

## A signature ID should not map to multiple standardized drugs.
l1000_sig_cross_drug <- l1000_pre %>%
  dplyr::group_by(sig_id) %>%
  dplyr::summarise(n_drugs = dplyr::n_distinct(drug_std), .groups = "drop") %>%
  dplyr::filter(n_drugs > 1)

readr::write_csv(l1000_sig_cross_drug, file.path(TAB_DIR, "00_L1000FWD_sig_id_cross_drug_QC.csv"))
if (nrow(l1000_sig_cross_drug) > 0) {
  stop("Some L1000FWD sig_id values map to multiple standardized drugs. See QC table.")
}

l1000_duplicate_signature_qc <- l1000_pre %>%
  dplyr::group_by(drug_std, sig_id) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_unique_scores = dplyr::n_distinct(similarity_scores),
    score_min = min(similarity_scores),
    score_max = max(similarity_scores),
    n_pert_ids = dplyr::n_distinct(pert_id, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_rows > 1)

readr::write_csv(
  l1000_duplicate_signature_qc,
  file.path(TAB_DIR, "00_L1000FWD_duplicate_signature_QC.csv")
)

l1000_duplicate_score_conflicts <- l1000_duplicate_signature_qc %>%
  dplyr::filter(n_unique_scores > 1 | n_pert_ids > 1)

if (nrow(l1000_duplicate_score_conflicts) > 0) {
  readr::write_csv(
    l1000_duplicate_score_conflicts,
    file.path(TAB_DIR, "00_L1000FWD_duplicate_signature_SCORE_CONFLICTS.csv")
  )
  if (STOP_ON_L1000_DUPLICATE_SCORE_CONFLICT) {
    stop(
      "Repeated L1000FWD drug-signature pairs contain conflicting similarity scores. ",
      "No arbitrary first row was retained. See SCORE_CONFLICTS.csv"
    )
  }
}

## Deterministic aggregation: identical duplicates remain identical;
## conflicting duplicates (only if explicitly allowed) use the predeclared median.
l1000 <- l1000_pre %>%
  dplyr::group_by(drug_std, sig_id) %>%
  dplyr::summarise(
    drug = first_nonmissing_chr(drug),
    pert_id = collapse_unique(pert_id, n = 20),
    similarity_scores = median(similarity_scores),
    duplicate_n_rows = dplyr::n(),
    .groups = "drop"
  )

l1000_drug_all <- l1000 %>%
  dplyr::group_by(drug_std) %>%
  dplyr::summarise(
    drug_l1000_name = first_nonmissing_chr(drug),
    n_l1000_signatures = dplyr::n_distinct(sig_id),
    n_l1000_negative = sum(similarity_scores < L1000_NEG_SIM_CUT),
    l1000_negative_fraction = mean(similarity_scores < L1000_NEG_SIM_CUT),
    min_similarity = min(similarity_scores),
    mean_similarity = mean(similarity_scores),
    median_similarity = median(similarity_scores),
    .groups = "drop"
  )

l1000_removed_lt3 <- l1000_drug_all %>%
  dplyr::filter(n_l1000_signatures < MIN_L1000_SIGNATURE_N) %>%
  dplyr::arrange(n_l1000_signatures, drug_std)
readr::write_csv(l1000_removed_lt3, file.path(TAB_DIR, "08_L1000FWD_removed_drugs_lt3_signatures.csv"))

l1000_drug <- l1000_drug_all %>%
  dplyr::filter(n_l1000_signatures >= MIN_L1000_SIGNATURE_N) %>%
  dplyr::mutate(
    l1000_score_mean = rank01(-mean_similarity),
    l1000_score_median = rank01(-median_similarity),
    ## Rank-transform the directional-consistency fraction so the
    ## primary score is genuinely rank-based across all components.
    l1000_score_neg_fraction = rank01(l1000_negative_fraction),
    L1000FWD_reversal_score =
      (l1000_score_mean + l1000_score_median + l1000_score_neg_fraction) / 3
  ) %>%
  dplyr::arrange(dplyr::desc(L1000FWD_reversal_score)) %>%
  dplyr::mutate(L1000FWD_rank = dplyr::row_number())

readr::write_csv(l1000_drug, file.path(TAB_DIR, "02_L1000FWD_drug_evidence_nge3.csv"))

l1000_qc <- tibble::tibble(
  metric = c(
    "Input rows after finite/nonempty filtering",
    "Unique drug-signature pairs after deterministic aggregation",
    "Duplicate rows collapsed",
    "Unique drugs before n>=3 filter",
    "Drugs retained with n>=3 unique signatures",
    "Drugs excluded with n<3 unique signatures"
  ),
  value = c(
    nrow(l1000_pre),
    nrow(l1000),
    nrow(l1000_pre) - nrow(l1000),
    nrow(l1000_drug_all),
    nrow(l1000_drug),
    nrow(l1000_removed_lt3)
  )
)
readr::write_csv(l1000_qc, file.path(TAB_DIR, "00_L1000FWD_input_QC.csv"))

############################################################
## 11. Resource-specific evidence and cross-resource concordance
##
## Critical change:
## - NO single overall rank mixes cMAP-only, L1000FWD-only, and shared drugs.
## - cMAP and L1000FWD ranks remain resource-specific.
## - CrossResourceConcordanceScore exists ONLY for same-drug compounds
##   present in both resources.
############################################################

integrated <- dplyr::full_join(cmap_drug, l1000_drug, by = "drug_std") %>%
  dplyr::mutate(
    drug_name = dplyr::coalesce(drug_cmap_name, drug_l1000_name, drug_std),
    has_cMAP = is.finite(cMAP_hit_score),
    has_L1000FWD = is.finite(L1000FWD_reversal_score),
    same_drug_cross_resource = has_cMAP & has_L1000FWD,
    evidence_scope = dplyr::case_when(
      same_drug_cross_resource ~ "Both L1000-based resources",
      has_cMAP ~ "cMAP only",
      has_L1000FWD ~ "L1000FWD only",
      TRUE ~ "No retained evidence"
    ),
    CrossResourceConcordanceScore = ifelse(
      same_drug_cross_resource,
      0.5 * cMAP_hit_score + 0.5 * L1000FWD_reversal_score,
      NA_real_
    ),
    CrossResourceBalanceScore = ifelse(
      same_drug_cross_resource,
      sqrt(cMAP_hit_score * L1000FWD_reversal_score),
      NA_real_
    ),
    CrossResourceRank = ifelse(
      same_drug_cross_resource,
      rank(-CrossResourceConcordanceScore, ties.method = "min", na.last = "keep"),
      NA_real_
    )
  ) %>%
  dplyr::left_join(
    BENZIMIDAZOLE_CLASS_TABLE %>% dplyr::select(drug_std, structural_class, curation_note),
    by = "drug_std"
  ) %>%
  dplyr::mutate(
    structural_class = dplyr::coalesce(structural_class, "Other / unclassified")
  ) %>%
  dplyr::arrange(
    dplyr::desc(same_drug_cross_resource),
    CrossResourceRank,
    cMAP_rank,
    L1000FWD_rank
  )

readr::write_csv(integrated, file.path(TAB_DIR, "03_resource_specific_and_cross_resource_evidence.csv"))

shared_drugs <- integrated %>%
  dplyr::filter(same_drug_cross_resource) %>%
  dplyr::arrange(CrossResourceRank)

readr::write_csv(shared_drugs, file.path(TAB_DIR, "04_same_drug_cross_resource_concordance.csv"))

############################################################
## 12. cMAP-defined mechanism-axis summaries
############################################################

cmap_axis_drug <- cmap_mechanism_drug_axis_long %>%
  dplyr::left_join(
    cmap_drug %>% dplyr::select(drug_std, drug_cmap_name, n_cmap_hits, cMAP_hit_score),
    by = "drug_std"
  ) %>%
  dplyr::left_join(
    l1000_drug %>% dplyr::select(
      drug_std, L1000FWD_reversal_score, n_l1000_signatures,
      l1000_negative_fraction, median_similarity
    ),
    by = "drug_std"
  )

axis_summary <- cmap_axis_drug %>%
  dplyr::group_by(cmap_mechanism_axis) %>%
  dplyr::summarise(
    n_cmap_drugs = dplyr::n_distinct(drug_std),
    n_cmap_hits = sum(n_cmap_hits, na.rm = TRUE),
    median_cmap_hit_score = safe_median(cMAP_hit_score),
    mean_cmap_hit_score = safe_mean(cMAP_hit_score),
    n_same_drug_L1000 = sum(is.finite(L1000FWD_reversal_score)),
    same_drug_L1000_fraction = n_same_drug_L1000 / n_cmap_drugs,
    median_L1000_score_same_drug = ifelse(
      n_same_drug_L1000 > 0,
      safe_median(L1000FWD_reversal_score),
      NA_real_
    ),
    n_same_drug_negative_median = sum(is.finite(median_similarity) & median_similarity < 0),
    representative_drugs = paste(head(unique(drug_cmap_name), 15), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(median_cmap_hit_score))

readr::write_csv(axis_summary, file.path(TAB_DIR, "05_cMAP_defined_mechanism_axis_summary.csv"))

############################################################
## 13. Benzimidazole-class summary after scoring
############################################################

benzo_l1000_summary <- l1000_drug %>%
  dplyr::inner_join(
    BENZIMIDAZOLE_CLASS_TABLE %>% dplyr::select(drug_std, structural_class, curation_note),
    by = "drug_std"
  ) %>%
  dplyr::mutate(selected_for_validation = drug_std %in% SELECTED_DRUGS) %>%
  dplyr::arrange(median_similarity)

readr::write_csv(
  benzo_l1000_summary,
  file.path(TAB_DIR, "06_benzimidazole_L1000_summary_nge3.csv")
)

benzo_l1000_raw <- l1000 %>%
  dplyr::semi_join(benzo_l1000_summary %>% dplyr::select(drug_std), by = "drug_std") %>%
  dplyr::left_join(
    benzo_l1000_summary %>% dplyr::select(
      drug_std, drug_l1000_name, n_l1000_signatures,
      l1000_negative_fraction, median_similarity, selected_for_validation
    ),
    by = "drug_std"
  )

############################################################
## 14. Selected compounds: post-scoring extraction only
############################################################

selected_validation_drugs <- integrated %>%
  dplyr::filter(drug_std %in% SELECTED_DRUGS) %>%
  dplyr::mutate(
    validation_status = "Selected post-scoring as experimentally feasible mechanism representative",
    evidence_boundary = dplyr::case_when(
      same_drug_cross_resource ~ "Retained same-drug evidence in both L1000-based resources",
      has_cMAP ~ "Retained direct evidence in cMAP only",
      has_L1000FWD ~ "Retained direct evidence in eligible L1000FWD only",
      TRUE ~ "No retained direct evidence in scored tables"
    ),
    selection_statement = paste0(
      "Perturbational screening highlighted microtubule/benzimidazole-related signals; ",
      "compound selected subsequently for experimental feasibility as a mechanism representative"
    )
  ) %>%
  dplyr::arrange(match(drug_std, SELECTED_DRUGS))

readr::write_csv(
  selected_validation_drugs,
  file.path(TAB_DIR, "07_selected_validation_drugs_post_scoring.csv")
)

############################################################
## 15. Sensitivity analysis
############################################################

cmap_score_sensitivity <- cmap_drug %>%
  dplyr::transmute(
    drug_std, drug_cmap_name, n_cmap_hits, n_cmap_cells,
    cMAP_hit_score_primary = cMAP_hit_score,
    cMAP_rank_primary = cMAP_rank,
    cMAP_hit_score_weighted_sensitivity,
    cMAP_rank_weighted_sensitivity,
    rank_shift = cMAP_rank_weighted_sensitivity - cMAP_rank
  ) %>%
  dplyr::arrange(cMAP_rank_primary)

readr::write_csv(cmap_score_sensitivity, file.path(TAB_DIR, "12_cMAP_primary_vs_weighted_sensitivity.csv"))

rho_cmap_primary_vs_sensitivity <- suppressWarnings(
  stats::cor(
    cmap_drug$cMAP_hit_score,
    cmap_drug$cMAP_hit_score_weighted_sensitivity,
    method = "spearman",
    use = "complete.obs"
  )
)

shared_sensitivity <- shared_drugs %>%
  dplyr::mutate(
    CrossResourceConcordanceScore_weighted_cMAP_sensitivity =
      0.5 * cMAP_hit_score_weighted_sensitivity + 0.5 * L1000FWD_reversal_score
  )

rho_shared_primary_vs_sensitivity <- suppressWarnings(
  stats::cor(
    shared_sensitivity$CrossResourceConcordanceScore,
    shared_sensitivity$CrossResourceConcordanceScore_weighted_cMAP_sensitivity,
    method = "spearman",
    use = "complete.obs"
  )
)

readr::write_csv(
  shared_sensitivity,
  file.path(TAB_DIR, "13_cross_resource_primary_vs_cMAP_sensitivity.csv")
)

############################################################
## 16. README / reproducibility note
############################################################

score_note <- c(
  "SCIENTIFIC BOUNDARIES AND SCORE DEFINITIONS",
  "",
  "1. cMAP and L1000FWD are described as two L1000-based perturbational resources.",
  "   They are not treated as independent data sources or independent cross-database validation.",
  "",
  "2. cMAP input is a prefiltered negative-hit table.",
  "   Signature IDs are deduplicated only after conflict checks.",
  "   -666 sentinels in numeric and categorical fields are converted to NA before summaries.",
  "",
  "3. Primary cMAP score is the equal available-component mean of:",
  "   rank01(-median_norm_cs), rank01(median_fdr_q_nlog10),",
  "   rank01(hiq_fraction), and rank01(mean_tas).",
  "",
  "4. L1000FWD negative_fraction is used only after a bidirectional-score QC.",
  paste0("   Drugs require n >= ", MIN_L1000_SIGNATURE_N, " unique signatures."),
  "   Repeated drug-signature pairs are never resolved by arbitrary first-row retention.",
  "   Conflicting duplicate scores stop the analysis by default.",
  "   Primary L1000FWD score = mean of rank01(-mean_similarity),",
  "   rank01(-median_similarity), and rank01(negative_fraction).",
  "",
  "5. cMAP-only, L1000FWD-only, and shared compounds are not placed into one overall rank.",
  "   Resource-specific ranks are retained separately.",
  "   CrossResourceConcordanceScore is calculated only for same-drug compounds present in both resources.",
  "",
  "6. Mechanism and structural-class annotations are post-scoring and descriptive.",
  "   Benzimidazole membership uses an a-priori curated table, not a drug-name regex.",
  "",
  "7. Albendazole and mebendazole are defined only after scoring.",
  "   The defensible statement is that perturbational screening highlighted microtubule/benzimidazole-related signals,",
  "   followed by selection of albendazole and mebendazole as experimentally feasible mechanism representatives.",
  "   Do not state that objective ranking directly selected both compounds or that both were dual-resource high-ranking hits."
)
writeLines(score_note, file.path(OUT_DIR, "README_scientific_boundaries_and_scores.txt"))

############################################################
## 17. Figure 3A: four-step workflow
############################################################

node_df <- tibble::tibble(
  x = c(1, 2, 3, 4),
  y = 1,
  step = c("1", "2", "3", "4"),
  title = c(
    "AML core",
    "Perturbational scoring",
    "Mechanism annotation",
    "Validation candidates"
  ),
  subtitle = c(
    "62-gene signature",
    "cMAP and L1000FWD",
    "Microtubule / benzimidazole axis",
    "Albendazole and mebendazole"
  ),
  fill = c(
    "#E6F0F5",
    "#ECE7F4",
    "#F5EBDD",
    "#F5DDDA"
  )
)

edge_df <- tibble::tibble(
  x = c(1.32, 2.32, 3.32),
  xend = c(1.68, 2.68, 3.68),
  y = 1,
  yend = 1
)

pA <- ggplot() +
  geom_segment(
    data = edge_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.45,
    color = "grey35",
    arrow = grid::arrow(
      length = grid::unit(0.12, "cm"),
      type = "closed"
    )
  ) +
  geom_label(
    data = node_df,
    aes(
      x = x,
      y = y,
      label = paste0(
        step,
        "\n",
        title,
        "\n",
        subtitle
      ),
      fill = fill
    ),
    size = 2.55,
    fontface = "bold",
    lineheight = 0.92,
    label.size = 0.18,
    label.padding = grid::unit(0.20, "lines"),
    label.r = grid::unit(0.08, "lines"),
    color = "black"
  ) +
  scale_fill_identity() +
  coord_cartesian(
    xlim = c(0.55, 4.45),
    ylim = c(0.72, 1.28),
    clip = "off"
  ) +
  theme_void(
    base_family = FIG_FONT
  ) +
  labs(
    title = "Perturbational screening and post-scoring candidate selection"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 9.2,
      hjust = 0
    ),
    plot.margin = margin(
      2,
      4,
      2,
      4
    )
  )

save_plot_both(
  pA,
  "Fig3A_four_step_workflow",
  width = 7.1,
  height = 1.25
)

save_plot_both(
  make_clean_plot(pA),
  "Fig3A_four_step_workflow_clean",
  width = 7.1,
  height = 1.25
)

############################################################
## 18. Supplementary Figure: same-drug cross-resource scatter
## Kept as supplementary, not used as main Figure 3 panel.
############################################################

plot_shared <- shared_drugs %>%
  dplyr::mutate(
    plot_group = dplyr::case_when(
      drug_std %in% SELECTED_DRUGS ~ "Selected validation compound",
      structural_class == "Benzimidazole compound class" ~
        "Benzimidazole compound",
      !is.na(cmap_mechanism_axes_all) &
        stringr::str_detect(
          cmap_mechanism_axes_all,
          stringr::fixed("Tubulin / microtubule axis")
        ) ~ "Tubulin / microtubule compound",
      TRUE ~ "Other shared compounds"
    )
  )

FIG3B_LABEL_DRUGS <- c(
  "nch-51",
  "albendazole",
  "nocodazole",
  "cyt-997"
)

label_shared <- plot_shared %>%
  dplyr::filter(
    drug_std %in% FIG3B_LABEL_DRUGS
  ) %>%
  dplyr::distinct(
    drug_std,
    .keep_all = TRUE
  )

rho_shared <- suppressWarnings(
  stats::cor(
    plot_shared$cMAP_hit_score,
    plot_shared$L1000FWD_reversal_score,
    method = "spearman",
    use = "complete.obs"
  )
)

mebendazole_note <- if (!"mebendazole" %in% plot_shared$drug_std) {
  "Mebendazole was not retained as a same-drug cMAP hit."
} else {
  NULL
}

pB <- ggplot(
  plot_shared,
  aes(
    x = cMAP_hit_score,
    y = L1000FWD_reversal_score
  )
) +
  geom_point(
    data = subset(
      plot_shared,
      plot_group == "Other shared compounds"
    ),
    color = "grey80",
    size = 1.30,
    alpha = 0.62
  ) +
  geom_point(
    data = subset(
      plot_shared,
      plot_group != "Other shared compounds"
    ),
    aes(
      color = plot_group
    ),
    size = 2.55,
    alpha = 0.95
  ) +
  geom_point(
    data = subset(
      plot_shared,
      drug_std %in% SELECTED_DRUGS
    ),
    color = "#B2182B",
    size = 3.25,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = label_shared,
    aes(
      label = drug_name
    ),
    size = 2.45,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.22,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.25,
    max.overlaps = Inf,
    seed = 20260710
  ) +
  scale_color_manual(
    breaks = c(
      "Benzimidazole compound",
      "Tubulin / microtubule compound",
      "Selected validation compound"
    ),
    values = c(
      "Benzimidazole compound" = FIG_COLORS[["benzimidazole"]],
      "Tubulin / microtubule compound" = FIG_COLORS[["microtubule"]],
      "Selected validation compound" = FIG_COLORS[["selected"]]
    ),
    labels = c(
      "Benzimidazole",
      "Tubulin / microtubule",
      "Selected"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Same-drug cross-resource concordance",
    subtitle = paste0(
      "Two L1000-based perturbational resources; n = ",
      nrow(plot_shared),
      "; Spearman rho = ",
      sprintf("%.2f", rho_shared)
    ),
    x = "cMAP qualified-hit evidence score",
    y = "L1000FWD reversal evidence score",
    color = NULL,
    caption = mebendazole_note
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        size = 2.2,
        alpha = 1
      )
    )
  ) +
  theme_sci(
    base_size = 8.5
  ) +
  theme(
    legend.position = "right",
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.key.height = grid::unit(0.13, "in"),
    legend.key.width = grid::unit(0.13, "in"),
    legend.text = element_text(size = 6.2),
    plot.caption = element_text(
      size = 6.4,
      color = "grey35",
      hjust = 0
    )
  )

save_plot_both(
  pB,
  "FigS_cross_resource_concordance",
  width = 3.6,
  height = 3.4
)

save_plot_both(
  make_clean_plot(pB),
  "FigS_cross_resource_concordance_clean",
  width = 3.6,
  height = 3.4
)

############################################################
## 19. Main Figure B:
## Candidate evidence and post-scoring annotation matrix
############################################################

benzo_support_order <- benzo_l1000_summary %>%
  dplyr::filter(
    !drug_std %in% SELECTED_DRUGS
  ) %>%
  dplyr::arrange(
    median_similarity
  ) %>%
  dplyr::pull(
    drug_std
  )

candidate_matrix_order <- unique(
  c(
    "albendazole",
    "mebendazole",
    "nocodazole",
    "flubendazole",
    "parbendazole",
    "oxfendazole",
    "fenbendazole",
    "cyt-997"
  )
)

candidate_matrix_base <- integrated %>%
  dplyr::filter(
    drug_std %in% candidate_matrix_order
  ) %>%
  dplyr::mutate(
    selected_for_validation = drug_std %in% SELECTED_DRUGS,
    benzimidazole_class =
      structural_class == "Benzimidazole compound class",
    cmap_tubulin_annotation =
      !is.na(cmap_mechanism_axes_all) &
      stringr::str_detect(
        cmap_mechanism_axes_all,
        stringr::fixed("Tubulin / microtubule axis")
      )
  )

candidate_matrix_order_present <- candidate_matrix_order[
  candidate_matrix_order %in% candidate_matrix_base$drug_std
]

candidate_name_map <- candidate_matrix_base %>%
  dplyr::select(
    drug_std,
    drug_name
  ) %>%
  dplyr::distinct()

candidate_name_order <- candidate_name_map$drug_name[
  match(
    candidate_matrix_order_present,
    candidate_name_map$drug_std
  )
]

candidate_score_long <- candidate_matrix_base %>%
  dplyr::select(
    drug_std,
    drug_name,
    cMAP_hit_score,
    L1000FWD_reversal_score
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      cMAP_hit_score,
      L1000FWD_reversal_score
    ),
    names_to = "evidence_source",
    values_to = "evidence_score"
  ) %>%
  dplyr::mutate(
    evidence_source = dplyr::recode(
      evidence_source,
      cMAP_hit_score = "cMAP",
      L1000FWD_reversal_score = "L1000FWD"
    ),
    evidence_source = factor(
      evidence_source,
      levels = c(
        "cMAP",
        "L1000FWD"
      )
    ),
    drug_name = factor(
      drug_name,
      levels = rev(candidate_name_order)
    ),
    score_label = ifelse(
      is.finite(evidence_score),
      sprintf("%.2f", evidence_score),
      "\u2013"
    )
  )

candidate_annotation_long <- candidate_matrix_base %>%
  dplyr::select(
    drug_std,
    drug_name,
    cmap_tubulin_annotation,
    benzimidazole_class,
    selected_for_validation
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      cmap_tubulin_annotation,
      benzimidazole_class,
      selected_for_validation
    ),
    names_to = "annotation_type",
    values_to = "annotation_present"
  ) %>%
  dplyr::mutate(
    annotation_type = dplyr::recode(
      annotation_type,
      cmap_tubulin_annotation = "cMAP tubulin/MT",
      benzimidazole_class = "Benzimidazole class",
      selected_for_validation = "Selected"
    ),
    annotation_type = factor(
      annotation_type,
      levels = c(
        "cMAP tubulin/MT",
        "Benzimidazole class",
        "Selected"
      )
    ),
    drug_name = factor(
      drug_name,
      levels = rev(candidate_name_order)
    )
  )

pS_candidate_scores <- ggplot(
  candidate_score_long,
  aes(
    x = evidence_source,
    y = drug_name,
    fill = evidence_score
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.65
  ) +
  geom_text(
    aes(
      label = score_label
    ),
    size = 2.65,
    color = "black",
    fontface = "bold"
  ) +
  scale_fill_gradientn(
    colors = c(
      "#F7FBFF",
      "#C6DBEF",
      "#6BAED6",
      "#2171B5",
      "#08306B"
    ),
    limits = c(0, 1),
    na.value = "grey92",
    name = "Evidence\nscore"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_sci(
    base_size = 8.5
  ) +
  theme(
    axis.text.x = element_text(
      face = "bold",
      size = 7.2
    ),
    axis.text.y = element_text(
      size = 7.2,
      color = "black"
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 6.6,
      face = "bold"
    ),
    legend.text = element_text(
      size = 6.0
    )
  )

pS_candidate_annotations <- ggplot(
  candidate_annotation_long,
  aes(
    x = annotation_type,
    y = drug_name
  )
) +
  geom_tile(
    fill = "white",
    color = "grey85",
    linewidth = 0.55
  ) +
  geom_point(
    data = subset(
      candidate_annotation_long,
      !annotation_present
    ),
    shape = 21,
    size = 3.1,
    fill = "white",
    color = "grey60",
    stroke = 0.65
  ) +
  geom_point(
    data = subset(
      candidate_annotation_long,
      annotation_present
    ),
    aes(
      fill = annotation_type
    ),
    shape = 21,
    size = 3.45,
    color = "grey20",
    stroke = 0.65
  ) +
  scale_fill_manual(
    values = c(
      "cMAP tubulin/MT" = "#2166AC",
      "Benzimidazole class" = "#EF8A62",
      "Selected" = "#B2182B"
    ),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_sci(
    base_size = 8.5
  ) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 6.8
    ),
    axis.text.y = element_text(
      size = 7.2,
      color = "black"
    ),
    axis.ticks.y = element_blank()
  )

pS_candidate_annotations_for_combined <- pS_candidate_annotations +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

FigS_candidate_matrix <- patchwork::wrap_plots(
  pS_candidate_scores,
  pS_candidate_annotations_for_combined,
  ncol = 2,
  widths = c(0.78, 1.15)
) +
  patchwork::plot_annotation(
    title = "Candidate evidence and post-scoring annotations",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 9.2,
        hjust = 0
      )
    )
  )

save_plot_both(
  pS_candidate_scores,
  "FigS_candidate_evidence_matrix_left",
  width = 3.4,
  height = 4.5
)

save_plot_both(
  make_clean_plot(pS_candidate_scores),
  "FigS_candidate_evidence_matrix_left_clean",
  width = 3.4,
  height = 4.5
)

save_plot_both(
  pS_candidate_annotations,
  "FigS_candidate_annotation_matrix_right",
  width = 4.8,
  height = 4.5
)

save_plot_both(
  make_clean_plot(pS_candidate_annotations),
  "FigS_candidate_annotation_matrix_right_clean",
  width = 4.8,
  height = 4.5
)

save_plot_both(
  FigS_candidate_matrix,
  "FigS_candidate_evidence_annotation_matrix",
  width = 7.1,
  height = 3.35
)

save_plot_both(
  make_clean_plot(FigS_candidate_matrix),
  "FigS_candidate_evidence_annotation_matrix_clean",
  width = 7.1,
  height = 3.35
)

############################################################
## 19B. Main Figure B:
## Perturbational reversal evidence landscapes
##
## IMPORTANT SCIENTIFIC BOUNDARY:
## - Primary composite scores remain the ranking statistics.
## - These panels are descriptive visualizations of major
##   underlying evidence dimensions.
## - cMAP input contains qualified negative-hit signatures only.
## - Mechanism/class annotations are post-scoring and descriptive.
############################################################

CMAP_LANDSCAPE_LABEL_N_PER_GROUP  <- 3
L1000_LANDSCAPE_LABEL_N_PER_GROUP <- 3


############################################################
## Helper: descriptive post-scoring annotation
############################################################

annotate_landscape_group <- function(df) {
  
  df %>%
    dplyr::mutate(
      
      is_selected =
        drug_std %in% SELECTED_DRUGS,
      
      is_tubulin =
        dplyr::coalesce(
          !is.na(cmap_mechanism_axes_all) &
            stringr::str_detect(
              cmap_mechanism_axes_all,
              stringr::fixed(
                "Tubulin / microtubule axis"
              )
            ),
          FALSE
        ),
      
      is_benzimidazole =
        structural_class ==
        "Benzimidazole compound class",
      
      display_group =
        dplyr::case_when(
          
          is_tubulin ~
            "Tubulin / microtubule",
          
          is_benzimidazole ~
            "Benzimidazole class",
          
          TRUE ~
            "Other qualified compounds"
        ),
      
      display_group = factor(
        display_group,
        levels = c(
          "Other qualified compounds",
          "Tubulin / microtubule",
          "Benzimidazole class"
        )
      )
    )
}


############################################################
## Helper: data-driven labels
##
## Always retain selected validation compounds.
## Additional labels are the strongest reversal compounds
## among highlighted mechanism / structural categories.
############################################################

make_landscape_labels <- function(
    df,
    n_per_group = 3
) {
  
  dplyr::bind_rows(
    
    ## Always label selected validation compounds
    df %>%
      dplyr::filter(
        is_selected
      ),
    
    ## Label strongest compounds within each highlighted category
    df %>%
      dplyr::filter(
        !is_selected,
        display_group !=
          "Other qualified compounds"
      ) %>%
      dplyr::group_by(
        display_group
      ) %>%
      dplyr::arrange(
        dplyr::desc(
          reversal_magnitude
        ),
        .by_group = TRUE
      ) %>%
      dplyr::slice_head(
        n = n_per_group
      ) %>%
      dplyr::ungroup()
    
  ) %>%
    
    dplyr::distinct(
      drug_std,
      .keep_all = TRUE
    )
}

############################################################
## Shared display colors
############################################################

LANDSCAPE_FILL_COLORS <- c(
  
  "Tubulin / microtubule" =
    FIG_COLORS[["microtubule"]],
  
  "Benzimidazole class" =
    FIG_COLORS[["benzimidazole"]]
)


SELECTED_RING_COLOR <-
  FIG_COLORS[["selected"]]


############################################################
## cMAP qualified negative-hit landscape
##
## x:
## stronger negative median normCS -> larger reversal magnitude
##
## y:
## median fdr_q_nlog10 among qualified negative-hit signatures
############################################################

cmap_landscape_data <- integrated %>%
  
  dplyr::filter(
    has_cMAP,
    is.finite(median_norm_cs),
    is.finite(median_fdr_q_nlog10)
  ) %>%
  
  annotate_landscape_group() %>%
  
  dplyr::mutate(
    
    reversal_magnitude =
      -median_norm_cs,
    
    fdr_support =
      median_fdr_q_nlog10
  )


cmap_landscape_labels <-
  make_landscape_labels(
    cmap_landscape_data,
    n_per_group =
      CMAP_LANDSCAPE_LABEL_N_PER_GROUP
  )


pB_cmap_landscape <- ggplot(
  cmap_landscape_data,
  aes(
    x = reversal_magnitude,
    y = fdr_support
  )
) +
  
  ## All qualified cMAP compounds
  geom_point(
    data = cmap_landscape_data %>%
      dplyr::filter(
        display_group ==
          "Other qualified compounds"
      ),
    color = "grey76",
    size = 1.05,
    alpha = 0.48
  ) +
  
  ## Post-scoring highlighted biological categories
  geom_point(
    data = cmap_landscape_data %>%
      dplyr::filter(
        display_group !=
          "Other qualified compounds"
      ),
    aes(
      fill = display_group
    ),
    shape = 21,
    size = 2.55,
    color = "grey20",
    stroke = 0.38,
    alpha = 0.96,
    show.legend = TRUE
  ) +
  
  ## Experimental validation compounds:
  ## red ring only; selection is post-scoring
  geom_point(
    data = cmap_landscape_data %>%
      dplyr::filter(
        is_selected
      ),
    aes(
      color = "Selected for validation"
    ),
    shape = 21,
    size = 3.80,
    fill = NA,
    stroke = 0.95
  ) +
  
  ggrepel::geom_text_repel(
    data = cmap_landscape_labels %>%
      dplyr::filter(
        !drug_std %in% c(
          "albendazole",
          "nocodazole",
          "dolastatin-10"
        )
      ),
    aes(
      label = drug_name
    ),
    size = 2.20,
    color = "black",
    box.padding = 0.34,
    point.padding = 0.40,
    force = 1.15,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 20260710
  ) +
  
  ## Keep the albendazole label clear of its point and selection ring.
  ggrepel::geom_text_repel(
    data = cmap_landscape_labels %>%
      dplyr::filter(
        drug_std == "albendazole"
      ),
    aes(
      label = drug_name
    ),
    size = 2.20,
    color = "black",
    nudge_x = 0,
    nudge_y = -0.70,
    hjust = 1,
    box.padding = 0.20,
    point.padding = 0.45,
    force = 1.10,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.25,
    max.overlaps = Inf,
    seed = 20260710
  ) +
  ggrepel::geom_text_repel(
    data = cmap_landscape_labels %>%
      dplyr::filter(
        drug_std == "nocodazole"
      ),
    aes(
      label = drug_name
    ),
    size = 2.20,
    color = "black",
    nudge_x = 0,
    nudge_y = 0.50,
    hjust = 1,
    direction = "y",
    box.padding = 0.18,
    point.padding = 0.38,
    force = 1.05,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.25,
    max.overlaps = Inf,
    seed = 20260710
  ) +
  ## Dolastatin-10: move label to the right of its point
  ggrepel::geom_text_repel(
    data = cmap_landscape_labels %>%
      dplyr::filter(
        drug_std == "dolastatin-10"
      ),
    aes(
      label = drug_name
    ),
    size = 2.20,
    color = "black",
    nudge_x = 0,
    nudge_y = 0.2,
    hjust = 0,
    direction = "y",
    box.padding = 0.18,
    point.padding = 0.35,
    force = 1.05,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.25,
    max.overlaps = Inf,
    seed = 20260710
  ) +
  scale_fill_manual(
    values = LANDSCAPE_FILL_COLORS,
    labels = c(
      "Tubulin / microtubule" = "MT-related",
      "Benzimidazole class" = "Benzimidazole"
    ),
    drop = FALSE,
    name = NULL
  ) +
  
  scale_color_manual(
    values = c(
      "Selected for validation" =
        SELECTED_RING_COLOR
    ),
    labels = c(
      "Selected for validation" = "Selected"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.04,
        0.12
      )
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.10
      ),
      add = c(
        0.30,
        0
      )
    )
  ) +
  
  labs(
    title = "cMAP qualified negative-hit landscape",
    x = "Reversal magnitude (-median normCS)",
    y = "Median FDR evidence (-log10 q)"
  ) +
  
  guides(
    
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        size = 2.4,
        fill = unname(
          LANDSCAPE_FILL_COLORS
        ),
        alpha = 1
      )
    ),
    
    color = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        size = 3.2,
        fill = NA,
        stroke = 0.9
      )
    )
  ) +
  
  theme_sci(
    base_size = 8.5
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 8.7,
      hjust = 0
    ),
    
    axis.title = element_text(
      size = 7.5
    ),
    
    axis.text = element_text(
      size = 6.8
    ),
    
    legend.position = "right",
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.text = element_text(
      size = 5.9
    ),
    legend.key.width = grid::unit(
      0.12,
      "in"
    ),
    legend.key.height = grid::unit(
      0.12,
      "in"
    ),
    legend.spacing.y = grid::unit(
      0.02,
      "in"
    ),
    
    plot.margin = margin(
      5,
      6,
      5,
      5
    )
  )


############################################################
## L1000FWD reversal landscape
##
## x:
## negative median similarity -> positive reversal magnitude
##
## y:
## fraction of eligible signatures with similarity < 0
############################################################

l1000_landscape_data <- integrated %>%
  
  dplyr::filter(
    has_L1000FWD,
    is.finite(median_similarity),
    is.finite(l1000_negative_fraction)
  ) %>%
  
  annotate_landscape_group() %>%
  
  dplyr::mutate(
    
    reversal_magnitude =
      -median_similarity,
    
    negative_fraction =
      l1000_negative_fraction
  )


l1000_landscape_labels <-
  make_landscape_labels(
    l1000_landscape_data,
    n_per_group =
      L1000_LANDSCAPE_LABEL_N_PER_GROUP
  )


pB_l1000_landscape <- ggplot(
  l1000_landscape_data,
  aes(
    x = reversal_magnitude,
    y = negative_fraction
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.38,
    color = "grey55"
  ) +
  
  ## All eligible L1000FWD compounds
  geom_point(
    data = l1000_landscape_data %>%
      dplyr::filter(
        display_group ==
          "Other qualified compounds"
      ),
    color = "grey76",
    size = 1.05,
    alpha = 0.48
  ) +
  
  ## Post-scoring biological categories
  geom_point(
    data = l1000_landscape_data %>%
      dplyr::filter(
        display_group !=
          "Other qualified compounds"
      ),
    aes(
      fill = display_group
    ),
    shape = 21,
    size = 2.55,
    color = "grey20",
    stroke = 0.38,
    alpha = 0.96,
    show.legend = TRUE
  ) +
  
  ## Selected validation compounds
  geom_point(
    data = l1000_landscape_data %>%
      dplyr::filter(
        is_selected
      ),
    aes(
      color = "Selected for validation"
    ),
    shape = 21,
    size = 3.80,
    fill = NA,
    stroke = 0.95
  ) +
  
  ggrepel::geom_text_repel(
    data = l1000_landscape_labels %>%
      dplyr::filter(
        drug_std != "parbendazole"
      ),
    aes(
      label = drug_name
    ),
    size = 2.20,
    color = "black",
    box.padding = 0.34,
    point.padding = 0.40,
    force = 1.15,
    xlim = c(-0.03, 0.13),
    min.segment.length = 0,
    segment.color = "grey55",
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 20260710
  ) +
  
  ## Keep parbendazole close to its observed data point.
  ggrepel::geom_text_repel(
    data = l1000_landscape_labels %>%
      dplyr::filter(
        drug_std == "parbendazole"
      ),
    aes(
      label = drug_name
    ),
    size = 2.20,
    color = "black",
    nudge_x = -0.012,
    nudge_y = 0.015,
    hjust = 1,
    box.padding = 0.18,
    point.padding = 0.35,
    force = 1.15,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 20260710
  ) +
  
  scale_fill_manual(
    values = LANDSCAPE_FILL_COLORS,
    labels = c(
      "Tubulin / microtubule" = "MT-related",
      "Benzimidazole class" = "Benzimidazole"
    ),
    drop = FALSE,
    name = NULL
  ) +
  
  scale_color_manual(
    values = c(
      "Selected for validation" =
        SELECTED_RING_COLOR
    ),
    labels = c(
      "Selected for validation" = "Selected"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.05,
        0.13
      )
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      0.25
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.04
      )
    )
  ) +
  
  labs(
    title = "L1000FWD reversal landscape",
    x = "Reversal magnitude (-median similarity)",
    y = "Negative-signature fraction"
  ) +
  
  guides(
    
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        size = 2.4,
        fill = unname(
          LANDSCAPE_FILL_COLORS
        ),
        alpha = 1
      )
    ),
    
    color = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        size = 3.2,
        fill = NA,
        stroke = 0.9
      )
    )
  ) +
  
  theme_sci(
    base_size = 8.5
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 8.7,
      hjust = 0
    ),
    
    axis.title = element_text(
      size = 7.5
    ),
    
    axis.text = element_text(
      size = 6.8
    ),
    
    legend.position = "right",
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.text = element_text(
      size = 5.9
    ),
    legend.key.width = grid::unit(
      0.12,
      "in"
    ),
    legend.key.height = grid::unit(
      0.12,
      "in"
    ),
    legend.spacing.y = grid::unit(
      0.02,
      "in"
    ),
    
    plot.margin = margin(
      5,
      6,
      5,
      5
    )
  )


############################################################
## Combined Figure 3B
############################################################

pB_main <- patchwork::wrap_plots(
  
  pB_cmap_landscape +
    theme(
      legend.position = "none"
    ),
  
  pB_l1000_landscape,
  
  ncol = 2,
  
  widths = c(
    1,
    1.12
  ),
  
  guides = "keep"
  
) +
  
  patchwork::plot_annotation(
    
    title =
      "Perturbational reversal evidence landscapes",
    
    theme = theme(
      
      plot.title = element_text(
        face = "bold",
        size = 9.2,
        hjust = 0
      )
    )
  )


############################################################
## Export Figure 3B
############################################################

save_plot_both(
  pB_cmap_landscape,
  "Fig3B_left_cMAP_negative_hit_landscape",
  width = 4.05,
  height = 3.45
)

save_plot_both(
  make_clean_plot(
    pB_cmap_landscape
  ),
  "Fig3B_left_cMAP_negative_hit_landscape_clean",
  width = 4.05,
  height = 3.45
)


save_plot_both(
  pB_l1000_landscape,
  "Fig3B_right_L1000FWD_reversal_landscape",
  width = 4.05,
  height = 3.45
)

save_plot_both(
  make_clean_plot(
    pB_l1000_landscape
  ),
  "Fig3B_right_L1000FWD_reversal_landscape_clean",
  width = 4.05,
  height = 3.45
)


save_plot_both(
  pB_main,
  "Fig3B_perturbational_reversal_landscapes",
  width = 7.1,
  height = 3.55
)

save_plot_both(
  make_clean_plot(
    pB_main
  ),
  "Fig3B_perturbational_reversal_landscapes_clean",
  width = 7.1,
  height = 3.55
)
############################################################
## 20. Main Figure C:
## Mechanism-axis lollipop plot
############################################################

axis_plot <- axis_summary %>%
  dplyr::filter(
    n_cmap_drugs >= MIN_AXIS_DRUGS_FOR_PLOT,
    n_same_drug_L1000 >= 3,
    !cmap_mechanism_axis %in% c(
      "Unannotated",
      "Other annotated mechanism"
    ),
    is.finite(median_cmap_hit_score),
    is.finite(median_L1000_score_same_drug)
  ) %>%
  dplyr::arrange(
    dplyr::desc(median_cmap_hit_score)
  ) %>%
  dplyr::slice_head(
    n = TOP_N_AXES_FOR_PLOT
  ) %>%
  dplyr::mutate(
    is_microtubule =
      cmap_mechanism_axis == "Tubulin / microtubule axis"
  )

if (nrow(axis_plot) == 0) {
  stop("No mechanism axes satisfy Figure 3C display criteria.")
}

axis_plot2 <- axis_plot %>%
  dplyr::mutate(
    axis_label = dplyr::recode(
      cmap_mechanism_axis,
      "Tubulin / microtubule axis" = "Tubulin / microtubule",
      "Glucocorticoid receptor axis" = "Glucocorticoid receptor",
      "Calcium channel blocker" = "Calcium channel",
      "Na/K-ATPase axis" = "Na/K-ATPase",
      "Retinoid receptor axis" = "Retinoid receptor",
      "Topoisomerase inhibitor" = "Topoisomerase inhibitor",
      "HDAC inhibitor" = "HDAC inhibitor",
      "HMGCR inhibitor" = "HMGCR inhibitor",
      "EGFR inhibitor" = "EGFR inhibitor",
      "CDK inhibitor" = "CDK inhibitor",
      "JAK inhibitor" = "JAK inhibitor",
      "RAF inhibitor" = "RAF inhibitor",
      "MEK inhibitor" = "MEK inhibitor",
      "mTOR inhibitor" = "mTOR inhibitor",
      "NF-kB axis" = "NF-kB axis"
    ),
    axis_label = ifelse(
      is.na(axis_label),
      cmap_mechanism_axis,
      axis_label
    ),
    axis_label = forcats::fct_reorder(
      axis_label,
      median_cmap_hit_score
    )
  )

pC <- ggplot(
  axis_plot2,
  aes(
    x = median_cmap_hit_score,
    y = axis_label
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = median_cmap_hit_score,
      y = axis_label,
      yend = axis_label
    ),
    linewidth = 0.42,
    color = "grey75"
  ) +
  geom_point(
    aes(
      size = n_same_drug_L1000,
      fill = median_L1000_score_same_drug
    ),
    shape = 21,
    color = "grey25",
    stroke = 0.30,
    alpha = 0.92
  ) +
  geom_point(
    data = axis_plot2 %>%
      dplyr::filter(
        is_microtubule
      ),
    aes(
      x = median_cmap_hit_score,
      y = axis_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = NA,
    color = "#B2182B",
    stroke = 1.05,
    size = 5.4
  ) +
  scale_fill_viridis_c(
    option = "C",
    direction = -1,
    limits = c(0, 1),
    name = "L1000FWD\nscore"
  ) +
  scale_size_continuous(
    range = c(2.4, 5.3),
    breaks = c(5, 10, 15),
    name = "Matched\ndrugs"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(
      mult = c(0.01, 0.06)
    )
  ) +
  labs(
    title = "Mechanism-axis evidence from\ncMAP and L1000FWD",
    x = "Median cMAP evidence score",
    y = NULL
  ) +
  guides(
    fill = guide_colorbar(
      order = 1,
      barheight = grid::unit(0.55, "in"),
      barwidth = grid::unit(0.08, "in")
    ),
    size = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        fill = "grey75",
        color = "grey25"
      )
    )
  ) +
  theme_sci(
    base_size = 8.5
  ) +
  theme(
    axis.text.y = element_text(
      size = 7.0,
      color = "black"
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 6.5,
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 8.7,
      hjust = 0,
      lineheight = 0.95,
      margin = margin(b = 4)
    ),
    legend.text = element_text(
      size = 6.0
    )
  )

save_plot_both(
  pC,
  "Fig3C_mechanism_axis_cross_resource",
  width = 3.65,
  height = 3.55
)

save_plot_both(
  make_clean_plot(pC),
  "Fig3C_mechanism_axis_cross_resource_clean",
  width = 3.65,
  height = 3.55
)

############################################################
## 21. Main Figure D:
## Horizontal L1000FWD boxplots for benzimidazole compounds
############################################################

if (nrow(benzo_l1000_summary) == 0) {
  warning("No eligible curated benzimidazole-class drugs with n>=3 signatures; Fig3D not generated.")
} else {
  
  benzo_other_order <- benzo_l1000_summary %>%
    dplyr::filter(
      !drug_std %in% SELECTED_DRUGS
    ) %>%
    dplyr::arrange(
      median_similarity
    ) %>%
    dplyr::pull(
      drug_std
    )
  
  benzo_display_order_top <- unique(
    c(
      "albendazole",
      "mebendazole",
      benzo_other_order
    )
  )
  
  benzo_display_order_top <- benzo_display_order_top[
    benzo_display_order_top %in% benzo_l1000_summary$drug_std
  ]
  
  benzo_label_map <- benzo_l1000_summary %>%
    dplyr::mutate(
      drug_label = paste0(
        drug_l1000_name,
        "\n(n=",
        n_l1000_signatures,
        "; neg=",
        scales::percent(
          l1000_negative_fraction,
          accuracy = 1
        ),
        ")"
      )
    ) %>%
    dplyr::select(
      drug_std,
      drug_label
    )
  
  benzo_l1000_raw_plot <- benzo_l1000_raw %>%
    dplyr::left_join(
      benzo_label_map,
      by = "drug_std"
    ) %>%
    dplyr::mutate(
      drug_std = factor(
        drug_std,
        levels = rev(benzo_display_order_top)
      ),
      plot_group = ifelse(
        selected_for_validation,
        "Albendazole / Mebendazole",
        "Other eligible benzimidazole compounds"
      )
    )
  
  pD <- ggplot(
    benzo_l1000_raw_plot,
    aes(
      y = drug_std,
      x = similarity_scores,
      color = plot_group
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.38,
      color = "grey40"
    ) +
    geom_boxplot(
      aes(
        fill = plot_group
      ),
      width = 0.52,
      outlier.shape = NA,
      alpha = 0.18,
      linewidth = 0.50,
      show.legend = FALSE
    ) +
    geom_jitter(
      height = 0.12,
      width = 0,
      size = 1.05,
      alpha = 0.62,
      show.legend = TRUE
    ) +
    scale_y_discrete(
      labels = setNames(
        benzo_label_map$drug_label,
        benzo_label_map$drug_std
      )
    ) +
    scale_color_manual(
      values = c(
        "Albendazole / Mebendazole" = "#B2182B",
        "Other eligible benzimidazole compounds" = "#EF8A62"
      ),
      labels = c(
        "Albendazole / Mebendazole" = "ABZ/MBZ",
        "Other eligible benzimidazole compounds" = "Other compounds"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Albendazole / Mebendazole" = "#B2182B",
        "Other eligible benzimidazole compounds" = "#EF8A62"
      ),
      guide = "none"
    ) +
    coord_cartesian(
      xlim = c(-0.22, 0.12),
      clip = "off"
    ) +
    labs(
      title = "L1000FWD profiles of curated benzimidazole compounds",
      x = "L1000FWD similarity score",
      y = NULL,
      color = NULL
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(
          shape = 16,
          size = 1.8,
          alpha = 1
        )
      )
    ) +
    theme_sci(
      base_size = 8.5
    ) +
    theme(
      legend.position = "right",
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      axis.text.y = element_text(
        size = 6.7,
        lineheight = 0.95
      ),
      legend.text = element_text(
        size = 5.9
      ),
      legend.key.width = grid::unit(
        0.12,
        "in"
      ),
      legend.key.height = grid::unit(
        0.12,
        "in"
      ),
      legend.spacing.y = grid::unit(
        0.02,
        "in"
      ),
      plot.margin = margin(
        5,
        7,
        5,
        5
      )
    )
  
  save_plot_both(
    pD,
    "Fig3D_benzimidazole_raw_points_median",
    width = 5.65,
    height = 3.7
  )
  
  save_plot_both(
    make_clean_plot(pD),
    "Fig3D_benzimidazole_raw_points_median_clean",
    width = 5.65,
    height = 3.7
  )
}

############################################################
## 22. Supplementary Figure:
## Score sensitivity plots
############################################################

pS_sens_cmap <- ggplot(
  cmap_score_sensitivity,
  aes(
    x = cMAP_hit_score_primary,
    y = cMAP_hit_score_weighted_sensitivity
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.42,
    color = "grey55"
  ) +
  geom_point(
    color = "grey55",
    alpha = 0.50,
    size = 1.15
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "cMAP score sensitivity",
    subtitle = paste0(
      "Spearman rho = ",
      sprintf("%.2f", rho_cmap_primary_vs_sensitivity)
    ),
    x = "Primary robust cMAP score",
    y = "Legacy weighted cMAP score"
  ) +
  theme_sci(
    base_size = 8.5
  )

pS_sens_shared <- ggplot(
  shared_sensitivity,
  aes(
    x = CrossResourceConcordanceScore,
    y = CrossResourceConcordanceScore_weighted_cMAP_sensitivity
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.42,
    color = "grey55"
  ) +
  geom_point(
    color = "grey55",
    alpha = 0.50,
    size = 1.15
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Cross-resource score sensitivity",
    subtitle = paste0(
      "Shared drugs only; Spearman rho = ",
      sprintf("%.2f", rho_shared_primary_vs_sensitivity)
    ),
    x = "Primary cross-resource concordance score",
    y = "Sensitivity concordance score"
  ) +
  theme_sci(
    base_size = 8.5
  )

FigS_score_sensitivity <- patchwork::wrap_plots(
  pS_sens_cmap,
  pS_sens_shared,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Sensitivity of perturbational evidence scores",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 11,
        hjust = 0
      )
    )
  )

save_plot_both(
  FigS_score_sensitivity,
  "FigS_score_sensitivity",
  width = 7.1,
  height = 3.3
)

save_plot_both(
  make_clean_plot(FigS_score_sensitivity),
  "FigS_score_sensitivity_clean",
  width = 7.1,
  height = 3.3
)
############################################################
## OPTIONAL APPEND BLOCK:
## Export the two score-sensitivity panels separately
##
## Paste this block AFTER Section 22 of:
##
## It does NOT recalculate or alter any score.
## It only exports the already-created objects:
##   pS_sens_cmap
##   pS_sens_shared
############################################################

required_objects <- c(
  "pS_sens_cmap",
  "pS_sens_shared",
  "save_plot_both",
  "make_clean_plot"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Run the drug-screening script through Section 22 first. ",
    "Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

save_plot_both(
  pS_sens_cmap,
  "FigS_sensitivity_primary_vs_legacy_cMAP",
  width = 3.55,
  height = 3.30
)

save_plot_both(
  make_clean_plot(pS_sens_cmap),
  "FigS_sensitivity_primary_vs_legacy_cMAP_clean",
  width = 3.55,
  height = 3.30
)

save_plot_both(
  pS_sens_shared,
  "FigS_sensitivity_cross_resource_concordance",
  width = 3.55,
  height = 3.30
)

save_plot_both(
  make_clean_plot(pS_sens_shared),
  "FigS_sensitivity_cross_resource_concordance_clean",
  width = 3.55,
  height = 3.30
)

message("Separate score-sensitivity panels exported to: ", FIG_DIR)
############################################################
## 23. Combined SCI Figure 3 at double-column width
############################################################

Fig3_bottom <- patchwork::wrap_plots(
  pC,
  pD,
  ncol = 2,
  widths = c(0.95, 1.25)
)

Fig3 <- patchwork::wrap_plots(
  pA,
  pB_main,
  Fig3_bottom,
  ncol = 1,
  heights = c(0.60, 1.55, 1.70)
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = paste0(
      "Perturbational reversal profiling highlights ",
      "microtubule-related signals and supports ",
      "post-scoring benzimidazole candidate selection"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 11.5,
        hjust = 0
      ),
      plot.tag = element_text(
        face = "bold",
        size = 11
      )
    )
  )

save_plot_both(
  Fig3,
  "Figure3_drug_screening_final",
  width = 7.1,
  height = 9.3,
  dpi = 600
)

Fig3_clean <- patchwork::wrap_plots(
  make_clean_plot(pA),
  make_clean_plot(pB_main),
  patchwork::wrap_plots(
    make_clean_plot(pC),
    make_clean_plot(pD),
    ncol = 2,
    widths = c(0.95, 1.25)
  ),
  ncol = 1,
  heights = c(0.60, 1.55, 1.70)
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 11
      )
    )
  )

save_plot_both(
  Fig3_clean,
  "Figure3_drug_screening_final_clean",
  width = 7.1,
  height = 9.3,
  dpi = 600
)
############################################################
## 24. Export figure data and session info
############################################################

readr::write_csv(plot_shared, file.path(TAB_DIR, "09_FigS_cross_resource_plot_data.csv"))
readr::write_csv(axis_plot, file.path(TAB_DIR, "10_Fig3C_mechanism_axis_plot_data.csv"))
readr::write_csv(benzo_l1000_summary, file.path(TAB_DIR, "11_Fig3D_benzimidazole_plot_summary.csv"))
readr::write_csv(candidate_score_long, file.path(TAB_DIR, "14_candidate_matrix_scores.csv"))
readr::write_csv(candidate_annotation_long, file.path(TAB_DIR, "15_candidate_matrix_annotations.csv"))
readr::write_csv(
  cmap_landscape_data,
  file.path(
    TAB_DIR,
    "17_Fig3B_cMAP_negative_hit_landscape.csv"
  )
)

readr::write_csv(
  l1000_landscape_data,
  file.path(
    TAB_DIR,
    "18_Fig3B_L1000FWD_reversal_landscape.csv"
  )
)

readr::write_csv(
  cmap_landscape_labels,
  file.path(
    TAB_DIR,
    "19_Fig3B_cMAP_landscape_labels.csv"
  )
)

readr::write_csv(
  l1000_landscape_labels,
  file.path(
    TAB_DIR,
    "20_Fig3B_L1000FWD_landscape_labels.csv"
  )
)


figure_export_qc <- .figure_export_log %>%
  dplyr::mutate(
    file_exists = file.exists(path),
    file_size_kb = ifelse(file_exists, round(file.info(path)$size / 1024, 1), NA_real_)
  )
readr::write_csv(figure_export_qc, file.path(TAB_DIR, "16_figure_export_QC.csv"))

sink(file.path(OUT_DIR, "sessionInfo.txt"))
sessionInfo()
sink()

############################################################
## 25. Console summary
############################################################

cat("\n============================================================\n")
cat("AML perturbational screening completed\n")
cat("============================================================\n")
cat("cMAP qualified unique signatures:", nrow(cmap_hit), "\n")
cat("cMAP unique drugs:", nrow(cmap_drug), "\n")
cat("L1000FWD eligible drugs (n>=3 unique signatures):", nrow(l1000_drug), "\n")
cat("Same-drug cross-resource compounds:", nrow(shared_drugs), "\n")
cat("Eligible curated benzimidazole-class drugs:", nrow(benzo_l1000_summary), "\n")
cat("\nSelected validation compounds (post-scoring):\n")
print(
  selected_validation_drugs %>%
    dplyr::select(
      drug_name, drug_std, evidence_scope,
      cMAP_hit_score, L1000FWD_reversal_score,
      CrossResourceConcordanceScore,
      structural_class, validation_status
    )
)
cat("\nOutput directory:\n", normalizePath(OUT_DIR), "\n")
cat("============================================================\n")
