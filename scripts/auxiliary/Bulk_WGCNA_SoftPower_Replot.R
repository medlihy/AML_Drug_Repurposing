############################################################
## GSE114868 WGCNA soft-power replot only
## Purpose:
##   Re-draw soft-threshold diagnostics without overlapping labels
##   and avoid corrupted/incomplete PDF devices.
##
## Inputs expected from previous WGCNA output:
##   02_soft_power/01_soft_threshold_fit_indices.tsv
##   02_soft_power/02_selected_soft_power.txt
##
## Outputs:
##   02_soft_power/03_soft_power_diagnostics_clean_no_numbers.pdf
##   02_soft_power/03b_soft_power_diagnostics_selected_label.pdf
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
POWER_DIR <- file.path(OUTDIR, "02_soft_power")

SFT_FILE <- file.path(POWER_DIR, "01_soft_threshold_fit_indices.tsv")
POWER_FILE <- file.path(POWER_DIR, "02_selected_soft_power.txt")

## scale-free topology fit threshold shown in the plot
R2_CUTOFF <- 0.85

## =========================
## 1. Helper functions
## =========================

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not installed: ", pkg)
  }
}

require_pkg("ggplot2")
suppressPackageStartupMessages(library(ggplot2))

safe_ggsave_pdf <- function(filename, plot, width, height) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
  message("Saved: ", filename)
}

find_col <- function(x, candidates, required = TRUE) {
  cn <- colnames(x)
  hit <- cn[tolower(cn) %in% tolower(candidates)]
  if (length(hit) > 0) return(hit[1])

  ## second pass: partial matching
  pattern <- paste(candidates, collapse = "|")
  hit <- grep(pattern, cn, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) return(hit[1])

  if (required) {
    stop("Cannot find required column. Candidates: ", paste(candidates, collapse = ", "),
         "\nAvailable columns: ", paste(cn, collapse = ", "))
  }
  NA_character_
}

## =========================
## 2. Read soft-power result table
## =========================

if (!file.exists(SFT_FILE)) {
  stop("Soft-threshold table not found: ", SFT_FILE)
}

sft_tab <- read.table(
  SFT_FILE,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

power_col <- find_col(sft_tab, c("Power", "power"))
r2_col    <- find_col(sft_tab, c("signed_R2", "signed.R2", "signed R2"), required = FALSE)
sft_r2_col <- find_col(sft_tab, c("SFT.R.sq", "SFT R sq", "SFT.R.sq."), required = FALSE)
slope_col <- find_col(sft_tab, c("slope"), required = FALSE)
mean_k_col <- find_col(sft_tab, c("mean.k.", "mean.k", "mean connectivity", "mean.k.."), required = FALSE)

sft_tab$Power <- as.numeric(sft_tab[[power_col]])

if (!is.na(r2_col)) {
  sft_tab$signed_R2 <- as.numeric(sft_tab[[r2_col]])
} else if (!is.na(sft_r2_col) && !is.na(slope_col)) {
  ## WGCNA classic signed scale-free fit:
  ## sign(slope) * R^2 is usually negative, so use -sign(slope) * R^2.
  sft_tab$signed_R2 <- -sign(as.numeric(sft_tab[[slope_col]])) * as.numeric(sft_tab[[sft_r2_col]])
} else if (!is.na(sft_r2_col)) {
  warning("slope column not found; using SFT.R.sq as signed_R2 approximation.")
  sft_tab$signed_R2 <- as.numeric(sft_tab[[sft_r2_col]])
} else {
  stop("Cannot construct signed_R2 from soft-threshold table.")
}

if (is.na(mean_k_col)) {
  stop("Cannot find mean connectivity column. Available columns: ",
       paste(colnames(sft_tab), collapse = ", "))
}
sft_tab$mean_connectivity <- as.numeric(sft_tab[[mean_k_col]])

sft_tab <- sft_tab[order(sft_tab$Power), , drop = FALSE]

## Selected power
if (file.exists(POWER_FILE)) {
  softPower <- suppressWarnings(as.numeric(readLines(POWER_FILE, warn = FALSE)[1]))
} else {
  warning("Selected soft-power file not found: ", POWER_FILE,
          "\nWill choose the first power with signed_R2 >= ", R2_CUTOFF)
  candidate <- sft_tab$Power[sft_tab$signed_R2 >= R2_CUTOFF]
  softPower <- if (length(candidate) > 0) candidate[1] else sft_tab$Power[which.max(sft_tab$signed_R2)]
}

if (is.na(softPower)) {
  softPower <- sft_tab$Power[which.max(sft_tab$signed_R2)]
  warning("Selected soft power is NA. Using power with max signed_R2: ", softPower)
}

selected_df <- sft_tab[sft_tab$Power == softPower, , drop = FALSE]
if (nrow(selected_df) == 0) {
  warning("Selected power ", softPower, " is not found in sft_tab. The vertical line will still be drawn.")
}

## =========================
## 3. Clean plot: no numeric labels
## =========================

p1 <- ggplot(sft_tab, aes(x = Power, y = signed_R2)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = R2_CUTOFF, linetype = "dashed", linewidth = 0.35) +
  geom_vline(xintercept = softPower, linetype = "dashed", linewidth = 0.35) +
  annotate(
    "text",
    x = softPower,
    y = min(1, max(sft_tab$signed_R2, na.rm = TRUE) + 0.04),
    label = paste0("power = ", softPower),
    hjust = -0.05,
    vjust = 0,
    size = 3.2
  ) +
  scale_x_continuous(
    breaks = sort(unique(c(seq(0, max(sft_tab$Power, na.rm = TRUE), by = 5), softPower))),
    minor_breaks = NULL
  ) +
  coord_cartesian(ylim = c(0, min(1, max(0.9, max(sft_tab$signed_R2, na.rm = TRUE) + 0.08)))) +
  labs(
    title = "Scale independence",
    x = "Soft threshold power",
    y = expression(paste("Signed scale-free topology fit ", R^2))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

p2 <- ggplot(sft_tab, aes(x = Power, y = mean_connectivity)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = softPower, linetype = "dashed", linewidth = 0.35) +
  annotate(
    "text",
    x = softPower,
    y = max(sft_tab$mean_connectivity, na.rm = TRUE),
    label = paste0("power = ", softPower),
    hjust = -0.05,
    vjust = 1.1,
    size = 3.2
  ) +
  scale_x_continuous(
    breaks = sort(unique(c(seq(0, max(sft_tab$Power, na.rm = TRUE), by = 5), softPower))),
    minor_breaks = NULL
  ) +
  labs(
    title = "Mean connectivity",
    x = "Soft threshold power",
    y = "Mean connectivity"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

## Use grid package to arrange two ggplots without extra dependencies
combined_file <- file.path(POWER_DIR, "03_soft_power_diagnostics_clean_no_numbers.pdf")
grDevices::cairo_pdf(combined_file, width = 11, height = 5.5)
tryCatch({
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2)))
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}, finally = {
  grDevices::dev.off()
})
message("Saved: ", combined_file)

## =========================
## 4. Optional selected-label version
##    Only selected power and threshold crossing are labelled.
## =========================

first_r2_power <- sft_tab$Power[which(sft_tab$signed_R2 >= R2_CUTOFF)[1]]
if (length(first_r2_power) == 0 || is.na(first_r2_power)) {
  first_r2_power <- softPower
}

label_df <- unique(data.frame(
  Power = c(first_r2_power, softPower),
  label = c(paste0("first R2≥", R2_CUTOFF, ": ", first_r2_power),
            paste0("selected: ", softPower)),
  stringsAsFactors = FALSE
))
label_df <- merge(label_df, sft_tab[, c("Power", "signed_R2", "mean_connectivity")],
                  by = "Power", all.x = TRUE)

p1_label <- p1 +
  geom_point(data = label_df, aes(x = Power, y = signed_R2), size = 2.5) +
  geom_text(
    data = label_df,
    aes(x = Power, y = signed_R2, label = label),
    nudge_y = 0.06,
    hjust = 0,
    size = 3.0,
    inherit.aes = FALSE
  )

p2_label <- p2 +
  geom_point(data = label_df, aes(x = Power, y = mean_connectivity), size = 2.5) +
  geom_text(
    data = label_df,
    aes(x = Power, y = mean_connectivity, label = label),
    nudge_y = max(sft_tab$mean_connectivity, na.rm = TRUE) * 0.04,
    hjust = 0,
    size = 3.0,
    inherit.aes = FALSE
  )

label_file <- file.path(POWER_DIR, "03b_soft_power_diagnostics_selected_label.pdf")
grDevices::cairo_pdf(label_file, width = 11, height = 5.5)
tryCatch({
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2)))
  print(p1_label, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2_label, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}, finally = {
  grDevices::dev.off()
})
message("Saved: ", label_file)

## =========================
## 5. Save cleaned soft-power table for checking
## =========================

write.table(
  sft_tab,
  file = file.path(POWER_DIR, "01b_soft_threshold_fit_indices_cleaned_for_plot.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Done. Soft-power figures saved in: ", normalizePath(POWER_DIR, mustWork = FALSE))
