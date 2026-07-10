#!/usr/bin/env Rscript
## Three NR/NV matrix schemes side by side: same VAF sparsity hexbins as analyze_nr_nv_pileup_depth.R
## (mean vs prevalence, max vs prevalence, mean vs max). No depth funnel, ECDFs, or comparison tables.
##
## Usage:
##   Rscript plot_nr_nv_sparsity_hexbin_three_schemes.R <group_id> <matrix_dir>
##
## Required in matrix_dir:
##   NR/NV_annotated_vcf_<group>_{unfiltered,pileup,pileup_hq_depth}.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_nr_nv_sparsity_hexbin_three_schemes.R <group_id> <matrix_dir>")
}

group_id   <- args[1]
matrix_dir <- args[2]

argv0 <- commandArgs(trailingOnly = FALSE)
file_arg <- argv0[grepl("^--file=", argv0)]
script_dir <- if (length(file_arg)) {
  dirname(sub("^--file=", "", file_arg[1]))
} else {
  "."
}
theme_candidates <- c(
  "/usr/local/bin/theme_ohchibi_pubr.R",
  file.path(script_dir, "theme_ohchibi_pubr.R"),
  "/home/ubuntu/projects/cursor_rules/theme_ohchibi_pubr.R"
)
theme_path <- theme_candidates[file.exists(theme_candidates)][1]
if (!is.na(theme_path) && nzchar(theme_path)) {
  source(theme_path)
} else {
  theme_ohchibi_pubr <- ggplot2::theme_bw
}

read_count_matrix <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, na.strings = c("", "NA", "."), data.table = TRUE)
  if (ncol(dt) < 2) stop("Matrix TSV must have variant column + >=1 sample column: ", path)
  vid <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  list(vid = vid, mat = mat, samples = colnames(dt)[-1])
}

sp_dt_from_paths <- function(nr_path, nv_path, scheme_label) {
  nr <- read_count_matrix(nr_path)
  nv <- read_count_matrix(nv_path)
  if (!identical(nr$vid, nv$vid)) stop("NR/NV variant IDs differ: ", nr_path)
  if (!identical(nr$samples, nv$samples)) stop("NR/NV sample columns differ.")
  nr_mat <- nr$mat
  nv_mat <- nv$mat
  n_samp <- ncol(nr_mat)
  vaf <- nv_mat / pmax(nr_mat, 1)
  vaf[nv_mat <= 0] <- NA
  mean_vaf <- rowMeans(vaf, na.rm = TRUE)
  max_vaf <- apply(vaf, 1L, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  prop_nv <- if (n_samp > 0) rowSums(nv_mat > 0) / n_samp else numeric(0)
  data.table(
    scheme = scheme_label,
    mean_vaf = mean_vaf,
    max_vaf = max_vaf,
    prop_cells_nv_gt0 = prop_nv
  )[is.finite(mean_vaf) & is.finite(max_vaf) & is.finite(prop_cells_nv_gt0)]
}

scheme_levels <- c("unfiltered", "pileup (Pass)", "pileup_hq_depth")
nr_unf <- file.path(matrix_dir, paste0("NR_annotated_vcf_", group_id, "_unfiltered.tsv"))
nv_unf <- file.path(matrix_dir, paste0("NV_annotated_vcf_", group_id, "_unfiltered.tsv"))
nr_pl <- file.path(matrix_dir, paste0("NR_annotated_vcf_", group_id, "_pileup.tsv"))
nv_pl <- file.path(matrix_dir, paste0("NV_annotated_vcf_", group_id, "_pileup.tsv"))
nr_hq <- file.path(matrix_dir, paste0("NR_annotated_vcf_", group_id, "_pileup_hq_depth.tsv"))
nv_hq <- file.path(matrix_dir, paste0("NV_annotated_vcf_", group_id, "_pileup_hq_depth.tsv"))

req <- c(nr_unf, nv_unf, nr_pl, nv_pl, nr_hq, nv_hq)
miss <- req[!file.exists(req)]
if (length(miss)) stop("Missing matrix files:\n", paste(miss, collapse = "\n"))

sp_all <- rbindlist(list(
  sp_dt_from_paths(nr_unf, nv_unf, scheme_levels[1]),
  sp_dt_from_paths(nr_pl, nv_pl, scheme_levels[2]),
  sp_dt_from_paths(nr_hq, nv_hq, scheme_levels[3])
))
sp_all[, scheme := factor(scheme, levels = scheme_levels)]

out_pdf <- paste0("pileup_sparsity_hexbin_three_schemes_", group_id, ".pdf")
if (nrow(sp_all) < 2) {
  pdf(out_pdf, width = 12, height = 4)
  plot.new()
  text(0.5, 0.5, "Insufficient variants for hexbin", cex = 1.2)
  dev.off()
  message("[plot_nr_nv_sparsity_hexbin_three_schemes] Wrote (empty): ", out_pdf)
  quit(save = "no")
}

use_hex <- requireNamespace("hexbin", quietly = TRUE)
kovesi_cols <- tryCatch(pals::kovesi.rainbow(256), error = function(e) NULL)
hex_fill_log10 <- aes(fill = after_stat(log10(pmax(count, 1))))
hex_legend_counts <- function(x) label_comma(accuracy = 1)(pmax(1, round(10^as.numeric(x))))

hex_fill_scale <- function() {
  if (!is.null(kovesi_cols) && length(kovesi_cols) > 1) {
    scale_fill_gradientn(
      colours = kovesi_cols,
      name = "Variants\nper hex",
      labels = hex_legend_counts,
      na.value = "grey92"
    )
  } else {
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#2166AC", "#08306B"),
      name = "Variants\nper hex",
      labels = hex_legend_counts,
      na.value = "grey92"
    )
  }
}

br01 <- seq(0, 1, 0.1)
scale_hex_xy <- function() {
  list(
    scale_x_continuous(breaks = br01, limits = c(0, 1), expand = c(0, 0)),
    scale_y_continuous(breaks = br01, limits = c(0, 1), expand = c(0, 0))
  )
}

facet_three <- facet_wrap(~scheme, nrow = 1L, scales = "fixed")

pdf(out_pdf, width = 14, height = 5)

p1 <- ggplot(sp_all, aes(mean_vaf, prop_cells_nv_gt0)) +
  { if (use_hex) geom_hex(hex_fill_log10, bins = 45) else geom_bin2d(hex_fill_log10, bins = 45) } +
  hex_fill_scale() +
  scale_hex_xy() +
  facet_three +
  labs(
    title = paste0("VAF sparsity by matrix scheme — ", group_id),
    subtitle = "Mean VAF (NV>0 cells only) vs proportion of cells with NV > 0",
    x = "Mean VAF (NV>0 cells only)",
    y = "Proportion of cells with NV > 0"
  ) +
  theme_ohchibi_pubr()
print(p1)

p2 <- ggplot(sp_all, aes(max_vaf, prop_cells_nv_gt0)) +
  { if (use_hex) geom_hex(hex_fill_log10, bins = 45) else geom_bin2d(hex_fill_log10, bins = 45) } +
  hex_fill_scale() +
  scale_hex_xy() +
  facet_three +
  labs(
    title = paste0("VAF sparsity (max VAF) — ", group_id),
    subtitle = "Max VAF over NV>0 cells vs proportion of cells with NV > 0",
    x = "Max VAF (NV>0 cells only)",
    y = "Proportion of cells with NV > 0"
  ) +
  theme_ohchibi_pubr()
print(p2)

p3 <- ggplot(sp_all, aes(mean_vaf, max_vaf)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
  { if (use_hex) geom_hex(hex_fill_log10, bins = 45) else geom_bin2d(hex_fill_log10, bins = 45) } +
  hex_fill_scale() +
  scale_hex_xy() +
  facet_three +
  labs(
    title = paste0("Mean vs max VAF (NV>0 cells only) — ", group_id),
    subtitle = "Dashed line y = x",
    x = "Mean VAF (NV>0 cells only)",
    y = "Max VAF (NV>0 cells only)"
  ) +
  theme_ohchibi_pubr()
print(p3)

dev.off()

message("[plot_nr_nv_sparsity_hexbin_three_schemes] Wrote: ", out_pdf)
