#!/usr/bin/env Rscript
## NR/NV matrix depth–coverage + VAF sparsity for one or more matrix schemes (e.g. unfiltered, pileup, pileup_hq_depth).
## Q1/Q2: cohort + per-SC counts (outputs include column `scheme`).
## Q3: hexbins — one PDF, four sheets: (1) 2×3 mean vs prop + max vs prop by scheme; (2) 1×3 mean vs max;
## (3) 1×3 max VAF vs max NR depth at max-VAF cell (NR maxed across schemes at tied cells);
## (4) 1×3 max VAF vs count of cells with NV ≥ 1.
## Shared legend per sheet where applicable; sheets 1–2 use x/y (0–1) for VAF axes; sheet 3 x = depth; sheet 4 y = cell count.
## Default 3 schemes.
##
## Usage:
##   Rscript analyze_nr_nv_pileup_depth.R <group_id> <min_sample_pct> <depth_thresholds_csv> <nr1> <nv1> <label1> [<nr2> <nv2> <label2> ...]
##   Each scheme is a triple: NR path, NV path, short label (e.g. unfiltered, pileup, pileup_hq_depth).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop(
    "Usage: Rscript analyze_nr_nv_pileup_depth.R <group_id> <min_sample_pct> ",
    "<depth_thresholds_csv> <nr1> <nv1> <label1> [<nr2> <nv2> <label2> ...]"
  )
}

group_id       <- args[1]
min_sample_pct <- as.numeric(args[2])
if (is.na(min_sample_pct) || min_sample_pct <= 0 || min_sample_pct > 100) {
  stop("min_sample_pct must be in (0, 100]")
}

parse_depth_thresholds <- function(s) {
  if (is.null(s) || !nzchar(trimws(s))) {
    return(seq_len(10L))
  }
  parts <- strsplit(trimws(s), "\\s*,\\s*")[[1]]
  parts <- parts[nzchar(parts)]
  d <- suppressWarnings(as.integer(parts))
  d <- d[!is.na(d) & d > 0L]
  if (length(d) == 0) {
    stop("pileup depth thresholds: no positive integers after parsing: ", s)
  }
  d[!duplicated(d)]
}

depth_thresholds_csv <- args[3]
depths <- parse_depth_thresholds(depth_thresholds_csv)

rest <- args[-(1:3)]
if (length(rest) %% 3L != 0L) {
  stop("After depth_thresholds_csv, arguments must be triples: nr_path nv_path label (got ", length(rest), " extra args)")
}

scheme_defs <- list()
for (k in seq_len(length(rest) / 3L)) {
  i <- (k - 1L) * 3L
  scheme_defs[[k]] <- list(nr = rest[i + 1L], nv = rest[i + 2L], label = rest[i + 3L])
}

argv0 <- commandArgs(trailingOnly = FALSE)
file_arg <- argv0[grepl("^--file=", argv0)]
script_dir <- if (length(file_arg)) {
  dirname(sub("^--file=", "", file_arg[1]))
} else {
  "."
}
theme_candidates <- c(
  "/usr/local/bin/theme_ohchibi_pubr.R",
  file.path(script_dir, "theme_ohchibi_pubr.R")
)
theme_path <- theme_candidates[file.exists(theme_candidates)][1]
if (!is.na(theme_path) && nzchar(theme_path)) {
  source(theme_path)
} else {
  theme_ohchibi_pubr <- ggplot2::theme_bw
}

read_count_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Matrix TSV not found: ", path)
  }
  fs <- file.info(path)$size
  if (is.na(fs) || fs == 0L) {
    stop(
      "Matrix TSV is empty (0 bytes): ", path,
      ". Fix upstream CREATE_NR_NV_MATRICES (unfiltered NR gate) or re-stage inputs."
    )
  }
  dt <- fread(path, sep = "\t", header = TRUE, na.strings = c("", "NA", "."), data.table = TRUE)
  if (is.null(dt) || !is.data.table(dt)) {
    stop("fread failed for matrix TSV: ", path)
  }
  if (ncol(dt) < 2) stop("Matrix TSV must have variant column + >=1 sample column: ", path)
  vid <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  list(vid = vid, mat = mat, samples = colnames(dt)[-1])
}

depth_funnel_tables <- function(nr_mat, nv_mat, samples, depths, min_pct, scheme_label) {
  n_samp <- ncol(nr_mat)
  cohort_rows <- data.table()
  per_sc_rows <- data.table()
  for (d in depths) {
    prop_ge <- if (n_samp > 0) rowSums(nr_mat >= d) / n_samp else numeric(0)
    row_ok <- prop_ge >= (min_pct / 100)
    nv_sub <- nv_mat[row_ok, , drop = FALSE]
    n_cohort <- if (nrow(nv_sub) == 0) 0L else sum(rowSums(nv_sub > 0) >= 1)
    cohort_rows <- rbind(cohort_rows, data.table(
      scheme = scheme_label,
      min_depth = d,
      cohort_variants_with_nv = as.integer(n_cohort)
    ))
    if (n_samp > 0) {
      n_per_sc <- if (nrow(nv_sub) == 0) integer(n_samp) else colSums(nv_sub > 0)
      per_sc_rows <- rbind(per_sc_rows, data.table(
        scheme = scheme_label,
        min_depth = d,
        sample = samples,
        n_variants_nv_gt0 = as.integer(n_per_sc)
      ))
    }
  }
  list(cohort = cohort_rows, per_sc = per_sc_rows)
}

## Focal scheme defines max VAF and which cell(s) tie for it; x = max NR at those cells across all schemes.
depth_max_nr_across_schemes_at_focal_max_vaf <- function(nr_mats, focal_nv_mat, focal_nr_mat) {
  vaf <- focal_nv_mat / pmax(focal_nr_mat, 1)
  vaf[focal_nv_mat <= 0] <- NA
  n_row <- nrow(focal_nv_mat)
  out <- rep(NA_real_, n_row)
  for (i in seq_len(n_row)) {
    rv <- vaf[i, , drop = TRUE]
    if (all(is.na(rv))) next
    mx <- max(rv, na.rm = TRUE)
    if (!is.finite(mx)) next
    pos <- focal_nv_mat[i, ] > 0 & !is.na(rv) & (rv >= mx - 1e-10)
    if (!any(pos)) next
    best_d <- -Inf
    for (j in which(pos)) {
      d_mult <- max(vapply(nr_mats, function(nm) nm[i, j], numeric(1)))
      if (d_mult > best_d) best_d <- d_mult
    }
    if (is.finite(best_d)) out[i] <- best_d
  }
  out
}

## Same tie logic, but NR depth from the focal scheme only (fallback when variant rows differ between schemes).
depth_focal_nr_at_max_vaf <- function(focal_nv_mat, focal_nr_mat) {
  vaf <- focal_nv_mat / pmax(focal_nr_mat, 1)
  vaf[focal_nv_mat <= 0] <- NA
  n_row <- nrow(focal_nv_mat)
  out <- rep(NA_real_, n_row)
  for (i in seq_len(n_row)) {
    rv <- vaf[i, , drop = TRUE]
    if (all(is.na(rv))) next
    mx <- max(rv, na.rm = TRUE)
    if (!is.finite(mx)) next
    pos <- focal_nv_mat[i, ] > 0 & !is.na(rv) & (rv >= mx - 1e-10)
    if (!any(pos)) next
    best_d <- max(focal_nr_mat[i, pos], na.rm = TRUE)
    if (is.finite(best_d)) out[i] <- best_d
  }
  out
}

build_sp_df <- function(nr_mat, nv_mat, nr_mats_all_for_cross = NULL) {
  n_samp <- ncol(nr_mat)
  vaf <- nv_mat / pmax(nr_mat, 1)
  vaf[nv_mat <= 0] <- NA
  mean_vaf <- rowMeans(vaf, na.rm = TRUE)
  max_vaf <- apply(vaf, 1L, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  prop_nv <- if (n_samp > 0) rowSums(nv_mat > 0) / n_samp else numeric(0)
  n_cells_nv_ge_1 <- as.integer(rowSums(nv_mat >= 1))
  dt <- data.table(
    mean_vaf = mean_vaf,
    max_vaf = max_vaf,
    prop_cells_nv_gt0 = prop_nv,
    n_cells_nv_ge_1 = n_cells_nv_ge_1
  )
  if (!is.null(nr_mats_all_for_cross)) {
    dt[, nr_depth_max_across_schemes_at_max_vaf := depth_max_nr_across_schemes_at_focal_max_vaf(
      nr_mats_all_for_cross, nv_mat, nr_mat
    )]
  } else {
    dt[, nr_depth_max_across_schemes_at_max_vaf := depth_focal_nr_at_max_vaf(nv_mat, nr_mat)]
  }
  dt <- dt[
    is.finite(mean_vaf) & is.finite(max_vaf) & is.finite(prop_cells_nv_gt0) &
      is.finite(nr_depth_max_across_schemes_at_max_vaf)
  ]
  dt
}

validation_from_sp <- function(sp_df, scheme_label) {
  n_sp <- nrow(sp_df)
  pct <- function(cond) if (n_sp == 0) NA_real_ else 100 * sum(cond) / n_sp
  data.table(
    scheme = scheme_label,
    metric = c(
      "n_variants_rows",
      "pct_prop_cells_nv_eq_0",
      "pct_prop_cells_nv_gt_0_and_lte_0.05",
      "pct_prop_cells_nv_gt_0.05_and_lte_0.10",
      "pct_prop_cells_nv_gt_0.10",
      "pct_mean_vaf_eq_0",
      "pct_mean_vaf_gt_0_and_lt_0.01",
      "pct_mean_vaf_gte_0.01_and_lt_0.05",
      "pct_mean_vaf_gte_0.05",
      "median_prop_cells_nv_gt0",
      "median_mean_vaf",
      "q90_prop_cells_nv_gt0",
      "q90_mean_vaf",
      "median_max_vaf",
      "q90_max_vaf",
      "interpretation_note"
    ),
    value = c(
      n_sp,
      pct(sp_df$prop_cells_nv_gt0 <= 0),
      pct(sp_df$prop_cells_nv_gt0 > 0 & sp_df$prop_cells_nv_gt0 <= 0.05),
      pct(sp_df$prop_cells_nv_gt0 > 0.05 & sp_df$prop_cells_nv_gt0 <= 0.10),
      pct(sp_df$prop_cells_nv_gt0 > 0.10),
      pct(sp_df$mean_vaf <= 0),
      pct(sp_df$mean_vaf > 0 & sp_df$mean_vaf < 0.01),
      pct(sp_df$mean_vaf >= 0.01 & sp_df$mean_vaf < 0.05),
      pct(sp_df$mean_vaf >= 0.05),
      if (n_sp) as.numeric(median(sp_df$prop_cells_nv_gt0)) else NA_real_,
      if (n_sp) as.numeric(median(sp_df$mean_vaf)) else NA_real_,
      if (n_sp) as.numeric(quantile(sp_df$prop_cells_nv_gt0, 0.90)) else NA_real_,
      if (n_sp) as.numeric(quantile(sp_df$mean_vaf, 0.90)) else NA_real_,
      if (n_sp) as.numeric(median(sp_df$max_vaf)) else NA_real_,
      if (n_sp) as.numeric(quantile(sp_df$max_vaf, 0.90)) else NA_real_,
      NA_real_
    )
  )
}

cohort_all <- data.table()
per_sc_all <- data.table()
sp_all_list <- list()
validation_all <- data.table()

nr_mats <- vector("list", length(scheme_defs))
nv_mats <- vector("list", length(scheme_defs))
vid_ref <- NULL
samples_ref <- NULL
schemes_aligned_for_cross_nr <- TRUE

for (ki in seq_along(scheme_defs)) {
  def <- scheme_defs[[ki]]
  lab <- def$label
  nr <- read_count_matrix(def$nr)
  nv <- read_count_matrix(def$nv)
  if (!identical(nr$vid, nv$vid)) stop("NR/NV variant IDs differ — ", lab)
  if (!identical(nr$samples, nv$samples)) stop("NR/NV sample columns differ — ", lab)
  if (ki == 1L) {
    vid_ref <- nr$vid
    samples_ref <- nr$samples
  } else {
    if (!identical(nr$vid, vid_ref) || !identical(nr$samples, samples_ref)) {
      schemes_aligned_for_cross_nr <- FALSE
    }
  }
  if (nrow(nr$mat) == 0) {
    warning("[analyze_nr_nv_pileup_depth] Empty matrix for scheme ", lab)
  }
  nr_mats[[ki]] <- nr$mat
  nv_mats[[ki]] <- nv$mat

  ft <- depth_funnel_tables(nr$mat, nv$mat, nr$samples, depths, min_sample_pct, lab)
  cohort_all <- rbind(cohort_all, ft$cohort)
  per_sc_all <- rbind(per_sc_all, ft$per_sc)
}

if (!schemes_aligned_for_cross_nr) {
  warning(
    "[analyze_nr_nv_pileup_depth] Variant IDs or sample columns differ between schemes; ",
    "hexbin sheet 3 x-axis uses focal-scheme NR depth at max-VAF cell only (not max NR across schemes)."
  )
}

for (ki in seq_along(scheme_defs)) {
  lab <- scheme_defs[[ki]]$label
  sp_df <- build_sp_df(
    nr_mats[[ki]],
    nv_mats[[ki]],
    nr_mats_all_for_cross = if (schemes_aligned_for_cross_nr) nr_mats else NULL
  )
  sp_df[, scheme := lab]
  sp_all_list[[length(sp_all_list) + 1L]] <- sp_df

  val <- validation_from_sp(sp_df[, .(mean_vaf, max_vaf, prop_cells_nv_gt0)], lab)
  val[, string_value := NA_character_]
  val[metric == "interpretation_note", string_value := paste0(
    "mean_vaf / max_vaf use NV/NR only where NV>0. max_vaf highlights the strongest single-sample ",
    "allele fraction; mean_vaf averages across positive cells. y-axis is always proportion of cells with NV>0."
  )]
  validation_all <- rbind(validation_all, val)
}

sp_all <- rbindlist(sp_all_list, use.names = TRUE)
scheme_lev <- vapply(scheme_defs, function(z) z$label, character(1))
sp_all[, scheme := factor(as.character(scheme), levels = scheme_lev)]

tsv_cohort <- paste0("pileup_depth_cohort_counts_", group_id, ".tsv")
tsv_per_sc <- paste0("pileup_depth_per_sc_", group_id, ".tsv")
fwrite(cohort_all, tsv_cohort, sep = "\t")
fwrite(per_sc_all, tsv_per_sc, sep = "\t")

per_sc_summary <- if (nrow(per_sc_all) > 0) {
  per_sc_all[, .(
    median_n_variants = as.numeric(median(n_variants_nv_gt0)),
    q1_n_variants = as.numeric(quantile(n_variants_nv_gt0, 0.25)),
    q3_n_variants = as.numeric(quantile(n_variants_nv_gt0, 0.75)),
    mean_n_variants = as.numeric(mean(n_variants_nv_gt0)),
    n_single_cells = .N
  ), by = .(scheme, min_depth)][order(scheme, min_depth)]
} else {
  data.table(
    scheme = character(), min_depth = integer(), median_n_variants = numeric(),
    q1_n_variants = numeric(), q3_n_variants = numeric(), mean_n_variants = numeric(),
    n_single_cells = integer()
  )
}
tsv_per_sc_med <- paste0("pileup_depth_per_sc_median_", group_id, ".tsv")
fwrite(per_sc_summary, tsv_per_sc_med, sep = "\t")

tsv_sp <- paste0("pileup_sparsity_vaf_", group_id, ".tsv")
fwrite(
  sp_all[, .(
    scheme,
    mean_vaf,
    max_vaf,
    prop_cells_nv_gt0,
    n_cells_nv_ge_1,
    nr_depth_max_across_schemes_at_max_vaf
  )],
  tsv_sp,
  sep = "\t"
)

tsv_val <- paste0("pileup_sparsity_validation_", group_id, ".tsv")
fwrite(validation_all[, .(scheme, metric, value, string_value)], tsv_val, sep = "\t", na = "")

# ── Depth PDF — REMOVED ──────────────────────────────────────────────────────
# pileup_depth_analysis_<group>.pdf is no longer generated. The cohort / per-SC
# depth TSVs (pileup_depth_cohort_counts, _per_sc, _per_sc_median) are still
# written above; only the depth figure was dropped.

# ── Sparsity: hexbin PDF — sheets 1–2 VAF prevalence; sheet 3 depth at max VAF; sheet 4 NV≥1 cell count vs max VAF ──
use_hex <- requireNamespace("hexbin", quietly = TRUE)
kovesi_cols <- tryCatch(pals::kovesi.rainbow(256), error = function(e) NULL)
hex_legend_counts <- function(x) label_comma(accuracy = 1)(pmax(1, round(10^as.numeric(x))))

hex_geom <- function() {
  if (use_hex) {
    geom_hex(aes(fill = after_stat(log10(pmax(count, 1)))), bins = 45)
  } else {
    geom_bin2d(aes(fill = after_stat(log10(pmax(count, 1)))), bins = 45)
  }
}

hex_scale_fixed <- function(max_log_fill) {
  if (!is.finite(max_log_fill) || max_log_fill <= 0) {
    max_log_fill <- 0.3
  }
  if (!is.null(kovesi_cols) && length(kovesi_cols) > 1) {
    scale_fill_gradientn(
      colours = kovesi_cols,
      name = "Variants per hex",
      labels = hex_legend_counts,
      limits = c(0, max_log_fill),
      oob = scales::squish,
      na.value = "grey92"
    )
  } else {
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#2166AC", "#08306B"),
      name = "Variants per hex",
      labels = hex_legend_counts,
      limits = c(0, max_log_fill),
      oob = scales::squish,
      na.value = "grey92"
    )
  }
}

# Max log10(count) across stat layers (shared fill cap per composite figure)
hex_max_log_from_plot <- function(p) {
  gb <- ggplot_build(p)
  mx <- 0
  for (dt in gb$data) {
    if (!is.null(dt) && "count" %in% names(dt)) {
      v <- suppressWarnings(max(dt$count, na.rm = TRUE))
      if (is.finite(v) && v > 0) {
        mx <- max(mx, log10(max(v, 1)))
      }
    }
  }
  mx
}

br01 <- seq(0, 1, 0.1)
scale_hex_xy <- function() {
  list(
    scale_x_continuous(breaks = br01, limits = c(0, 1), expand = c(0, 0)),
    scale_y_continuous(breaks = br01, limits = c(0, 1), expand = c(0, 0))
  )
}

pdf_hex <- paste0("pileup_sparsity_hexbin_", group_id, ".pdf")
pdf(pdf_hex, width = 12, height = 9)

if (nrow(sp_all) < 2) {
  plot.new()
  text(0.5, 0.5, "Insufficient variants for hexbin", cex = 1.1)
} else {
  pr_mean <- factor("Mean VAF vs prevalence", levels = c("Mean VAF vs prevalence", "Max VAF vs prevalence"))
  pr_max <- factor("Max VAF vs prevalence", levels = c("Mean VAF vs prevalence", "Max VAF vs prevalence"))
  long_prop <- rbind(
    sp_all[, .(scheme, x = mean_vaf, y = prop_cells_nv_gt0, panel_row = pr_mean)],
    sp_all[, .(scheme, x = max_vaf, y = prop_cells_nv_gt0, panel_row = pr_max)]
  )

  # Sheet 1: rows = mean vs prop, max vs prop; cols = scheme order (unfiltered, pileup, …)
  p1_pre <- ggplot(long_prop, aes(x = x, y = y)) +
    hex_geom() +
    facet_grid(panel_row ~ scheme)
  mx1 <- hex_max_log_from_plot(p1_pre)

  p_sheet1 <- ggplot(long_prop, aes(x = x, y = y)) +
    hex_geom() +
    facet_grid(panel_row ~ scheme) +
    scale_hex_xy() +
    hex_scale_fixed(mx1) +
    labs(
      title = paste0("VAF vs cellular prevalence — ", group_id),
      subtitle = paste0(
        "Columns: ", paste(levels(sp_all$scheme), collapse = ", "),
        ". Top: mean VAF (NV>0); bottom: max VAF (NV>0). One color scale (log10 variants per hex)."
      ),
      x = "VAF (NV > 0 cells only)",
      y = "Proportion of cells with NV > 0"
    ) +
    theme_ohchibi_pubr() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1.2, "cm"),
      strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 0.5)
    ) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = unit(12, "cm")))
  print(p_sheet1)

  # Sheet 2: mean vs max, one row × schemes
  p2_pre <- ggplot(sp_all, aes(mean_vaf, max_vaf)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    hex_geom() +
    facet_wrap(~scheme, nrow = 1L, ncol = length(scheme_lev))
  mx2 <- hex_max_log_from_plot(p2_pre)

  p_sheet2 <- ggplot(sp_all, aes(mean_vaf, max_vaf)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    hex_geom() +
    facet_wrap(~scheme, nrow = 1L, ncol = length(scheme_lev)) +
    scale_hex_xy() +
    hex_scale_fixed(mx2) +
    labs(
      title = paste0("Mean vs max VAF (NV>0 cells) — ", group_id),
      subtitle = "Dashed line y = x. One color scale (log10 variants per hex).",
      x = "Mean VAF (NV>0 cells only)",
      y = "Max VAF (NV>0 cells only)"
    ) +
    theme_ohchibi_pubr() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1.2, "cm")
    ) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = unit(12, "cm")))
  print(p_sheet2)

  # Sheet 3: max VAF (focal scheme) vs max NR depth at cell(s) tied for max VAF (NR maxed across schemes)
  p3_pre <- ggplot(sp_all, aes(x = nr_depth_max_across_schemes_at_max_vaf, y = max_vaf)) +
    hex_geom() +
    facet_wrap(~scheme, nrow = 1L, ncol = length(scheme_lev))
  mx3 <- hex_max_log_from_plot(p3_pre)

  sub3_cross <- paste0(
    "Per column (scheme): y = max VAF in that scheme; among NV>0 cells tied for that max, ",
    "x = max NR at those cells, taking the maximum across schemes (",
    paste(scheme_lev, collapse = ", "), "). One color scale (log10 variants per hex)."
  )
  sub3_focal <- paste0(
    "Per column (scheme): y = max VAF; x = NR depth in that scheme at cell(s) tied for max VAF ",
    "(variant rows or samples differ between schemes, so cross-scheme NR max is not applied). ",
    "One color scale (log10 variants per hex)."
  )
  p_sheet3 <- ggplot(sp_all, aes(x = nr_depth_max_across_schemes_at_max_vaf, y = max_vaf)) +
    hex_geom() +
    facet_wrap(~scheme, nrow = 1L, ncol = length(scheme_lev)) +
    scale_x_continuous(labels = label_comma(), expand = expansion(mult = c(0.02, 0.04))) +
    scale_y_continuous(breaks = br01, limits = c(0, 1), expand = c(0, 0)) +
    hex_scale_fixed(mx3) +
    labs(
      title = paste0("Max VAF vs NR depth at max-VAF cell — ", group_id),
      subtitle = if (schemes_aligned_for_cross_nr) sub3_cross else sub3_focal,
      x = if (schemes_aligned_for_cross_nr) {
        "NR depth (max across schemes at max-VAF cell)"
      } else {
        "NR depth (focal scheme at max-VAF cell)"
      },
      y = "Max VAF (NV>0 cells, focal scheme)"
    ) +
    theme_ohchibi_pubr() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1.2, "cm")
    ) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = unit(12, "cm")))
  print(p_sheet3)

  # Sheet 4: max VAF vs number of cells with NV ≥ 1 (one row × schemes)
  p4_pre <- ggplot(sp_all, aes(x = max_vaf, y = n_cells_nv_ge_1)) +
    hex_geom() +
    facet_wrap(~scheme, nrow = 1L, ncol = length(scheme_lev))
  mx4 <- hex_max_log_from_plot(p4_pre)

  p_sheet4 <- ggplot(sp_all, aes(x = max_vaf, y = n_cells_nv_ge_1)) +
    hex_geom() +
    facet_wrap(~scheme, nrow = 1L, ncol = length(scheme_lev)) +
    scale_x_continuous(breaks = br01, limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.02, 0.04))) +
    hex_scale_fixed(mx4) +
    labs(
      title = paste0("Cells with NV ≥ 1 vs max VAF — ", group_id),
      subtitle = paste0(
        "Per column (scheme): x = max VAF (NV>0 cells); y = count of cells with NV ≥ 1 (any depth). ",
        "One color scale (log10 variants per hex)."
      ),
      x = "Max VAF (NV > 0 cells only)",
      y = "Number of cells with NV ≥ 1"
    ) +
    theme_ohchibi_pubr() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1.2, "cm")
    ) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = unit(12, "cm")))
  print(p_sheet4)
}

dev.off()

message(
  "[analyze_nr_nv_pileup_depth] schemes: ", paste(vapply(scheme_defs, function(z) z$label, ""), collapse = ", "),
  " | Wrote: ", tsv_cohort, ", ", tsv_per_sc, ", ", tsv_per_sc_med, ", ",
  tsv_sp, ", ", tsv_val, ", ", pdf_hex
)
