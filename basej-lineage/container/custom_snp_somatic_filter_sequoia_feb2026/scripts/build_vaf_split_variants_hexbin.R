#!/usr/bin/env Rscript
# build_vaf_split_variants_hexbin.R
# ---------------------------------------------------------------------------
# Reproducible pipeline that, given an NV count matrix, an NR count matrix,
# and a variant master filter table, will:
#
#   1. Build a VAF matrix as NV / NR, mapping 0/0 entries to VAF = 0.
#   2. Subset the master table to variants present in the NR matrix.
#   3. From the NV matrix, classify every variant as singleton (present in
#      exactly one cell, NV > 0) or shared (present in more than one cell)
#      and write the two corresponding ID lists.
#   4. Build a hexbin plot of SecondRunSequoia_Rho (x) vs
#      log10(SecondRunSequoia_Germline_qval) (y), using *only* the shared
#      variants, with red longdash threshold lines at Rho = 0.2 and
#      log10(qval) = -10.
#   5. Retain the shared variants that pass both thresholds
#      (SecondRunSequoia_Rho >= rho_thr  AND
#       log10(SecondRunSequoia_Germline_qval) <= log10q_thr) and write the
#      retained ID list plus the matching master-table subset.
#
#   Downstream, the filtered (retained-shared + singleton) matrices are explored
#   and binarized SEPARATELY for the two variant classes:
#     * Threshold exploration (step 7) emits one figure for shared variants and
#       one for singletons.
#     * Binarization (step 8) uses class-specific NV/VAF rules:
#         - singletons : NV >= --bin-nv        & VAF >= --bin-vaf-singleton
#         - shared     : anchor cell = NV >= --bin-nv & VAF >= --bin-vaf-shared-anchor;
#                        if >=1 anchor cell exists, the variant's other cells are
#                        called at the relaxed gate NV >= --bin-nv-shared &
#                        VAF >= --bin-vaf; shared variants with no anchor cell
#                        require the strict anchor gate in every cell.
#       (The NR depth gate --bin-nr applies to every rule.)
#
# Inputs are TSVs:
#   --nv      : NV count matrix (first column = VariantId, remaining = cells)
#   --nr      : NR count matrix (same shape and order as NV)
#   --master  : variant master filter table (must have a 'VariantId' column,
#               plus 'SecondRunSequoia_Rho' and 'SecondRunSequoia_Germline_qval'
#               for the hexbin plot)
#   --outdir  : output directory (created if missing)
#   --prefix  : optional file-name prefix for outputs
#
# Outputs (under --outdir):
#   <prefix>VAF_matrix.tsv
#   <prefix>variant_master_filter_table_in_NR.tsv
#   <prefix>singleton_variants.txt
#   <prefix>shared_variants.txt
#   <prefix>hexbin_SecondRunSequoia_Rho_vs_GermlineQval_shared.pdf
#   <prefix>hexbin_SecondRunSequoia_Rho_vs_GermlineQval_shared.png
#   <prefix>shared_variants_retained.txt
#   <prefix>variant_master_filter_table_shared_retained.tsv
#
# Usage:
#   Rscript build_vaf_split_variants_hexbin.R \
#     --nv     NV_annotated_vcf_Rodrigues_unfiltered.tsv \
#     --nr     NR_annotated_vcf_Rodrigues_unfiltered.tsv \
#     --master variant_master_filter_table_Rodrigues.tsv \
#     --outdir /home/ubuntu/projects/rodrigues_artefacts \
#     --prefix Rodrigues_

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(hexbin)
  library(scales)
  library(egg)
  library(optparse)
})

# Source the shared ggplot theme. Baked into the container at /usr/local/bin;
# fall back to the script directory, then to theme_bw if neither is found.
argv0 <- commandArgs(trailingOnly = FALSE)
file_arg <- argv0[grepl("^--file=", argv0)]
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg[1])) else "."
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

# ---- CLI --------------------------------------------------------------------
option_list <- list(
  make_option("--nv",     type = "character", help = "NV matrix TSV (VariantId + cells)"),
  make_option("--nr",     type = "character", help = "NR matrix TSV (VariantId + cells)"),
  make_option("--master", type = "character", help = "Variant master filter table TSV"),
  make_option("--outdir", type = "character", default = ".", help = "Output directory [default: %default]"),
  make_option("--prefix", type = "character", default = "SHFQC_pipeline_",  help = "Output file-name prefix [default: %default]"),
  make_option("--rho-col",  type = "character", default = "SecondRunSequoia_Rho",
              help = "Master-table column for Rho [default: %default]"),
  make_option("--qval-col", type = "character", default = "SecondRunSequoia_Germline_qval",
              help = "Master-table column for germline qval [default: %default]"),
  make_option("--rho-thr",     type = "double", default = 0.2,
              help = "Rho threshold drawn as red longdash [default: %default]"),
  make_option("--log10q-thr",  type = "double", default = -10,
              help = "log10(qval) threshold drawn as red longdash [default: %default]"),
  make_option("--first-rho-col",  type = "character", default = "Binom_Rho",
              help = "Master-table column for first-round (binom/betabinom) Rho [default: %default]"),
  make_option("--first-qval-col", type = "character", default = "Binom_Germline_qval",
              help = "Master-table column for first-round germline qval [default: %default]"),
  make_option("--first-rho-thr",     type = "double", default = 0.2,
              help = "First-round Rho threshold (red longdash) [default: %default]"),
  make_option("--first-log10q-thr",  type = "double", default = -5,
              help = "First-round log10(qval) threshold (red longdash) [default: %default]"),
  make_option("--y-lo",        type = "double", default = -50,
              help = "Lower limit of the y-axis (log10 qval) [default: %default]"),
  make_option("--y-hi",        type = "double", default = 0,
              help = "Upper limit of the y-axis (log10 qval) [default: %default]"),
  make_option("--y-step",      type = "double", default = 5,
              help = "Step (in log10 units) for y-axis breaks [default: %default]"),
  make_option("--bin-nv",      type = "integer", default = 3L,
              help = paste("Anchor/strict minimum NV (>=): used for singleton calls and for",
                           "shared anchor detection / unanchored shared cells [default: %default]")),
  make_option("--bin-nv-shared", type = "integer", default = 2L,
              help = paste("Relaxed minimum NV (>=) for the companion cells of an anchored",
                           "shared variant [default: %default]")),
  make_option("--bin-nr",      type = "integer", default = 10L,
              help = "Binarization threshold: minimum NR (>=) [default: %default]"),
  make_option("--bin-vaf",     type = "double",  default = 0.10,
              help = paste("Shared-variant relaxed VAF threshold (>=): applied to the other",
                           "cells of a shared variant once it has an anchor cell",
                           "[default: %default]")),
  make_option("--bin-vaf-shared-anchor", type = "double", default = 0.30,
              help = paste("Shared-variant anchor VAF threshold (>=): a shared variant needs",
                           ">=1 cell at this VAF to unlock the relaxed threshold for its other",
                           "cells; shared variants with no anchor cell require this threshold",
                           "in every cell [default: %default]")),
  make_option("--bin-vaf-singleton", type = "double", default = 0.30,
              help = "Singleton-variant VAF threshold (>=) [default: %default]"),
  make_option("--mandatory",   type = "character", default = NULL,
              help = paste("Optional TSV with mandatory variants (first column = VariantId,",
                           "e.g. merged_priority_variants_*.tsv). Variants not already in the",
                           "QCFiltered / ForPhylogeny cascades are appended into",
                           "*PlusMandatory.tsv outputs; binarization for the appended rows is",
                           "NV > 0 -> 1 [default: %default]"))
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("nv", "nr", "master")
miss <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(miss) > 0) {
  stop("Missing required argument(s): ", paste(paste0("--", miss), collapse = ", "))
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(name) file.path(opt$outdir, paste0(opt$prefix, name))

out_vaf       <- out("VAF_matrix.tsv")
out_master    <- out("variant_master_filter_table_in_NR.tsv")
out_single    <- out("singleton_variants.txt")
out_shared    <- out("shared_variants.txt")
out_pdf       <- out("hexbin_SecondRunSequoia_Rho_vs_GermlineQval_shared.pdf")
out_png       <- out("hexbin_SecondRunSequoia_Rho_vs_GermlineQval_shared.png")
out_first_pdf <- out("hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.pdf")
out_first_png <- out("hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png")
out_retained        <- out("shared_variants_retained.txt")
out_master_retained <- out("variant_master_filter_table_shared_retained.tsv")
out_final_nv  <- out("NV_HQRoundStatisticalFiltered.tsv")
out_final_nr  <- out("NR_HQRoundStatisticalFiltered.tsv")
out_final_vaf <- out("VAF_HQRoundStatisticalFiltered.tsv")

cat("Inputs:\n")
cat(sprintf("  NV     : %s\n", opt$nv))
cat(sprintf("  NR     : %s\n", opt$nr))
cat(sprintf("  master : %s\n", opt$master))
cat(sprintf("Outputs in: %s (prefix='%s')\n", opt$outdir, opt$prefix))

# ---- 1. Load NV / NR --------------------------------------------------------
cat("\n[1/4] Loading NV / NR matrices...\n")
nv_dt <- fread(opt$nv, sep = "\t", header = TRUE, check.names = FALSE)
nr_dt <- fread(opt$nr, sep = "\t", header = TRUE, check.names = FALSE)
setnames(nv_dt, 1, "VariantId")
setnames(nr_dt, 1, "VariantId")

if (!identical(nv_dt$VariantId, nr_dt$VariantId)) {
  stop("VariantId order differs between NV and NR matrices. Aborting.")
}
if (!identical(colnames(nv_dt), colnames(nr_dt))) {
  stop("Cell column order differs between NV and NR matrices. Aborting.")
}

variant_ids <- nv_dt$VariantId
cell_cols   <- colnames(nv_dt)[-1]
cat(sprintf("  %d variants x %d cells\n", length(variant_ids), length(cell_cols)))

nv <- as.matrix(nv_dt[, -1]); mode(nv) <- "numeric"
nr <- as.matrix(nr_dt[, -1]); mode(nr) <- "numeric"

# ---- Drop variants with NV==0 in every cell (no variant-read evidence) -------
# Applied at the source so ALL outputs (unfiltered VAF/binary matrices, the
# master-in-NR subset, and the singleton/shared cascade) operate on evidence-
# bearing variants only. Keeps the two matrix lineages (VAF_SPLIT group-level and
# CREATE_NR_NV per-sample) consistent. Cascade stage counts are unchanged, since
# all-NV-zero variants are NV>0 in 0 cells and were already excluded downstream.
nv_evidence <- rowSums(nv > 0) >= 1
n_drop_nv0  <- sum(!nv_evidence)
if (n_drop_nv0 > 0) {
  cat(sprintf("  Dropping %d variant(s) with NV==0 in every cell (of %d); %d retained.\n",
              n_drop_nv0, length(variant_ids), sum(nv_evidence)))
  nv_dt       <- nv_dt[nv_evidence]
  nr_dt       <- nr_dt[nv_evidence]
  nv          <- nv[nv_evidence, , drop = FALSE]
  nr          <- nr[nv_evidence, , drop = FALSE]
  variant_ids <- variant_ids[nv_evidence]
}

# ---- Build VAF matrix (NV / NR; 0/0 -> 0) -----------------------------------
cat("\n[1/4] Building VAF matrix (NV / NR; 0/0 -> 0)...\n")
vaf <- matrix(0.0, nrow = nrow(nv), ncol = ncol(nv),
              dimnames = list(NULL, cell_cols))
covered <- nr > 0
vaf[covered] <- nv[covered] / nr[covered]

n_capped <- sum(vaf > 1, na.rm = TRUE)
if (n_capped > 0) {
  cat(sprintf("  Note: %d entries had NV>NR; capping VAF at 1.0\n", n_capped))
  vaf[vaf > 1] <- 1
}

vaf_out <- data.table(VariantId = variant_ids, as.data.table(round(vaf, 6)))
fwrite(vaf_out, out_vaf, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s\n", out_vaf))

# Binary matrix for unfiltered data: presence = NV > 0
out_binary_unfiltered <- out("binary_matrix_unfiltered.tsv")
binary_unfilt <- as.integer(nv > 0)
dim(binary_unfilt) <- dim(nv)
colnames(binary_unfilt) <- cell_cols
binary_unfilt_dt <- data.table(VariantId = variant_ids, as.data.table(binary_unfilt))
fwrite(binary_unfilt_dt, out_binary_unfiltered, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_binary_unfiltered, nrow(binary_unfilt_dt), length(cell_cols)))

# ---- 2. Subset master table to NR variants ---------------------------------
cat("\n[2/4] Subsetting master table to variants present in NR matrix...\n")
master <- fread(opt$master, sep = "\t", header = TRUE, check.names = FALSE, quote = "")
if (!"VariantId" %in% colnames(master)) {
  stop("Master table does not have a 'VariantId' column.")
}
cat(sprintf("  master: %d rows x %d cols\n", nrow(master), ncol(master)))

master_sub <- master[VariantId %in% variant_ids]
n_unique_kept    <- length(unique(master_sub$VariantId))
n_nr_in_master   <- sum(variant_ids %in% master$VariantId)
n_nr_missing     <- length(variant_ids) - n_nr_in_master
cat(sprintf("  master rows kept              : %d\n", nrow(master_sub)))
cat(sprintf("  unique master VariantIds kept : %d\n", n_unique_kept))
cat(sprintf("  NR variants found in master   : %d / %d\n",
            n_nr_in_master, length(variant_ids)))
cat(sprintf("  NR variants NOT in master     : %d\n", n_nr_missing))

fwrite(master_sub, out_master, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Wrote %s\n", out_master))

# ---- 3. Singleton vs shared variants from NV matrix -------------------------
cat("\n[3/4] Classifying variants as singleton vs shared (NV > 0 across cells)...\n")
n_pos_per_var <- rowSums(nv > 0)
is_singleton  <- n_pos_per_var == 1
is_shared     <- n_pos_per_var >  1
n_zero_cells  <- sum(n_pos_per_var == 0)

singleton_ids <- variant_ids[is_singleton]
shared_ids    <- variant_ids[is_shared]

cat(sprintf("  variants with NV>0 in 0 cells  : %d (excluded from both lists)\n", n_zero_cells))
cat(sprintf("  singletons (exactly 1 cell)    : %d\n", length(singleton_ids)))
cat(sprintf("  shared    (more than 1 cell)   : %d\n", length(shared_ids)))

writeLines(singleton_ids, out_single)
writeLines(shared_ids,    out_shared)
cat(sprintf("  Wrote %s\n", out_single))
cat(sprintf("  Wrote %s\n", out_shared))

# ---- 4. Hexbin plots: Rho vs germline q-value -------------------------------
# Reusable builder, called twice:
#   - round 2 (second-pass Sequoia): SecondRunSequoia_* columns, SHARED variants
#   - round 1 (first-pass binom/betabinom): Binom_* columns, ALL variants
cat("\n[4/11] Building Rho-vs-germline-qval hexbin plots...\n")

y_lo <- opt[["y-lo"]]; y_hi <- opt[["y-hi"]]; y_step <- opt[["y-step"]]
if (y_lo >= y_hi)  stop("--y-lo must be < --y-hi")
if (y_step <= 0)   stop("--y-step must be > 0")
y_breaks     <- seq(y_hi, y_lo, by = -y_step)
rainbow_cols <- c("#4B0082", "#0000FF", "#00FFFF", "#00FF00",
                  "#FFFF00", "#FFA500", "#FF0000")

# Hexbin of <rho_col> (x, restricted to [0,1]) vs log10(<qval_col>) (y), with red
# longdash threshold lines at rho_thr / log10q_thr. Skips gracefully if columns
# are missing or no variant has both values.
make_rho_qval_hexbin <- function(dt_in, rho_col, qval_col, rho_thr, log10q_thr,
                                 title_main, subset_label, out_pdf, out_png) {
  if (!rho_col %in% colnames(dt_in) || !qval_col %in% colnames(dt_in)) {
    cat(sprintf("  [%s] column(s) '%s' / '%s' not in master table; skipping.\n",
                subset_label, rho_col, qval_col)); return(invisible(NULL))
  }
  dt <- dt_in[, c("VariantId", rho_col, qval_col), with = FALSE]
  setnames(dt, c("VariantId", "Rho", "qval"))
  dt[, Rho  := suppressWarnings(as.numeric(Rho))]
  dt[, qval := suppressWarnings(as.numeric(qval))]
  cat(sprintf("  [%s] variants: %d ; with both Rho and qval: %d\n",
              subset_label, nrow(dt), sum(!is.na(dt$Rho) & !is.na(dt$qval))))

  dt_plot <- dt[!is.na(Rho) & !is.na(qval)]
  floor_q <- 1e-300
  dt_plot[qval <= 0,       qval := floor_q]
  dt_plot[qval <  floor_q, qval := floor_q]
  dt_plot[, log10_qval := log10(qval)]
  oor <- dt_plot[Rho < 0 | Rho > 1, .N]
  if (oor > 0) cat(sprintf("  [%s] Rho outside [0,1] (excluded): %d\n", subset_label, oor))
  dt_plot <- dt_plot[Rho >= 0 & Rho <= 1]
  if (nrow(dt_plot) == 0) {
    warning(sprintf("[%s] no variants with both Rho and qval; skipping hexbin.", subset_label))
    return(invisible(NULL))
  }

  p <- ggplot(dt_plot, aes(x = Rho, y = log10_qval)) +
    geom_hex(bins = 50) +
    scale_fill_gradientn("Variants\n(log scale)", trans = "log10", colours = rainbow_cols) +
    scale_x_continuous(name = rho_col, breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    scale_y_continuous(name = bquote(log[10](.(qval_col))), limits = c(y_lo, y_hi),
                       breaks = y_breaks, oob = scales::squish) +
    geom_vline(xintercept = rho_thr,    linetype = "longdash", colour = "red", linewidth = 0.5) +
    geom_hline(yintercept = log10q_thr, linetype = "longdash", colour = "red", linewidth = 0.5) +
    annotate("text", x = rho_thr, y = y_hi, label = sprintf("Rho = %g", rho_thr),
             hjust = -0.05, vjust = 1.3, colour = "red", size = 3) +
    annotate("text", x = 1, y = log10q_thr, label = sprintf("log10(qval) = %g", log10q_thr),
             hjust = 1.05, vjust = -0.4, colour = "red", size = 3) +
    ggtitle(title_main, subtitle = sprintf("n = %s %s with both values",
                                            format(nrow(dt_plot), big.mark = ","), subset_label)) +
    theme_ohchibi_pubr()

  cairo_pdf(out_pdf, width = 7, height = 5); print(p); dev.off()
  ggsave(out_png, p, width = 7, height = 5, dpi = 200)
  cat(sprintf("  [%s] Wrote %s\n", subset_label, out_pdf))
  cat(sprintf("  [%s] Wrote %s\n", subset_label, out_png))
  invisible(NULL)
}

# Round 2 — second-pass Sequoia, shared variants only (the original figure).
make_rho_qval_hexbin(
  master_sub[VariantId %in% shared_ids],
  opt[["rho-col"]], opt[["qval-col"]], opt[["rho-thr"]], opt[["log10q-thr"]],
  "Sequoia second-run: Rho vs germline q-value (shared variants)", "shared variants",
  out_pdf, out_png)

# Round 1 — first-pass binom/betabinom, all variants in the NR matrix.
make_rho_qval_hexbin(
  master_sub,
  opt[["first-rho-col"]], opt[["first-qval-col"]], opt[["first-rho-thr"]], opt[["first-log10q-thr"]],
  "Sequoia first-run: Rho vs germline q-value (all variants)", "variants",
  out_first_pdf, out_first_png)

# ---- 5. Retain shared variants that pass the thresholds --------------------
# Rebuild the shared-variant second-pass Rho/qval table that step 5 consumes
# (the hexbin builder above works on a local copy, so reconstruct it here).
rho_thr    <- opt[["rho-thr"]]
log10q_thr <- opt[["log10q-thr"]]
floor_q    <- 1e-300
dt <- master_sub[VariantId %in% shared_ids,
                 c("VariantId", opt[["rho-col"]], opt[["qval-col"]]), with = FALSE]
setnames(dt, c("VariantId", "Rho", "qval"))
dt[, Rho  := suppressWarnings(as.numeric(Rho))]
dt[, qval := suppressWarnings(as.numeric(qval))]

cat("\n[5/11] Retaining shared variants with Rho >= ",
    rho_thr, " AND log10(qval) <= ", log10q_thr, "...\n", sep = "")

# Use the full shared-variant table (`dt`), not the plot-clipped one, so we
# don't drop variants just because their Rho fell outside the [0,1] plot range.
dt_filt <- copy(dt)
dt_filt[, qval_filt := qval]
dt_filt[!is.na(qval_filt) & qval_filt <= 0,        qval_filt := floor_q]
dt_filt[!is.na(qval_filt) & qval_filt <  floor_q,  qval_filt := floor_q]
dt_filt[, log10_qval := suppressWarnings(log10(qval_filt))]

n_missing_rho_or_qval <- sum(is.na(dt_filt$Rho) | is.na(dt_filt$log10_qval))
cat(sprintf("  shared variants                       : %d\n", nrow(dt_filt)))
cat(sprintf("  with NA Rho or qval (excluded)        : %d\n", n_missing_rho_or_qval))

retained <- dt_filt[!is.na(Rho) & !is.na(log10_qval) &
                    Rho >= rho_thr & log10_qval <= log10q_thr]

cat(sprintf("  passing Rho >= %g AND log10(qval) <= %g : %d\n",
            rho_thr, log10q_thr, nrow(retained)))

retained_ids <- retained$VariantId
writeLines(retained_ids, out_retained)
cat(sprintf("  Wrote %s\n", out_retained))

master_retained <- master_sub[VariantId %in% retained_ids]
fwrite(master_retained, out_master_retained, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Wrote %s (%d rows)\n", out_master_retained, nrow(master_retained)))

# ---- 6. Union of retained shared + singletons → subset NV, NR, VAF ---------
cat("\n[6/11] Building union of retained shared variants + singletons and subsetting matrices...\n")

union_ids <- unique(c(retained_ids, singleton_ids))
cat(sprintf("  retained shared variants : %d\n", length(retained_ids)))
cat(sprintf("  singleton variants       : %d\n", length(singleton_ids)))
cat(sprintf("  union (unique)           : %d\n", length(union_ids)))

# Row indices in the original matrices
keep_idx <- which(variant_ids %in% union_ids)
cat(sprintf("  rows matched in matrices : %d\n", length(keep_idx)))

# Subset NV
nv_sub <- data.table(VariantId = variant_ids[keep_idx],
                     as.data.table(nv[keep_idx, , drop = FALSE]))
fwrite(nv_sub, out_final_nv, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_final_nv, nrow(nv_sub), length(cell_cols)))

# Subset NR
nr_sub <- data.table(VariantId = variant_ids[keep_idx],
                     as.data.table(nr[keep_idx, , drop = FALSE]))
fwrite(nr_sub, out_final_nr, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_final_nr, nrow(nr_sub), length(cell_cols)))

# Subset VAF
vaf_sub <- data.table(VariantId = variant_ids[keep_idx],
                      as.data.table(round(vaf[keep_idx, , drop = FALSE], 6)))
fwrite(vaf_sub, out_final_vaf, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_final_vaf, nrow(vaf_sub), length(cell_cols)))

# Binary matrix for HQRoundStatisticalFiltered: presence = NV > 0
out_binary_hqstat <- out("binary_matrix_HQRoundStatisticalFiltered.tsv")
nv_sub_mat <- as.matrix(nv_sub[, -1]); mode(nv_sub_mat) <- "numeric"
binary_hqstat <- as.integer(nv_sub_mat > 0)
dim(binary_hqstat) <- dim(nv_sub_mat)
colnames(binary_hqstat) <- cell_cols
binary_hqstat_dt <- data.table(VariantId = nv_sub$VariantId, as.data.table(binary_hqstat))
fwrite(binary_hqstat_dt, out_binary_hqstat, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_binary_hqstat, nrow(binary_hqstat_dt), length(cell_cols)))

# ---- 7. Threshold exploration figures (shared + singleton) ------------------
cat("\n[7/11] Building threshold exploration figures (shared and singleton)...\n")

# Numeric matrices for the filtered subset
nv_f  <- as.matrix(nv_sub[, -1]); mode(nv_f) <- "numeric"
nr_f  <- as.matrix(nr_sub[, -1]); mode(nr_f) <- "numeric"
vaf_f <- as.matrix(vaf_sub[, -1]); mode(vaf_f) <- "numeric"

n_cell_f <- ncol(nv_f)

# Split the filtered rows into shared (retained) vs singleton
row_ids_f       <- nv_sub$VariantId
shared_row_f    <- row_ids_f %in% retained_ids
singleton_row_f <- row_ids_f %in% singleton_ids
cat(sprintf("  filtered rows: %d shared + %d singleton (of %d total)\n",
            sum(shared_row_f), sum(singleton_row_f), length(row_ids_f)))

# ---- Reusable builder: 6-panel threshold exploration for one variant class --
build_threshold_exploration <- function(nv_f, nr_f, vaf_f, tag,
                                        out_thresh_pdf, out_thresh_png) {
  n_var_f  <- nrow(nv_f)
  n_cell_f <- ncol(nv_f)
  if (n_var_f == 0) {
    cat(sprintf("  [%s] no variants in this class; skipping exploration figure.\n", tag))
    return(invisible(NULL))
  }
  cat(sprintf("  [%s] %d variants x %d cells\n", tag, n_var_f, n_cell_f))

  # Flatten
  nr_v  <- as.vector(nr_f)
  nv_v  <- as.vector(nv_f)
  vaf_v <- as.vector(vaf_f)

  dt_all_f <- data.table(NR = nr_v, NV = nv_v, VAF = vaf_v)
  dt_cov_f <- dt_all_f[NR > 0]

  # Down-sample helper for plotting speed
  set.seed(42)
  sample_for_plot <- function(x, n_max = 2e5) {
    if (length(x) <= n_max) return(x)
    x[sample.int(length(x), n_max)]
  }

  # ---- Panel 1: NR distribution (NR > 0) on log10 scale ----
  nr_pos_f <- dt_cov_f$NR
  p_nr_f <- ggplot(data.table(NR = sample_for_plot(nr_pos_f)),
                   aes(x = log10(NR + 1))) +
    geom_histogram(bins = 80, fill = "#3B7DD8", colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = log10(c(5, 10, 20) + 1),
               linetype = c("dashed", "dashed", "dotted"),
               colour   = c("red", "darkorange", "grey40"),
               linewidth = 0.4) +
    annotate("text", x = log10(5  + 1), y = Inf, label = "NR=5",  hjust = -0.05, vjust = 1.3, colour = "red",        size = 3) +
    annotate("text", x = log10(10 + 1), y = Inf, label = "NR=10", hjust = -0.05, vjust = 2.6, colour = "darkorange", size = 3) +
    annotate("text", x = log10(20 + 1), y = Inf, label = "NR=20", hjust = -0.05, vjust = 3.9, colour = "grey40",     size = 3) +
    scale_x_continuous(
      name = expression(log[10](NR + 1)),
      breaks = log10(c(0, 1, 5, 10, 50, 100, 500, 1000) + 1),
      labels = c(0, 1, 5, 10, 50, 100, 500, 1000)
    ) +
    ylab("Cell-variant entries (NR > 0)") +
    ggtitle(sprintf("NR (read depth) per cell-variant entry [%s]", tag),
            subtitle = sprintf("median=%g  \u00b7  p99=%g",
                               median(nr_pos_f), quantile(nr_pos_f, 0.99))) +
    theme_ohchibi_pubr()

  # ---- Panel 2: NV distribution (NV > 0) ----
  nv_pos_f <- dt_all_f[NV > 0]$NV
  p_nv_f <- ggplot(data.table(NV = sample_for_plot(nv_pos_f)),
                   aes(x = log10(NV + 1))) +
    geom_histogram(bins = 80, fill = "#E07B39", colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = log10(c(2, 3, 6) + 1),
               linetype = c("dashed", "dotted", "dotted"),
               colour   = c("red", "grey40", "grey40"),
               linewidth = 0.4) +
    annotate("text", x = log10(2 + 1), y = Inf, label = "NV=2", hjust = -0.05, vjust = 1.3, colour = "red",    size = 3) +
    annotate("text", x = log10(3 + 1), y = Inf, label = "NV=3", hjust = -0.05, vjust = 2.6, colour = "grey40", size = 3) +
    annotate("text", x = log10(6 + 1), y = Inf, label = "NV=6", hjust = -0.05, vjust = 3.9, colour = "grey40", size = 3) +
    scale_x_continuous(
      name = expression(log[10](NV + 1)),
      breaks = log10(c(0, 1, 2, 5, 10, 50, 100, 500) + 1),
      labels = c(0, 1, 2, 5, 10, 50, 100, 500)
    ) +
    ylab("Cell-variant entries (NV > 0)") +
    ggtitle(sprintf("NV (variant reads) per cell-variant entry [%s]", tag),
            subtitle = sprintf("median=%g  \u00b7  p99=%g",
                               median(nv_pos_f), quantile(nv_pos_f, 0.99))) +
    theme_ohchibi_pubr()

  # ---- Panel 3: VAF distribution by NR coverage bin ----
  nr_breaks_f <- c(1, 5, 10, 20, Inf)
  nr_labels_f <- c("NR=1-4", "NR=5-9", "NR=10-19", "NR>=20")
  dt_cov_f[, NR_bin := cut(NR, breaks = nr_breaks_f, right = FALSE,
                            labels = nr_labels_f, include.lowest = TRUE)]
  dt_plot_vaf_f <- dt_cov_f[, .SD[sample(.N, min(.N, 100000))], by = NR_bin]

  p_vaf_f <- ggplot(dt_plot_vaf_f, aes(x = VAF, fill = NR_bin)) +
    geom_histogram(bins = 60, position = "identity", alpha = 0.55,
                   colour = "white", linewidth = 0.1) +
    geom_vline(xintercept = c(0.05, 0.10, 0.20, 0.30),
               linetype = c("dotted", "dashed", "dotted", "dotted"),
               colour   = c("grey40", "red", "grey40", "grey40"),
               linewidth = 0.4) +
    annotate("text", x = 0.05, y = Inf, label = "0.05", hjust = -0.05, vjust = 1.3, colour = "grey40", size = 3) +
    annotate("text", x = 0.10, y = Inf, label = "0.10", hjust = -0.05, vjust = 2.6, colour = "red",    size = 3) +
    annotate("text", x = 0.20, y = Inf, label = "0.20", hjust = -0.05, vjust = 3.9, colour = "grey40", size = 3) +
    annotate("text", x = 0.30, y = Inf, label = "0.30", hjust = -0.05, vjust = 5.2, colour = "grey40", size = 3) +
    scale_x_continuous(name = "VAF (NV / NR, NR > 0)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    ylab("Cell-variant entries (sampled)") +
    scale_fill_brewer("NR bin", palette = "Set2") +
    ggtitle(sprintf("VAF distribution by coverage bin [%s]", tag),
            subtitle = "low-depth entries dominate the noise floor near 0 and 0.5/1.0 spikes") +
    theme_ohchibi_pubr()

  # ---- Panel 4: VAF density for callable entries (NR>=5 & NV>=1) ----
  dt_call_f <- dt_cov_f[NR >= 5 & NV >= 1]
  p_vaf_callable_f <- ggplot(dt_call_f[sample(.N, min(.N, 200000))], aes(x = VAF)) +
    geom_density(fill = "#3B7DD8", alpha = 0.5, colour = "#1F4E8C") +
    geom_vline(xintercept = c(0.05, 0.10, 0.20, 0.30),
               linetype = c("dotted", "dashed", "dotted", "dotted"),
               colour   = c("grey40", "red", "grey40", "grey40"),
               linewidth = 0.4) +
    annotate("text", x = 0.05, y = Inf, label = "0.05", hjust = -0.05, vjust = 1.3, colour = "grey40", size = 3) +
    annotate("text", x = 0.10, y = Inf, label = "0.10", hjust = -0.05, vjust = 2.6, colour = "red",    size = 3) +
    annotate("text", x = 0.20, y = Inf, label = "0.20", hjust = -0.05, vjust = 3.9, colour = "grey40", size = 3) +
    annotate("text", x = 0.30, y = Inf, label = "0.30", hjust = -0.05, vjust = 5.2, colour = "grey40", size = 3) +
    scale_x_continuous(name = "VAF", limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    ylab("Density") +
    ggtitle(sprintf("VAF density among 'callable' entries  (NR >= 5 & NV >= 1) [%s]", tag),
            subtitle = sprintf("n = %s entries", format(nrow(dt_call_f), big.mark = ","))) +
    theme_ohchibi_pubr()

  # ---- Panel 5: NV vs NR 2D density (NV >= 1) ----
  dt_2d_f <- dt_cov_f[NV >= 1]
  dt_2d_s_f <- dt_2d_f[sample(.N, min(.N, 250000))]

  p_2d_f <- ggplot(dt_2d_s_f, aes(x = log10(NR + 1), y = log10(NV + 1))) +
    geom_bin2d(bins = 80) +
    scale_fill_gradient("Entries\n(log scale)", trans = "log10",
                        low = "#FFF5F0", high = "#67000D") +
    geom_abline(slope = 1, intercept = log10(0.10), linetype = "dashed",
                colour = "blue", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = log10(0.20), linetype = "dotted",
                colour = "blue", linewidth = 0.5) +
    geom_vline(xintercept = log10(5 + 1),  linetype = "dashed", colour = "red", linewidth = 0.4) +
    geom_hline(yintercept = log10(2 + 1),  linetype = "dashed", colour = "red", linewidth = 0.4) +
    annotate("text", x = log10(5 + 1),  y = Inf, label = "NR=5",
             hjust = -0.1, vjust = 1.3, colour = "red", size = 3) +
    annotate("text", x = Inf, y = log10(2 + 1),
             label = "NV=2", hjust = 1.1, vjust = -0.5, colour = "red", size = 3) +
    scale_x_continuous(
      name = expression(log[10](NR + 1)),
      breaks = log10(c(0, 1, 5, 10, 50, 100, 500) + 1),
      labels = c(0, 1, 5, 10, 50, 100, 500)
    ) +
    scale_y_continuous(
      name = expression(log[10](NV + 1)),
      breaks = log10(c(0, 1, 2, 5, 10, 50, 100) + 1),
      labels = c(0, 1, 2, 5, 10, 50, 100)
    ) +
    ggtitle(sprintf("Joint NR vs NV density (entries with NV >= 1) [%s]", tag),
            subtitle = "Blue lines: VAF = 0.10 (dashed) and 0.20 (dotted)") +
    theme_ohchibi_pubr()

  # ---- Panel 6: Per-variant cell-call count under candidate thresholds ----
  combos_f <- list(
    "VAF>=0.10 NR>=5 NV>=2"  = list(vaf_t = 0.10, nr_t = 5,  nv_t = 2),
    "VAF>=0.20 NR>=5 NV>=2"  = list(vaf_t = 0.20, nr_t = 5,  nv_t = 2),
    "VAF>=0.30 NR>=10 NV>=3" = list(vaf_t = 0.30, nr_t = 10, nv_t = 3)
  )

  count_calls_per_var_f <- function(vaf_t, nr_t, nv_t) {
    bin <- (vaf_f >= vaf_t) & (nr_f >= nr_t) & (nv_f >= nv_t)
    rowSums(bin)
  }

  n_binary_list_f <- lapply(combos_f, function(p) count_calls_per_var_f(p$vaf_t, p$nr_t, p$nv_t))

  per_var_f <- rbindlist(lapply(names(n_binary_list_f), function(nm) {
    data.table(combo = nm, n_binary = n_binary_list_f[[nm]])
  }))

  p_per_var_f <- ggplot(per_var_f, aes(x = n_binary, fill = combo)) +
    geom_histogram(binwidth = 1, position = "dodge", colour = "white", linewidth = 0.1) +
    scale_fill_brewer("Threshold combo", palette = "Set1") +
    scale_x_continuous("Cells calling variant (binary = 1)",
                       breaks = seq(0, n_cell_f, by = 4), limits = c(-0.5, n_cell_f + 0.5)) +
    scale_y_continuous("Variants", trans = "log1p",
                       breaks = c(0, 10, 100, 1000, 10000, 1e5)) +
    ggtitle(sprintf("Per-variant cell-call count under candidate thresholds [%s]", tag),
            subtitle = "Look for combos that suppress mass at n=0 / n=1 (likely noise)") +
    theme_ohchibi_pubr()

  # ---- Assemble composition (3 rows x 2 cols) ----
  p_comp_f <- egg::ggarrange(
    p_nr_f,            p_nv_f,
    p_vaf_f,           p_vaf_callable_f,
    p_2d_f,            p_per_var_f,
    nrow = 3, ncol = 2,
    widths  = c(1, 1),
    heights = c(1, 1, 1),
    draw = FALSE
  )

  # Multi-page PDF with individual panels
  cairo_pdf(out_thresh_pdf, width = 9, height = 6, onefile = TRUE)
  print(p_nr_f)
  print(p_nv_f)
  print(p_vaf_f)
  print(p_vaf_callable_f)
  print(p_2d_f)
  print(p_per_var_f)
  dev.off()

  # Composition PNG
  ggsave(out_thresh_png, p_comp_f, width = 14, height = 12, dpi = 180)

  cat(sprintf("  [%s] Wrote %s\n", tag, out_thresh_pdf))
  cat(sprintf("  [%s] Wrote %s\n", tag, out_thresh_png))
}

# ---- Shared-variant exploration ----
build_threshold_exploration(
  nv_f[shared_row_f, , drop = FALSE],
  nr_f[shared_row_f, , drop = FALSE],
  vaf_f[shared_row_f, , drop = FALSE],
  tag           = "shared",
  out_thresh_pdf = out("threshold_exploration_filtered_shared.pdf"),
  out_thresh_png = out("threshold_exploration_filtered_shared.png")
)

# ---- Singleton-variant exploration ----
build_threshold_exploration(
  nv_f[singleton_row_f, , drop = FALSE],
  nr_f[singleton_row_f, , drop = FALSE],
  vaf_f[singleton_row_f, , drop = FALSE],
  tag           = "singleton",
  out_thresh_pdf = out("threshold_exploration_filtered_singleton.pdf"),
  out_thresh_png = out("threshold_exploration_filtered_singleton.png")
)

# ---- 8. Binarize the filtered matrices (split: shared vs singleton) --------
out_binary           <- out("binary_matrix.tsv")
out_binary_shared    <- out("binary_matrix_shared.tsv")
out_binary_singleton <- out("binary_matrix_singleton.tsv")

bin_nv_thr            <- opt[["bin-nv"]]               # anchor / strict NV
bin_nv_shared         <- opt[["bin-nv-shared"]]        # relaxed NV (shared companions)
bin_nr_thr            <- opt[["bin-nr"]]
bin_vaf_relaxed       <- opt[["bin-vaf"]]               # shared, anchored companion cells
bin_vaf_shared_anchor <- opt[["bin-vaf-shared-anchor"]] # shared anchor / strict
bin_vaf_singleton     <- opt[["bin-vaf-singleton"]]     # singleton

cat(sprintf("\n[8/11] Binarizing filtered matrices (NR>=%d throughout)...\n", bin_nr_thr))
cat(sprintf("  singleton rule : NV >= %d & VAF >= %g\n", bin_nv_thr, bin_vaf_singleton))
cat(sprintf("  shared anchor  : NV >= %d & VAF >= %g (>=1 cell unlocks relaxed)\n",
            bin_nv_thr, bin_vaf_shared_anchor))
cat(sprintf("  shared relaxed : NV >= %d & VAF >= %g (companion cells of anchored variants)\n",
            bin_nv_shared, bin_vaf_relaxed))
cat(sprintf("  shared (no anchor) -> every cell requires NV >= %d & VAF >= %g\n",
            bin_nv_thr, bin_vaf_shared_anchor))

# Per-rule gates (NR gate is shared by all rules)
nr_ok           <- nr_f >= bin_nr_thr
depth_strict    <- (nv_f >= bin_nv_thr)    & nr_ok   # singleton + anchor + unanchored shared
depth_relaxed   <- (nv_f >= bin_nv_shared) & nr_ok   # anchored shared companion cells

# Logical call matrix, filled per variant class
called <- matrix(FALSE, nrow = nrow(nv_f), ncol = ncol(nv_f),
                 dimnames = list(NULL, cell_cols))

# --- Singleton variants: strict NV + single strict VAF threshold ---
if (any(singleton_row_f)) {
  call_single <- depth_strict & (vaf_f >= bin_vaf_singleton)
  called[singleton_row_f, ] <- call_single[singleton_row_f, , drop = FALSE]
}

# --- Shared variants: anchor-aware relaxed VAF/NV ---
#   * anchor cell : strict NV AND VAF >= anchor thr;
#   * a variant with >=1 anchor cell calls its other cells at the relaxed
#     NV/VAF thresholds (the anchor cell also passes the relaxed rule);
#   * a variant with no anchor cell must meet the strict anchor rule in every
#     called cell.
if (any(shared_row_f)) {
  anchor_pass  <- depth_strict  & (vaf_f >= bin_vaf_shared_anchor)
  call_relaxed <- depth_relaxed & (vaf_f >= bin_vaf_relaxed)
  has_anchor   <- rowSums(anchor_pass) >= 1

  shared_call <- anchor_pass                                   # default (no anchor): strict
  shared_call[has_anchor, ] <- call_relaxed[has_anchor, , drop = FALSE]
  called[shared_row_f, ] <- shared_call[shared_row_f, , drop = FALSE]

  cat(sprintf("  shared variants with >=1 anchor cell : %d / %d\n",
              sum(has_anchor & shared_row_f), sum(shared_row_f)))
}

binary_mat <- matrix(as.integer(called), nrow = nrow(called), ncol = ncol(called),
                     dimnames = list(NULL, cell_cols))

binary_dt <- data.table(VariantId = nv_sub$VariantId, as.data.table(binary_mat))
fwrite(binary_dt, out_binary, sep = "\t", quote = FALSE)

# Per-class binary matrices (transparency / downstream inspection)
binary_shared_dt    <- binary_dt[shared_row_f]
binary_singleton_dt <- binary_dt[singleton_row_f]
fwrite(binary_shared_dt,    out_binary_shared,    sep = "\t", quote = FALSE)
fwrite(binary_singleton_dt, out_binary_singleton, sep = "\t", quote = FALSE)

n_ones  <- sum(binary_mat == 1)
n_zeros <- sum(binary_mat == 0)
n_vars_with_call <- sum(rowSums(binary_mat) >= 1)
cat(sprintf("  matrix size: %d variants x %d cells\n", nrow(binary_mat), ncol(binary_mat)))
cat(sprintf("  1s (present): %d  |  0s (absent): %d\n", n_ones, n_zeros))
cat(sprintf("  variants with >= 1 call: %d\n", n_vars_with_call))
cat(sprintf("    shared    : %d / %d\n",
            sum(rowSums(binary_mat) >= 1 & shared_row_f),    sum(shared_row_f)))
cat(sprintf("    singleton : %d / %d\n",
            sum(rowSums(binary_mat) >= 1 & singleton_row_f), sum(singleton_row_f)))
cat(sprintf("  Wrote %s\n", out_binary))
cat(sprintf("  Wrote %s\n", out_binary_shared))
cat(sprintf("  Wrote %s\n", out_binary_singleton))

# ---- 9. Subset to variants with Pileup_Verdict == "Pass" (QC) --------------
# Cascade order: QC (Pileup_Verdict) is applied BEFORE the drop-all-zero depth
# step, so QC runs on the full HQ-statistical set (binary_dt / nv_sub / ...).
cat("\n[9/11] Subsetting HQ-statistical variants to those with Pileup_Verdict == 'Pass' (QC)...\n")

out_binary_qc  <- out("binary_matrix_HQRoundStatisticalFilteredPlusQCFiltered.tsv")
out_nv_qc      <- out("NV_HQRoundStatisticalFilteredPlusQCFiltered.tsv")
out_nr_qc      <- out("NR_HQRoundStatisticalFilteredPlusQCFiltered.tsv")
out_vaf_qc     <- out("VAF_HQRoundStatisticalFilteredPlusQCFiltered.tsv")

if (!"Pileup_Verdict" %in% colnames(master_sub)) {
  stop("Column 'Pileup_Verdict' not found in master table.")
}

hqstat_row_ids <- binary_dt$VariantId
pass_ids       <- master_sub[Pileup_Verdict == "Pass", unique(VariantId)]
keep_qc        <- hqstat_row_ids %in% pass_ids

cat(sprintf("  HQ-statistical variants total   : %d\n", length(hqstat_row_ids)))
cat(sprintf("  with Pileup_Verdict == 'Pass'   : %d\n", sum(keep_qc)))
cat(sprintf("  without (removed)               : %d\n", sum(!keep_qc)))

# Subset binary / NV / NR / VAF (all row-aligned with binary_dt)
binary_qc_dt <- binary_dt[keep_qc]
nv_qc_dt     <- nv_sub[keep_qc]
nr_qc_dt     <- nr_sub[keep_qc]
vaf_qc_dt    <- vaf_sub[keep_qc]
fwrite(binary_qc_dt, out_binary_qc, sep = "\t", quote = FALSE)
fwrite(nv_qc_dt,     out_nv_qc,     sep = "\t", quote = FALSE)
fwrite(nr_qc_dt,     out_nr_qc,     sep = "\t", quote = FALSE)
fwrite(vaf_qc_dt,    out_vaf_qc,    sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_binary_qc, nrow(binary_qc_dt), length(cell_cols)))
cat(sprintf("  Wrote %s\n", out_nv_qc))
cat(sprintf("  Wrote %s\n", out_nr_qc))
cat(sprintf("  Wrote %s\n", out_vaf_qc))

# ---- 10. Remove all-zero variants (depth) from the QC-filtered set ---------
# Drop variants whose binarized row has no called cell, now applied to the
# QC-passed set (step 9). Output: ...PlusQCFilteredPlusDepthFiltered.
cat("\n[10/11] Removing all-zero variants from the QC-filtered binary matrix and subsetting NV/NR/VAF...\n")

out_binary_depth <- out("binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered.tsv")
out_nv_depth     <- out("NV_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered.tsv")
out_nr_depth     <- out("NR_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered.tsv")
out_vaf_depth    <- out("VAF_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered.tsv")

binary_qc_mat <- as.matrix(binary_qc_dt[, -1, with = FALSE]); mode(binary_qc_mat) <- "integer"
has_call <- rowSums(binary_qc_mat) >= 1
cat(sprintf("  QC-passed variants total                  : %d\n", nrow(binary_qc_dt)))
cat(sprintf("  variants with 0 across all cells (removed): %d\n", sum(!has_call)))
cat(sprintf("  variants retained (>= 1 call)             : %d\n", sum(has_call)))

# Subset binary / NV / NR / VAF (all row-aligned with binary_qc_dt)
binary_depth_dt <- binary_qc_dt[has_call]
nv_depth_dt     <- nv_qc_dt[has_call]
nr_depth_dt     <- nr_qc_dt[has_call]
vaf_depth_dt    <- vaf_qc_dt[has_call]
fwrite(binary_depth_dt, out_binary_depth, sep = "\t", quote = FALSE)
fwrite(nv_depth_dt,     out_nv_depth,     sep = "\t", quote = FALSE)
fwrite(nr_depth_dt,     out_nr_depth,     sep = "\t", quote = FALSE)
fwrite(vaf_depth_dt,    out_vaf_depth,    sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_binary_depth, nrow(binary_depth_dt), length(cell_cols)))
cat(sprintf("  Wrote %s\n", out_nv_depth))
cat(sprintf("  Wrote %s\n", out_nr_depth))
cat(sprintf("  Wrote %s\n", out_vaf_depth))

# ---- 11. Coverage filter for phylogeny: NR >= bin-nr in >= 70% of cells -----
cat(sprintf("\n[11/11] Retaining variants with NR >= %d in >= 70%% of cells (for phylogeny)...\n",
            bin_nr_thr))

out_binary_phylo <- out("binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv")
out_nv_phylo     <- out("NV_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv")
out_nr_phylo     <- out("NR_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv")
out_vaf_phylo    <- out("VAF_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv")

# Work from the QC+depth NR matrix
nr_depth_mat <- as.matrix(nr_depth_dt[, -1]); mode(nr_depth_mat) <- "numeric"
n_cells_pass <- ncol(nr_depth_mat)
cov_threshold <- 0.70

# For each variant, count cells with NR >= bin_nr_thr
cells_covered <- rowSums(nr_depth_mat >= bin_nr_thr)
keep_phylo <- cells_covered >= (cov_threshold * n_cells_pass)

n_keep_phylo <- sum(keep_phylo)
n_drop_phylo <- sum(!keep_phylo)
cat(sprintf("  QC+depth variants                : %d\n", nrow(nr_depth_dt)))
cat(sprintf("  cells required (70%% of %d)       : %d\n",
            n_cells_pass, ceiling(cov_threshold * n_cells_pass)))
cat(sprintf("  variants passing coverage filter : %d\n", n_keep_phylo))
cat(sprintf("  variants removed                 : %d\n", n_drop_phylo))

# Subset binary
binary_phylo_dt <- binary_depth_dt[keep_phylo]
fwrite(binary_phylo_dt, out_binary_phylo, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_binary_phylo, nrow(binary_phylo_dt), length(cell_cols)))

# Subset NV
nv_phylo_dt <- nv_depth_dt[keep_phylo]
fwrite(nv_phylo_dt, out_nv_phylo, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_nv_phylo, nrow(nv_phylo_dt), length(cell_cols)))

# Subset NR
nr_phylo_dt <- nr_depth_dt[keep_phylo]
fwrite(nr_phylo_dt, out_nr_phylo, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_nr_phylo, nrow(nr_phylo_dt), length(cell_cols)))

# Subset VAF
vaf_phylo_dt <- vaf_depth_dt[keep_phylo]
fwrite(vaf_phylo_dt, out_vaf_phylo, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
            out_vaf_phylo, nrow(vaf_phylo_dt), length(cell_cols)))

# ---- 12. Append mandatory variants (optional) ------------------------------
if (!is.null(opt$mandatory)) {
  cat("\n[12/12] Appending mandatory variants from --mandatory file...\n")

  out_binary_pass_mand <- out("binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatory.tsv")
  out_nv_pass_mand     <- out("NV_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatory.tsv")
  out_nr_pass_mand     <- out("NR_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatory.tsv")
  out_vaf_pass_mand    <- out("VAF_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatory.tsv")

  out_binary_phylo_mand <- out("binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatory.tsv")
  out_nv_phylo_mand     <- out("NV_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatory.tsv")
  out_nr_phylo_mand     <- out("NR_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatory.tsv")
  out_vaf_phylo_mand    <- out("VAF_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatory.tsv")

  mand_dt <- fread(opt$mandatory, sep = "\t", header = TRUE, check.names = FALSE)
  if (ncol(mand_dt) < 1) stop("Mandatory TSV appears to have no columns.")
  setnames(mand_dt, 1, "VariantId")
  mand_ids_all <- unique(mand_dt$VariantId)
  cat(sprintf("  mandatory variants in file       : %d\n", length(mand_ids_all)))

  # Only those that exist in the original NV/NR matrices can be appended
  in_matrix <- mand_ids_all %in% variant_ids
  cat(sprintf("  present in NV/NR matrices        : %d\n", sum(in_matrix)))
  cat(sprintf("  not in NV/NR matrices (skipped)  : %d\n", sum(!in_matrix)))
  mand_ids <- mand_ids_all[in_matrix]

  append_mandatory <- function(binary_existing_dt, nv_existing_dt, nr_existing_dt,
                               vaf_existing_dt, label) {
    existing_ids <- binary_existing_dt$VariantId
    add_ids <- setdiff(mand_ids, existing_ids)
    cat(sprintf("  [%s] existing variants : %d\n", label, length(existing_ids)))
    cat(sprintf("  [%s] already present   : %d\n", label,
                length(intersect(mand_ids, existing_ids))))
    cat(sprintf("  [%s] to append         : %d\n", label, length(add_ids)))

    if (length(add_ids) == 0) {
      return(list(binary = binary_existing_dt,
                  nv     = nv_existing_dt,
                  nr     = nr_existing_dt,
                  vaf    = vaf_existing_dt))
    }

    # Pull rows from the original (unfiltered) NV / NR / VAF matrices
    add_idx <- match(add_ids, variant_ids)
    nv_add  <- nv[add_idx, , drop = FALSE]
    nr_add  <- nr[add_idx, , drop = FALSE]
    vaf_add <- vaf[add_idx, , drop = FALSE]

    # Binarization for mandatory rows: NV > 0 -> 1
    binary_add <- matrix(as.integer(nv_add > 0),
                         nrow = nrow(nv_add), ncol = ncol(nv_add),
                         dimnames = list(NULL, cell_cols))

    nv_add_dt     <- data.table(VariantId = add_ids, as.data.table(nv_add))
    nr_add_dt     <- data.table(VariantId = add_ids, as.data.table(nr_add))
    vaf_add_dt    <- data.table(VariantId = add_ids, as.data.table(round(vaf_add, 6)))
    binary_add_dt <- data.table(VariantId = add_ids, as.data.table(binary_add))

    list(binary = rbind(binary_existing_dt, binary_add_dt, use.names = TRUE),
         nv     = rbind(nv_existing_dt,     nv_add_dt,     use.names = TRUE),
         nr     = rbind(nr_existing_dt,     nr_add_dt,     use.names = TRUE),
         vaf    = rbind(vaf_existing_dt,    vaf_add_dt,    use.names = TRUE))
  }

  # --- Append onto the QC+Depth cascade
  qc <- append_mandatory(binary_depth_dt, nv_depth_dt, nr_depth_dt, vaf_depth_dt,
                         label = "QCFilteredPlusDepthFiltered")
  fwrite(qc$binary, out_binary_pass_mand, sep = "\t", quote = FALSE)
  fwrite(qc$nv,     out_nv_pass_mand,     sep = "\t", quote = FALSE)
  fwrite(qc$nr,     out_nr_pass_mand,     sep = "\t", quote = FALSE)
  fwrite(qc$vaf,    out_vaf_pass_mand,    sep = "\t", quote = FALSE)
  cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
              out_binary_pass_mand, nrow(qc$binary), length(cell_cols)))
  cat(sprintf("  Wrote %s\n", out_nv_pass_mand))
  cat(sprintf("  Wrote %s\n", out_nr_pass_mand))
  cat(sprintf("  Wrote %s\n", out_vaf_pass_mand))

  # --- Append onto the ForPhylogeny cascade
  ph <- append_mandatory(binary_phylo_dt, nv_phylo_dt, nr_phylo_dt, vaf_phylo_dt,
                         label = "ForPhylogeny")
  fwrite(ph$binary, out_binary_phylo_mand, sep = "\t", quote = FALSE)
  fwrite(ph$nv,     out_nv_phylo_mand,     sep = "\t", quote = FALSE)
  fwrite(ph$nr,     out_nr_phylo_mand,     sep = "\t", quote = FALSE)
  fwrite(ph$vaf,    out_vaf_phylo_mand,    sep = "\t", quote = FALSE)
  cat(sprintf("  Wrote %s (%d variants x %d cells)\n",
              out_binary_phylo_mand, nrow(ph$binary), length(cell_cols)))
  cat(sprintf("  Wrote %s\n", out_nv_phylo_mand))
  cat(sprintf("  Wrote %s\n", out_nr_phylo_mand))
  cat(sprintf("  Wrote %s\n", out_vaf_phylo_mand))
} else {
  cat("\n[12/12] --mandatory not provided; skipping mandatory-variant append step.\n")
}

# ---- 13. Per-cell SNV / indel counts across filtering stages ---------------
cat("\n[13/14] Building per-cell SNV and indel count plots across filtering stages...\n")

out_counts_tsv <- out("per_cell_variant_counts_by_stage.tsv")
out_snv_pdf    <- out("per_cell_SNV_counts_by_stage.pdf")
out_snv_png    <- out("per_cell_SNV_counts_by_stage.png")
out_indel_pdf  <- out("per_cell_INDEL_counts_by_stage.pdf")
out_indel_png  <- out("per_cell_INDEL_counts_by_stage.png")
out_combo_pdf  <- out("per_cell_SNV_INDEL_counts_by_stage.pdf")
out_combo_png  <- out("per_cell_SNV_INDEL_counts_by_stage.png")

# Helper: parse VariantId (chr_pos_REF_ALT) → SNV (REF and ALT both length 1).
classify_snv_indel <- function(variant_ids) {
  parts <- tstrsplit(variant_ids, "_", fixed = TRUE)
  if (length(parts) < 4) {
    stop("VariantId does not look like 'chr_pos_REF_ALT'; cannot classify SNV/indel.")
  }
  ref <- parts[[3]]
  alt <- parts[[4]]
  nchar(ref) == 1L & nchar(alt) == 1L
}

# Helper: per-cell counts (SNV and INDEL) for a binary data.table that has
# VariantId as its first column and cells in the remaining columns.
per_cell_counts <- function(binary_dt, stage_label) {
  if (is.null(binary_dt) || nrow(binary_dt) == 0) {
    return(data.table(stage = character(), cell = character(),
                      VariantType = character(), Count = integer()))
  }
  vids <- binary_dt$VariantId
  is_snv <- classify_snv_indel(vids)
  mat <- as.matrix(binary_dt[, -1, with = FALSE]); mode(mat) <- "integer"
  snv_counts   <- colSums(mat[is_snv,  , drop = FALSE], na.rm = TRUE)
  indel_counts <- colSums(mat[!is_snv, , drop = FALSE], na.rm = TRUE)
  rbind(
    data.table(stage = stage_label, cell = colnames(mat),
               VariantType = "SNV",   Count = as.integer(snv_counts)),
    data.table(stage = stage_label, cell = colnames(mat),
               VariantType = "INDEL", Count = as.integer(indel_counts))
  )
}

# Assemble the stages (PlusMandatory ones only if --mandatory was given).
stages <- list(
  list(label = "HQRoundStatisticalFiltered",
       dt    = binary_hqstat_dt),
  list(label = "HQRoundStatisticalFilteredPlusQCFiltered",
       dt    = binary_qc_dt),
  list(label = "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered",
       dt    = binary_depth_dt),
  list(label = "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny",
       dt    = binary_phylo_dt)
)
if (!is.null(opt$mandatory) && exists("qc")) {
  stages[[length(stages) + 1]] <- list(
    label = "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatory",
    dt    = qc$binary)
  stages[[length(stages) + 1]] <- list(
    label = "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatory",
    dt    = ph$binary)
}

counts_dt <- rbindlist(lapply(stages, function(s) per_cell_counts(s$dt, s$label)))
stage_levels <- vapply(stages, function(s) s$label, character(1))
counts_dt[, stage       := factor(stage, levels = stage_levels)]
counts_dt[, VariantType := factor(VariantType, levels = c("SNV", "INDEL"))]

# Drop zeros so log10 is well-defined (record how many we drop, per stage).
n_zero <- counts_dt[Count == 0, .N, by = .(stage, VariantType)]
if (nrow(n_zero) > 0) {
  cat("  cells with 0 counts (excluded from log10 axis):\n")
  print(n_zero)
}

fwrite(counts_dt, out_counts_tsv, sep = "\t", quote = FALSE)
cat(sprintf("  Wrote %s\n", out_counts_tsv))

# Compact axis labels on the x-axis (one per stage).
short_label_map <- c(
  "HQRoundStatisticalFiltered"                                                          = "HQStat",
  "HQRoundStatisticalFilteredPlusQCFiltered"                                            = "HQStat\n+QC",
  "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered"                           = "HQStat\n+QC\n+Depth",
  "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny"               = "HQStat\n+QC\n+Depth\n+Phylo",
  "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatory"              = "HQStat\n+QC\n+Depth\n+Mand",
  "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatory"  = "HQStat\n+QC\n+Depth\n+Phylo\n+Mand"
)
short_labels <- short_label_map[stage_levels]

stage_palette <- c(
  "HQRoundStatisticalFiltered"                                                          = "#9D9D9D",
  "HQRoundStatisticalFilteredPlusQCFiltered"                                            = "#5DADE2",
  "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered"                           = "#1F77B4",
  "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny"               = "#1A5276",
  "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatory"              = "#E59866",
  "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatory"  = "#BA4A00"
)

make_stage_plot <- function(variant_type, ylab, palette) {
  d <- counts_dt[VariantType == variant_type & Count > 0]
  if (nrow(d) == 0) {
    return(ggplot() + theme_void() +
             ggtitle(sprintf("No %s data to plot", variant_type)))
  }
  med <- d[, .(med = median(Count)), by = stage]
  med[, label := format(med, big.mark = ",", scientific = FALSE)]

  ggplot(d, aes(x = stage, y = Count, fill = stage, colour = stage)) +
    geom_jitter(width = 0.18, height = 0, alpha = 0.55, size = 1.0,
                shape = 16) +
    geom_boxplot(width = 0.45, alpha = 0.7, outlier.shape = NA,
                 colour = "black", linewidth = 0.4) +
    geom_text(data = med, aes(x = stage, y = med, label = label),
              inherit.aes = FALSE, vjust = -0.6, size = 2.8) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_colour_manual(values = palette, guide = "none") +
    scale_x_discrete(labels = short_labels) +
    scale_y_log10(labels = scales::label_comma(accuracy = 1)) +
    annotation_logticks(sides = "l", short = unit(0.05, "cm"),
                        mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    labs(x = NULL, y = ylab,
         subtitle = sprintf("n = %d cells per stage", uniqueN(d$cell))) +
    theme_ohchibi_pubr() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                     size = 8))
}

p_snv <- make_stage_plot(
  "SNV",
  expression("Number of SNVs per cell (log"[10]*" scale)"),
  stage_palette
) + ggtitle("Per-cell SNV burden across filtering stages")

p_indel <- make_stage_plot(
  "INDEL",
  expression("Number of indels per cell (log"[10]*" scale)"),
  stage_palette
) + ggtitle("Per-cell indel burden across filtering stages")

# Width scales with number of stages (each stage gets ~1.1 inches)
plot_width <- max(6, 1.6 * length(stage_levels) + 1)

cairo_pdf(out_snv_pdf, width = plot_width, height = 5)
print(p_snv)
dev.off()
ggsave(out_snv_png, p_snv, width = plot_width, height = 5, dpi = 200, bg = "white")
cat(sprintf("  Wrote %s\n", out_snv_pdf))
cat(sprintf("  Wrote %s\n", out_snv_png))

cairo_pdf(out_indel_pdf, width = plot_width, height = 5)
print(p_indel)
dev.off()
ggsave(out_indel_png, p_indel, width = plot_width, height = 5, dpi = 200, bg = "white")
cat(sprintf("  Wrote %s\n", out_indel_pdf))
cat(sprintf("  Wrote %s\n", out_indel_png))

p_combo <- egg::ggarrange(p_snv, p_indel, nrow = 1, ncol = 2, draw = FALSE)
cairo_pdf(out_combo_pdf, width = plot_width * 2, height = 5)
print(p_combo)
dev.off()
ggsave(out_combo_png, p_combo, width = plot_width * 2, height = 5, dpi = 200, bg = "white")
cat(sprintf("  Wrote %s\n", out_combo_pdf))
cat(sprintf("  Wrote %s\n", out_combo_png))

# ---- 14. Drop all-zero variants from the PlusMandatory matrices -------------
# The mandatory-append step (step 12) force-includes priority variants and
# binarizes them as NV>0 -> 1, so a mandatory variant with no variant reads in
# any cell enters as an all-zero row. Here we remove such empty variants from
# both PlusMandatory cascades, using the binary matrix as the signal and
# applying the same row mask to the binary / NV / NR / VAF matrices.
if (!is.null(opt$mandatory) && exists("qc") && exists("ph")) {
  cat("\n[14/14] Removing all-zero variants (by binary signal) from PlusMandatory matrices...\n")

  drop_all_zero <- function(mats, label, out_binary, out_nv, out_nr, out_vaf) {
    bmat <- as.matrix(mats$binary[, -1, with = FALSE]); mode(bmat) <- "integer"
    keep <- rowSums(bmat) >= 1
    cat(sprintf("  [%s] variants total           : %d\n", label, nrow(mats$binary)))
    cat(sprintf("  [%s] all-zero variants removed : %d\n", label, sum(!keep)))
    cat(sprintf("  [%s] variants retained         : %d\n", label, sum(keep)))
    fwrite(mats$binary[keep], out_binary, sep = "\t", quote = FALSE)
    fwrite(mats$nv[keep],     out_nv,     sep = "\t", quote = FALSE)
    fwrite(mats$nr[keep],     out_nr,     sep = "\t", quote = FALSE)
    fwrite(mats$vaf[keep],    out_vaf,    sep = "\t", quote = FALSE)
    cat(sprintf("  [%s] Wrote %s (%d variants x %d cells)\n",
                label, out_binary, sum(keep), length(cell_cols)))
    cat(sprintf("  [%s] Wrote %s\n", label, out_nv))
    cat(sprintf("  [%s] Wrote %s\n", label, out_nr))
    cat(sprintf("  [%s] Wrote %s\n", label, out_vaf))
  }

  drop_all_zero(
    qc, "QCFiltered+Mandatory",
    out("binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatoryNonEmpty.tsv"),
    out("NV_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatoryNonEmpty.tsv"),
    out("NR_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatoryNonEmpty.tsv"),
    out("VAF_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatoryNonEmpty.tsv")
  )

  drop_all_zero(
    ph, "ForPhylogeny+Mandatory",
    out("binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatoryNonEmpty.tsv"),
    out("NV_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatoryNonEmpty.tsv"),
    out("NR_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatoryNonEmpty.tsv"),
    out("VAF_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogenyPlusMandatoryNonEmpty.tsv")
  )
} else {
  cat("\n[14/14] --mandatory not provided; no PlusMandatory matrices to clean.\n")
}

cat("\nDone.\n")
