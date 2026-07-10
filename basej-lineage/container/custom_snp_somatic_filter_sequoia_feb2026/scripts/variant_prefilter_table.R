#!/usr/bin/env Rscript
## Build the variant PREFILTER table: per-variant Pileup_* aggregation + SEQUOIA
## second-pass columns. This was sections 4-5 of master_table.R, extracted so it can
## run BEFORE VAF_SPLIT_VARIANTS_HEXBIN (which needs SecondRunSequoia_Rho/_Germline_qval
## + Pileup_Verdict for --master) while CUSTOM_VARIANT_FILTER_PROVENANCE runs AFTER
## VAF_SPLIT and consumes this table (so master_table.R no longer recomputes 4-5).
##
## MERGE step: the per-variant Pileup_* aggregation is now produced PER-CHROMOSOME by
## variant_prefilter_table_perchr.R (fanned out, so the full ~140M-row pileup pool is
## never loaded in one task). This step rbinds those small partials and joins the
## group-level SEQUOIA second-pass output. Output is identical to the former single pass.
##
## Usage: Rscript variant_prefilter_table.R <group_id>
## Reads from the current (Nextflow work) dir:
##   prefilter_pileup_agg_<group>_*.tsv   (per-chr aggregated Pileup_* partials)
##   *_filtering_all.txt                  (SEQUOIA second-pass merged output)
## Writes:
##   variant_prefilter_table_<group>.tsv   (VariantId + Pileup_* + SecondRunSequoia_*)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript variant_prefilter_table.R <group_id>")
group_id <- args[1]

suppressPackageStartupMessages(library(data.table))

agg_files <- list.files(".", pattern = paste0("^prefilter_pileup_agg_", group_id, "_.*\\.tsv$"))
seq_file  <- list.files(".", pattern = "_filtering_all\\.txt$")[1]
if (length(agg_files) == 0) stop("No prefilter_pileup_agg_", group_id, "_*.tsv partials found")
if (is.na(seq_file))        stop("No *_filtering_all.txt (SEQUOIA) file found")

cat("Files found:\n")
cat("  Per-chr pileup aggregates : ", length(agg_files), "files\n")
cat("  SEQUOIA: ", seq_file, "\n\n")

# ── Pileup: rbind the per-chr aggregated partials (each variant lives on one chr, so the
# partials are disjoint; the final sort by VariantId makes row order identical to the
# former single-pass aggregation over the pooled pile). ──────────────────────────────
pileup_agg <- rbindlist(lapply(agg_files, fread, sep = "\t", header = TRUE,
                               data.table = TRUE, showProgress = FALSE), use.names = TRUE)
cat("Per-variant pileup rows (all chr):", nrow(pileup_agg), "\n")

# ── SEQUOIA second-pass filter ──────────────────────────────────────────────────
seq_df <- read.table(seq_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(seq_df) > 0 && length(grep("_.*_.*_", seq_df[[1]][1])) > 0) {
  seq_df$VariantId <- seq_df[[1]]
} else {
  seq_df$VariantId <- rownames(seq_df)
}
if (!"SecondPassFilter" %in% colnames(seq_df)) {
  seq_df$SecondPassFilter <- ifelse(
    !is.na(seq_df$Depth_filter)  & seq_df$Depth_filter == TRUE &
    !is.na(seq_df$Germline)      & seq_df$Germline == 1L &
    !is.na(seq_df$Beta_binomial) & seq_df$Beta_binomial == 1L,
    "Pass", "Fail"
  )
}
seq_sel <- seq_df[, c("VariantId", setdiff(colnames(seq_df), "VariantId")), drop = FALSE]
colnames(seq_sel)[-1] <- paste0("SecondRunSequoia_", colnames(seq_sel)[-1])

# ── Join pileup_agg ⋈ sequoia (full outer on VariantId) and write ──────────────
prefilter <- merge(pileup_agg, as.data.table(seq_sel), by = "VariantId", all = TRUE)
prefilter <- as.data.frame(prefilter, stringsAsFactors = FALSE)
prefilter <- prefilter[order(prefilter$VariantId), ]

write.table(prefilter, paste0("variant_prefilter_table_", group_id, ".tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Prefilter table written:", nrow(prefilter), "variants x", ncol(prefilter), "columns\n")
