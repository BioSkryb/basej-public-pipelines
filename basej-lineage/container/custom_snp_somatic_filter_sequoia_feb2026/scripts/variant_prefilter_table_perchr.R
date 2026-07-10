#!/usr/bin/env Rscript
## variant_prefilter_table_perchr.R <group_id> <chr>
## Per-chromosome HALF of variant_prefilter_table.R, fanned out so the full ~140M-row
## pileup pool is never rbind'd in a single task. Reads ONE chromosome's
## res_pileup_all_group_<group>_<chr>.tsv, performs the per-variant Pileup_* aggregation
## (every VariantId lives on exactly one chromosome, so the aggregation is fully
## chromosome-independent), and writes the partial table
## prefilter_pileup_agg_<group>_<chr>.tsv. The merge step (variant_prefilter_table.R)
## rbinds these partials (cheap) and joins the group-level SEQUOIA second-pass output.
## Producing identical results to the former single monolithic pass.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript variant_prefilter_table_perchr.R <group_id> <chr>")
group_id <- args[1]; chr <- args[2]

suppressPackageStartupMessages(library(data.table))

pileup_file <- sprintf("res_pileup_all_group_%s_%s.tsv", group_id, chr)
if (!file.exists(pileup_file)) stop("Pileup file not found: ", pileup_file)

num_cols   <- c("NUM_FRAGMENTS_ALLQ_MQ_BQ_F", "NUM_FRAGMENTS_ALLQ_MQ_BQ_R",
                "NUM_FRAGMENTS_HQ_MQ_BQ_F", "NUM_FRAGMENTS_HQ_MQ_BQ_R",
                "MEDIAN_AS_VARIANT_READS", "PROP_BASES_CLIPPED",
                "PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F", "PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R",
                "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F", "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
                "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F", "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
                "NUM_FRAGMENTS_ALLQ_POSITION", "NUM_FRAGMENTS_HQ_POSITION",
                "PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F", "PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R")
filt_cols  <- c("AS_Filter", "PropClipped_Filter", "BPPos_Filter", "Depth_Filter", "Verdict")
const_cols <- c("CHROM", "POS", "REF", "ALT")
yn_cols    <- c("PresentVCF_NotInBAM", "PresentBAM_NotInVCF")

# Column projection: parse only the fields the aggregation uses (drops SampleId + any
# unused columns), so the per-chr pileup parse is cheaper. Aggregation result is identical
# (Pileup_NumSamples = .N counts rows, independent of which columns are read).
sel <- c("VariantId", const_cols, num_cols, filt_cols, yn_cols)
dt <- fread(pileup_file, sep = "\t", header = TRUE, data.table = TRUE,
            showProgress = FALSE, select = sel)
cat("Pileup rows loaded (", chr, "):", nrow(dt), "\n")

# Grouped passes (all keyby=VariantId → identical group order, so column-bind is safe).
agg_n     <- dt[, .(Pileup_NumSamples = .N), keyby = VariantId]
agg_const <- dt[, lapply(.SD, function(x) x[1L]),                          keyby = VariantId, .SDcols = const_cols]
agg_mean  <- dt[, lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 4)), keyby = VariantId, .SDcols = num_cols]
agg_pass  <- dt[, lapply(.SD, function(x) as.integer(sum(x == "Pass"))),   keyby = VariantId, .SDcols = filt_cols]
agg_yn    <- dt[, lapply(.SD, function(x) x[1L]),                          keyby = VariantId, .SDcols = yn_cols]
setnames(agg_const, const_cols, paste0("Pileup_", const_cols))
setnames(agg_mean,  num_cols,   paste0("Pileup_Mean_", num_cols))
setnames(agg_pass,  filt_cols,  paste0("Pileup_", filt_cols, "_PassCount"))
setnames(agg_yn,    yn_cols,    paste0("Pileup_", yn_cols))

pileup_agg <- cbind(agg_n,
                    agg_const[, -1, with = FALSE],
                    agg_mean[,  -1, with = FALSE],
                    agg_pass[,  -1, with = FALSE],
                    agg_yn[,    -1, with = FALSE])
pileup_agg[, Pileup_Verdict := ifelse(Pileup_Verdict_PassCount > 0, "Pass", "Fail")]

out <- sprintf("prefilter_pileup_agg_%s_%s.tsv", group_id, chr)
fwrite(pileup_agg, out, sep = "\t", na = "NA")
cat("Per-chr pileup aggregation written:", out, "-", nrow(pileup_agg), "variants\n")
