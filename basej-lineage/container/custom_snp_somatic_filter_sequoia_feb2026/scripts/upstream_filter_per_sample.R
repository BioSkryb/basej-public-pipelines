#!/usr/bin/env Rscript
## Compute per-sample variant counts at each upstream filter stage.
##
## All four stages are counted from the same source (the pre-bulk VCF-derived
## df_nv matrix from MERGE_PROCESSED_VCF), so counts are directly comparable
## and form a strict funnel.  Filter membership at each stage is determined by
## the corresponding column in variant_master_filter_table.
##
## Reads:
##   1. variant_master_filter_table_<group>.tsv  — per-variant filter status
##   2. df_nv_prebulk  — full pre-bulk NV matrix from MERGE_PROCESSED_VCF
##                        (rows = VariantId CHROM_POS_REF_ALT, cols = SampleId,
##                         values = AD in "REF,ALT" format)
##   3. mat_nv_files   — per-chromosome pileup NV matrices (no longer used for
##                        stage counting; retained as positional args for
##                        backwards compatibility with the Nextflow module)
##
## Writes:
##   upstream_filter_per_sample_<group>.tsv
##   Columns: sample | scheme | NumberOfSNVs | NumberOfIndels
##   Stages (pipeline order):
##     pre_bulk   — all variants with NV > 0 in the pre-bulk matrix.
##     post_vep   — variants where VEP_FILTER_STATUS == "PASS".
##                  When the column is absent (VEP disabled) equals pre_bulk.
##     post_bulk  — variants where Bulk_RemainingAfterBulk == "Pass".
##                  When the column is absent (no bulk VCF) equals post_vep.
##     post_binom — variants where Binom_BinomialBetabinomialFilter == "Pass".
##
## Called from CUSTOM_VARIANT_FILTER_PROVENANCE after master_table.R.
## Usage: Rscript upstream_filter_per_sample.R <group_id> <master_table_tsv> <df_nv_prebulk> [mat_nv_files...]

library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) stop(
    "Usage: upstream_filter_per_sample.R <group_id> <master_table_tsv> <df_nv_prebulk> [mat_nv_files...]"
)

group_id      <- args[1]
master_tsv    <- args[2]
df_nv_prebulk <- args[3]

cat("[upstream_filter_per_sample] Group:", group_id, "\n")
cat("[upstream_filter_per_sample] Pre-bulk NV file:", df_nv_prebulk, "\n")

# ── 0. Load pre-bulk NV matrix ───────────────────────────────────────────────
prebulk_mat <- read.delim(df_nv_prebulk, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
cat("[upstream_filter_per_sample] Pre-bulk variants:", nrow(prebulk_mat), "\n")
cat("[upstream_filter_per_sample] Pre-bulk samples:",  ncol(prebulk_mat), "\n")

# df_nv stores AD values as "REF,ALT"; extract ALT count (NV).
prebulk_nv <- matrix(
    suppressWarnings(as.integer(sub(".*,", "", as.matrix(prebulk_mat)))),
    nrow = nrow(prebulk_mat), ncol = ncol(prebulk_mat),
    dimnames = dimnames(prebulk_mat)
)
prebulk_nv[is.na(prebulk_nv)] <- 0L

prebulk_long <- reshape2::melt(prebulk_nv, varnames=c("VariantId","SampleId"), value.name="NV")
prebulk_long$VariantId <- as.character(prebulk_long$VariantId)
prebulk_long$SampleId  <- as.character(prebulk_long$SampleId)
prebulk_long           <- prebulk_long[prebulk_long$NV > 0, ]
cat("[upstream_filter_per_sample] Pre-bulk variant-sample pairs (NV>0):", nrow(prebulk_long), "\n")

prebulk_samples <- unique(prebulk_long$SampleId)

# Classify SNV vs INDEL from VariantId (CHROM_POS_REF_ALT).
pb_vtype <- data.frame(VariantId = unique(prebulk_long$VariantId), stringsAsFactors=FALSE)
pb_vtype$vtype <- sapply(pb_vtype$VariantId, function(vid) {
    parts <- strsplit(vid, "_", fixed=TRUE)[[1]]
    n     <- length(parts)
    if (n < 4) return("SNV")
    if (nchar(parts[n]) == nchar(parts[n-1])) "SNV" else "INDEL"
})
prebulk_df         <- merge(prebulk_long, pb_vtype, by="VariantId")
prebulk_df$is_snv  <- prebulk_df$vtype == "SNV"
prebulk_df$is_indel <- prebulk_df$vtype == "INDEL"

# ── 1. Load master table filter columns ─────────────────────────────────────
mt <- read.delim(master_tsv, stringsAsFactors=FALSE, check.names=FALSE)
required_cols <- c("VariantId", "Binom_BinomialBetabinomialFilter")
missing_req   <- setdiff(required_cols, colnames(mt))
if (length(missing_req) > 0) {
    stop("Missing columns in master_table: ", paste(missing_req, collapse=", "))
}
has_vep  <- "VEP_FILTER_STATUS"       %in% colnames(mt)
has_bulk <- "Bulk_RemainingAfterBulk" %in% colnames(mt)

join_cols <- c("VariantId", "Binom_BinomialBetabinomialFilter",
               if (has_vep)  "VEP_FILTER_STATUS"       else NULL,
               if (has_bulk) "Bulk_RemainingAfterBulk" else NULL)
pb_joined <- merge(prebulk_df, mt[, join_cols], by="VariantId", all.x=TRUE)

# ── 2. Helper: count SNV/INDEL per sample from a prebulk_df subset ───────────
count_prebulk <- function(df_sub, stage_name) {
    if (nrow(df_sub) == 0) {
        return(data.frame(sample=prebulk_samples, scheme=stage_name,
                          NumberOfSNVs=0L, NumberOfIndels=0L,
                          stringsAsFactors=FALSE))
    }
    agg <- aggregate(cbind(NumberOfSNVs=is_snv, NumberOfIndels=is_indel) ~ SampleId,
                     data=df_sub, FUN=sum)
    missing_smp <- setdiff(prebulk_samples, agg$SampleId)
    if (length(missing_smp) > 0) {
        agg <- rbind(agg, data.frame(SampleId=missing_smp, NumberOfSNVs=0L, NumberOfIndels=0L))
    }
    data.frame(sample=agg$SampleId, scheme=stage_name,
               NumberOfSNVs=as.integer(agg$NumberOfSNVs),
               NumberOfIndels=as.integer(agg$NumberOfIndels),
               stringsAsFactors=FALSE)
}

# ── 3. Compute each stage ────────────────────────────────────────────────────
# post_vep: only strict PASS (NA = not processed by VEP → treated as REMOVED).
pb_post_vep <- if (has_vep) {
    pb_joined[!is.na(pb_joined$VEP_FILTER_STATUS) & pb_joined$VEP_FILTER_STATUS == "PASS", ]
} else {
    pb_joined
}

# post_bulk: variants surviving bulk filter; equals post_vep when no bulk VCF.
pb_post_bulk <- if (has_bulk) {
    pb_joined[!is.na(pb_joined$Bulk_RemainingAfterBulk) & pb_joined$Bulk_RemainingAfterBulk == "Pass", ]
} else {
    pb_post_vep
}

# post_binom: variants passing 1st-pass binomial/beta-binomial filter.
pb_post_binom <- pb_joined[!is.na(pb_joined$Binom_BinomialBetabinomialFilter) &
                             pb_joined$Binom_BinomialBetabinomialFilter == "Pass", ]

# ── 4. Assemble output in pipeline order ────────────────────────────────────
rows <- rbind(
    count_prebulk(pb_joined,    "pre_bulk"),
    count_prebulk(pb_post_vep,  "post_vep"),
    count_prebulk(pb_post_bulk, "post_bulk"),
    count_prebulk(pb_post_binom,"post_binom")
)

out_file <- paste0("upstream_filter_per_sample_", group_id, ".tsv")
write.table(rows, out_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("[upstream_filter_per_sample] Written:", out_file, "(", nrow(rows), "rows)\n")
