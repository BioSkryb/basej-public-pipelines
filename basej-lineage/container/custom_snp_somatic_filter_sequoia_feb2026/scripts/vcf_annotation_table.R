#!/usr/bin/env Rscript
## Build a per-variant annotation table for downstream VCF annotation.
## Called from CUSTOM_VARIANT_FILTER_PROVENANCE (or a dedicated annotation process).
## Usage: Rscript vcf_annotation_table.R <group_id>
## Input:  variant_master_filter_table_<group_id>.tsv  (current working directory)
## Output: vcf_annotation_table_<group_id>.tsv         (current working directory)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript vcf_annotation_table.R <group_id>")
group_id <- args[1]

master_file <- paste0("variant_master_filter_table_", group_id, ".tsv")
if (!file.exists(master_file))
    stop("Master table not found: ", master_file)

cat("Reading master table:", master_file, "\n")
master <- read.table(master_file, sep = "\t", header = TRUE, quote = "",
                     comment.char = "", stringsAsFactors = FALSE, na.strings = "NA")
cat("  Total variants:", nrow(master), "\n\n")

# ── 1. Subset: variants that passed all three sequential filters OR are priority ─
# Pipeline order: VEP → Bulk → Statistical (Binom/Betabinom)
passed_all <- !is.na(master$VEP_FILTER_STATUS)          & master$VEP_FILTER_STATUS          == "PASS" &
              !is.na(master$Bulk_RemainingAfterBulk)     & master$Bulk_RemainingAfterBulk     == "Pass" &
              !is.na(master$Binom_BinomialBetabinomialFilter) & master$Binom_BinomialBetabinomialFilter == "Pass"
is_priority <- !is.na(master$PriorityVariant) & master$PriorityVariant == "Yes"

vep_pass <- master[passed_all | is_priority, ]
cat("Passed all filters (VEP+Bulk+Binom):  ", sum(passed_all), "variants\n")
cat("Priority variants (rescued):           ", sum(is_priority & !passed_all), "additional\n")
cat("Total in annotation table:             ", nrow(vep_pass), "variants\n\n")

# ── 2. Select and rename columns for downstream VCF annotation ────────────────
# Sequoia columns were renamed from Sequoia_ to SecondRunSequoia_ in master_table.R;
# support both names so the script works on old and new runs alike.
seq_qval_src <- if ("SecondRunSequoia_Germline_qval" %in% colnames(vep_pass))
                    "SecondRunSequoia_Germline_qval" else "Sequoia_Germline_qval"
seq_rho_src  <- if ("SecondRunSequoia_Rho" %in% colnames(vep_pass))
                    "SecondRunSequoia_Rho" else "Sequoia_Rho"

keep_cols <- c(
    "VariantId",
    "PriorityVariant",
    "Bulk_RemainingAfterBulk",
    "Binom_Mean_Depth", "Binom_Depth_filter",
    "Binom_Germline_pval", "Binom_Germline_pval_log10",
    "Binom_Rho", "Binom_Gender",
    "Binom_Germline_qval", "Binom_Germline_qval_log10",
    "Binom_Germline_filter", "Binom_Betabinomial_filter",
    "Binom_Verdict", "Binom_BinomialBetabinomialFilter",
    "VEP_Existing_variation", "VEP_AF", "VEP_MAX_AF",
    "VEP_FILTER_STATUS", "VEP_FILTER_REASON",
    "Pileup_AS_Filter_PassCount", "Pileup_PropClipped_Filter_PassCount",
    "Pileup_BPPos_Filter_PassCount",
    "Pileup_Depth_Filter_PassCount", "Pileup_Verdict_PassCount",
    "Pileup_PresentVCF_NotInBAM", "Pileup_PresentBAM_NotInVCF",
    "Pileup_Verdict",
    seq_qval_src, seq_rho_src
)

# Verify all columns are present before subsetting
missing <- setdiff(keep_cols, colnames(vep_pass))
if (length(missing) > 0)
    stop("Missing columns in master table: ", paste(missing, collapse = ", "))

annot <- vep_pass[, keep_cols]

# Rename Sequoia source columns to the canonical output names
colnames(annot)[colnames(annot) == seq_qval_src] <- "SEQUOIA_SecondPass_Germline_qval"
colnames(annot)[colnames(annot) == seq_rho_src]  <- "SEQUOIA_SecondPass_Rho"

# ── Add log10 of SEQUOIA_SecondPass_Germline_qval ─────────────────────────────
# log10(0) = -Inf is kept as-is; NA propagates naturally
annot$SEQUOIA_SecondPass_Germline_qval_log10 <- log10(annot$SEQUOIA_SecondPass_Germline_qval)

cat("Output columns (", ncol(annot), "):\n")
cat(paste0("  ", colnames(annot), "\n"), sep = "")
cat("\nVariants in annotation table:", nrow(annot), "\n\n")

# ── 3. Write output ────────────────────────────────────────────────────────────
out_file <- paste0("vcf_annotation_table_", group_id, ".tsv")
write.table(annot, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Written:", out_file, "\n")
