#!/usr/bin/env Rscript
# mandatory_variants_qc_status.R
# ---------------------------------------------------------------------------
# For every user-supplied mandatory variant, report whether it WOULD have passed
# the HQ-statistical filter and the pileup-QC filter on its own merits (mandatory
# variants are force-kept / rescued downstream, so this shows which are "real"
# passes vs rescued).
#   HQStat : Scheme_HQStat == 1        (from scheme_membership)
#   QC     : Pileup_Verdict == "Pass"  (from the variant master filter table)
# A mandatory variant absent from the NR/NV matrices (e.g. no read evidence) gets
# HQStat = NA; one absent from the master gets QC = NA.
#
# Usage: Rscript mandatory_variants_qc_status.R <group> <master_table.tsv> <scheme_membership.tsv> <mandatory_variants.txt>

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: mandatory_variants_qc_status.R <group> <master> <scheme_membership> <mandatory_variants>")
group <- args[1]; master_f <- args[2]; memb_f <- args[3]; mand_f <- args[4]

out_tsv <- sprintf("mandatory_variants_qc_status_%s.tsv", group)
hdr <- c("VariantId", "InMatrices", "HQStat", "QC_Pileup_Verdict", "PassesBoth")

# ── mandatory variant list (one VariantId per line; tolerate header / extra cols) ──
mand <- character(0)
if (file.exists(mand_f) && file.info(mand_f)$size > 0) {
  ln <- readLines(mand_f, warn = FALSE)
  ln <- trimws(ln)
  ln <- ln[nzchar(ln)]
  ln <- sub("[ \t].*$", "", ln)          # first whitespace/tab-delimited token
  ln <- ln[ln != "VariantId"]            # drop a header if present
  mand <- unique(ln)
}

if (length(mand) == 0) {
  write.table(data.frame(matrix(ncol = length(hdr), nrow = 0, dimnames = list(NULL, hdr))),
              out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[mandatory_variants_qc_status] %s: no mandatory variants provided; wrote empty table.\n", group))
  quit(save = "no", status = 0)
}

# ── HQStat membership ──
memb <- fread(memb_f, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("VariantId", "Scheme_HQStat") %in% colnames(memb)))
  stop("scheme_membership missing VariantId / Scheme_HQStat")
hq_lookup <- setNames(memb$Scheme_HQStat, memb$VariantId)

# ── QC (Pileup_Verdict) ──
master <- fread(master_f, sep = "\t", header = TRUE, data.table = FALSE,
                select = c("VariantId", "Pileup_Verdict"))
pv_lookup <- setNames(as.character(master$Pileup_Verdict), master$VariantId)

# ── per mandatory variant ──
in_mat <- mand %in% names(hq_lookup)
in_mst <- mand %in% names(pv_lookup)
# Treat both "absent from the index" and "present but value is NA" (e.g. Pileup_Verdict
# read as NA) as the "NA" category, so no real NA leaks into the vectors / summary counts.
hqstat <- ifelse(!in_mat | is.na(hq_lookup[mand]), "NA", ifelse(hq_lookup[mand] >= 1, "Pass", "Fail"))
qc     <- ifelse(!in_mst | is.na(pv_lookup[mand]), "NA", ifelse(pv_lookup[mand] == "Pass", "Pass", "Fail"))
passes_both <- ifelse(hqstat == "Pass" & qc == "Pass", "Yes", "No")

res <- data.frame(VariantId = mand,
                  InMatrices = ifelse(in_mat, "Yes", "No"),
                  HQStat = hqstat,
                  QC_Pileup_Verdict = qc,
                  PassesBoth = passes_both,
                  stringsAsFactors = FALSE)
write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# ── summary to stdout ──
n <- nrow(res)
cat(sprintf("[mandatory_variants_qc_status] %s: %d mandatory variants\n", group, n))
cat(sprintf("  passed both (HQStat & QC) : %d\n", sum(res$PassesBoth == "Yes")))
cat(sprintf("  HQStat Pass / Fail / NA   : %d / %d / %d\n",
            sum(hqstat == "Pass"), sum(hqstat == "Fail"), sum(hqstat == "NA")))
cat(sprintf("  QC     Pass / Fail / NA   : %d / %d / %d\n",
            sum(qc == "Pass"), sum(qc == "Fail"), sum(qc == "NA")))
cat(sprintf("  absent from NR/NV matrices: %d\n", sum(!in_mat)))
