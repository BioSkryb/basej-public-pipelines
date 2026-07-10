#!/usr/bin/env Rscript
# compare_hqstat_qc_redundancy.R
# ---------------------------------------------------------------------------
# Is the HQ-statistical filter redundant with the pileup-QC filter?
# Treats each as an independent per-variant classifier over the candidate universe
# (all variants in scheme_membership):
#   HQStat-pass : Scheme_HQStat == 1        (singleton OR Sequoia-retained shared)
#   QC-pass     : Pileup_Verdict == "Pass"  (from the variant master filter table)
# Cross-tabulates them into a 2x2 contingency with overlap metrics, and writes:
#   hqstat_vs_qc_contingency_<group>.tsv  (metric/value)
#   hqstat_vs_qc_redundancy_<group>.md    (2x2 table + verdict)
#   hqstat_vs_qc_heatmap_<group>.pdf/.png (2x2 heatmap)
#
# Usage: Rscript compare_hqstat_qc_redundancy.R <group> <master_table.tsv> <scheme_membership.tsv>

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
if (file.exists("/usr/local/bin/theme_ohchibi_pubr.R")) {
  source("/usr/local/bin/theme_ohchibi_pubr.R")
} else {
  theme_ohchibi_pubr <- ggplot2::theme_bw
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: compare_hqstat_qc_redundancy.R <group> <master_table> <scheme_membership>")
group <- args[1]; master_f <- args[2]; memb_f <- args[3]

# ── QC-pass set (Pileup_Verdict == "Pass") from the master table ─────────────
master <- fread(master_f, sep = "\t", header = TRUE, data.table = FALSE,
                select = c("VariantId", "Pileup_Verdict"))
if (!all(c("VariantId", "Pileup_Verdict") %in% colnames(master)))
  stop("master table missing VariantId / Pileup_Verdict")
qc_set <- unique(master$VariantId[!is.na(master$Pileup_Verdict) & master$Pileup_Verdict == "Pass"])

# ── HQStat membership from scheme_membership ─────────────────────────────────
memb <- fread(memb_f, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("VariantId", "Scheme_HQStat") %in% colnames(memb)))
  stop("scheme_membership missing VariantId / Scheme_HQStat")
vid <- memb$VariantId
hq  <- memb$Scheme_HQStat >= 1
qc  <- vid %in% qc_set

# ── 2x2 contingency + overlap metrics ────────────────────────────────────────
a <- sum(hq & qc); b <- sum(hq & !qc); cc <- sum(!hq & qc); dd <- sum(!hq & !qc)
n_total <- length(vid); hq_pass <- a + b; qc_pass <- a + cc; uni <- a + b + cc
jacc     <- if (uni > 0)     a / uni     else 0
qc_of_hq <- if (hq_pass > 0) a / hq_pass else 0   # frac of HQStat survivors also kept by QC
hq_of_qc <- if (qc_pass > 0) a / qc_pass else 0   # frac of QC survivors also kept by HQStat

ct <- data.frame(
  metric = c("n_variants", "HQStat_pass", "QC_pass", "both_pass(a)", "HQStat_only(b)",
             "QC_only(c)", "both_fail(d)", "jaccard_HQStat_QC",
             "frac_HQStat_kept_by_QC", "frac_QC_kept_by_HQStat"),
  value  = c(n_total, hq_pass, qc_pass, a, b, cc, dd,
             sprintf("%.4f", jacc), sprintf("%.4f", qc_of_hq), sprintf("%.4f", hq_of_qc)),
  stringsAsFactors = FALSE)
write.table(ct, sprintf("hqstat_vs_qc_contingency_%s.tsv", group),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── verdict ───────────────────────────────────────────────────────────────────
verdict <- if (b == 0 && cc == 0) "Fully redundant: HQStat and QC keep exactly the same variants." else
           if (b == 0)            "QC is redundant after HQStat: QC drops none of HQStat's survivors (every HQStat-pass variant is also QC-pass)." else
           if (cc == 0)           "HQStat is redundant after QC: every QC-pass variant is also HQStat-pass (HQStat drops nothing QC keeps)." else
                                  "Complementary: the two filters drop different variants."
short <- if (b == 0 && cc == 0) "Fully redundant" else
         if (b == 0)            "QC redundant after HQStat" else
         if (cc == 0)           "HQStat redundant after QC" else
                                "Complementary"

md <- c(
  sprintf("# HQStatistical vs QC filter redundancy - %s", group), "",
  sprintf("Universe: %d candidate variants.", n_total), "",
  "| | QC pass | QC fail | total |", "|---|---|---|---|",
  sprintf("| **HQStat pass** | %d | %d | %d |", a, b, hq_pass),
  sprintf("| **HQStat fail** | %d | %d | %d |", cc, dd, cc + dd),
  sprintf("| **total** | %d | %d | %d |", qc_pass, b + dd, n_total), "",
  sprintf("- Jaccard(HQStat, QC) = %.4f", jacc),
  sprintf("- HQStat survivors also kept by QC: %.1f%% (QC drops %d of %d)", qc_of_hq * 100, b, hq_pass),
  sprintf("- QC survivors also kept by HQStat: %.1f%% (HQStat drops %d of %d)", hq_of_qc * 100, cc, qc_pass), "",
  sprintf("**Verdict:** %s", verdict))
writeLines(md, sprintf("hqstat_vs_qc_redundancy_%s.md", group))

# ── 2x2 heatmap ───────────────────────────────────────────────────────────────
df <- data.frame(
  HQStat = factor(c("pass", "pass", "fail", "fail"), levels = c("fail", "pass")),
  QC     = factor(c("pass", "fail", "pass", "fail"), levels = c("pass", "fail")),
  n      = c(a, b, cc, dd))
p <- ggplot(df, aes(QC, HQStat, fill = n)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = format(n, big.mark = ",")), fontface = "bold", size = 6) +
  scale_fill_gradient(low = "#EEF3FB", high = "#2166AC") +
  labs(title = sprintf("HQStatistical vs QC filter  |  %s", group),
       subtitle = sprintf("Jaccard = %.2f   |   QC-only drops %d   ·   HQStat-only drops %d", jacc, b, cc),
       caption = sprintf("Verdict: %s", short),
       x = "QC  (Pileup_Verdict == Pass)", y = "HQStatistical filter", fill = "variants") +
  theme_ohchibi_pubr()
ggsave(sprintf("hqstat_vs_qc_heatmap_%s.png", group), p, width = 8.5, height = 5.2, dpi = 200, bg = "white")
ggsave(sprintf("hqstat_vs_qc_heatmap_%s.pdf", group), p, width = 8.5, height = 5.2)

cat(sprintf("[compare_hqstat_qc_redundancy] %s: a=%d b=%d c=%d d=%d  jaccard=%.4f  [%s]\n",
            group, a, b, cc, dd, jacc, short))
