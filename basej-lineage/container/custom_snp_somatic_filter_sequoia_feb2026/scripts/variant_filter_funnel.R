#!/usr/bin/env Rscript
## Variant filter funnel: tracking table + markdown report + funnel plot.
## Split out of master_table.R so it can run AFTER VAF_SPLIT_VARIANTS_HEXBIN.
##
## Continuous 8-stage funnel:
##   MERGE -> VEP -> Bulk -> Binom -> HQ Statistical (ref) -> +Depth -> +Pileup QC -> +Phylogeny coverage
## Stages 1-4 are upstream per-variant filters (read from the master filter table);
## stages 5-8 are the VAF-split phylogeny-matrix cascade (read from VAF_SPLIT binary
## matrices), computed on the *unfiltered* NR/NV matrices (a different universe from
## Binom: includes singletons, excludes no-evidence variants). Stage 4 -> 5 is therefore
## a universe change, not a strict subset. HQ Statistical reports its TRUE input = the
## unfiltered matrix count; its removed = shared variants dropped by 2nd-pass Sequoia.
##
## Usage: Rscript variant_filter_funnel.R <group_id>
## Reads from the current (Nextflow work) dir:
##   variant_master_filter_table_<group>.tsv   (VEP / Bulk / Binom stages)
##   all_variants_*.txt                         (n_total reference count)
##   <group>_binary_matrix_HQRoundStatisticalFiltered.tsv
##   <group>_binary_matrix_HQRoundStatisticalFilteredPlusQCFiltered.tsv
##   <group>_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered.tsv
##   <group>_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv
## Writes:
##   variant_filter_tracking_<group>.tsv
##   variant_filter_report_<group>.md
##   variant_filter_plot_<group>.pdf

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript variant_filter_funnel.R <group_id>")
group_id <- args[1]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
# Shared ggplot theme (baked into the container); fall back to theme_bw.
if (file.exists("/usr/local/bin/theme_ohchibi_pubr.R")) {
  source("/usr/local/bin/theme_ohchibi_pubr.R")
} else {
  theme_ohchibi_pubr <- ggplot2::theme_bw
}

master_file     <- paste0("variant_master_filter_table_", group_id, ".tsv")
merged_vcf_file <- list.files(".", pattern = "^all_variants_.*\\.txt$")[1]
if (!file.exists(master_file)) stop("Master table not found: ", master_file)
if (is.na(merged_vcf_file))    stop("all_variants_*.txt (n_total reference) not found")

master <- fread(master_file, sep = "\t", header = TRUE, na.strings = "NA",
                data.table = FALSE, showProgress = FALSE)

pct <- function(n, d) if (is.na(d) || d == 0) NA else round(n / d * 100, 2)
fmt_pct <- function(n, d) { p <- pct(n, d); if (is.na(p)) "NA" else sprintf("%.2f", p) }
# Count data rows (variants) of a binary matrix without loading all cell columns.
count_rows <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NA_integer_)
  as.integer(nrow(fread(path, sep = "\t", header = TRUE, select = 1L,
                        data.table = FALSE, showProgress = FALSE)))
}
find1 <- function(pattern) { f <- list.files(".", pattern = pattern); if (length(f)) f[1] else NA_character_ }

# ── n_total: variants in the merged VCF (reference universe) ──────────────────
merged_vcf <- readLines(merged_vcf_file)
merged_vcf <- merged_vcf[nzchar(merged_vcf)]
n_total    <- length(merged_vcf)

# ── A. Upstream master-table stages (MERGE -> VEP -> Bulk -> Binom) ───────────
n_vep_in     <- n_total
n_vep_pass   <- sum(master$VEP_FILTER_STATUS == "PASS", na.rm = TRUE)
n_bulk_in    <- n_vep_pass
n_bulk_pass  <- sum(master$Bulk_RemainingAfterBulk == "Pass", na.rm = TRUE)
n_binom_in   <- sum(!is.na(master$Binom_BinomialBetabinomialFilter))
n_binom_pass <- sum(master$Binom_BinomialBetabinomialFilter == "Pass", na.rm = TRUE)

# ── B. VAF-split cascade stages (HQ Statistical ref -> +QC -> +Depth -> +Phylo) ─
g <- group_id
n_vaf_unfilt <- count_rows(find1(paste0("^", g, "_binary_matrix_unfiltered\\.tsv$")))
n_vaf_hqstat <- count_rows(find1(paste0("^", g, "_binary_matrix_HQRoundStatisticalFiltered\\.tsv$")))
n_vaf_qc     <- count_rows(find1(paste0("^", g, "_binary_matrix_HQRoundStatisticalFilteredPlusQCFiltered\\.tsv$")))
n_vaf_depth  <- count_rows(find1(paste0("^", g, "_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered\\.tsv$")))
n_vaf_phylo  <- count_rows(find1(paste0("^", g, "_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny\\.tsv$")))
have_vaf <- !any(is.na(c(n_vaf_hqstat, n_vaf_depth, n_vaf_qc, n_vaf_phylo)))
if (!have_vaf)
  cat("WARNING: one or more VAF-split binary matrices not found; VAF-split stages will be NA.\n")
# Fallback if the unfiltered matrix is unavailable: use HQStat count (old reset behaviour).
if (is.na(n_vaf_unfilt)) n_vaf_unfilt <- n_vaf_hqstat

# ── C. Tracking table (8 continuous stages) ───────────────────────────────────
# The VAF-split cascade (stages 5-8) is computed on the *unfiltered* NR/NV matrices,
# a different universe from the Binom stage (it includes singletons and excludes
# no-evidence variants), so stage 4 -> 5 is a universe change, not a strict subset.
# HQ Statistical now reports its TRUE input = the unfiltered matrix (singletons + shared);
# n_in - n_pass is the count of shared variants dropped by the 2nd-pass Sequoia (Rho/qval).
stages <- list(
  list(stage = "MERGE_PROCESSED_VCF",     n_in = n_total,      n_pass = n_total),
  list(stage = "VEP_Germline_Filter",     n_in = n_vep_in,     n_pass = n_vep_pass),
  list(stage = "Bulk_Filter",             n_in = n_bulk_in,    n_pass = n_bulk_pass),
  list(stage = "Binom_Betabinom_Filter",  n_in = n_binom_in,   n_pass = n_binom_pass),
  list(stage = "VAFSplit_HQStatistical",  n_in = n_vaf_unfilt, n_pass = n_vaf_hqstat),
  list(stage = "VAFSplit_QCFiltered",     n_in = n_vaf_hqstat, n_pass = n_vaf_qc),
  list(stage = "VAFSplit_DepthFiltered",  n_in = n_vaf_qc,     n_pass = n_vaf_depth),
  list(stage = "VAFSplit_ForPhylogeny",   n_in = n_vaf_depth,  n_pass = n_vaf_phylo)
)

tracking <- do.call(rbind, lapply(stages, function(s) {
  data.frame(
    Stage                  = s$stage,
    N_Input                = s$n_in,
    N_Passing              = s$n_pass,
    N_Filtered_Out         = ifelse(is.na(s$n_in) | is.na(s$n_pass), NA, s$n_in - s$n_pass),
    Pct_Retained_vs_Input  = pct(s$n_pass, s$n_in),
    Pct_Retained_vs_Total  = pct(s$n_pass, n_total),
    stringsAsFactors       = FALSE
  )
}))

write.table(tracking, paste0("variant_filter_tracking_", group_id, ".tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Filter tracking table:\n"); print(tracking)

# ── D. Markdown report ────────────────────────────────────────────────────────
bulk_pass_ids <- master$VariantId[!is.na(master$Bulk_RemainingAfterBulk) & master$Bulk_RemainingAfterBulk == "Pass"]
binom_all_ids <- master$VariantId[!is.na(master$Binom_BinomialBetabinomialFilter)]
n_vcf_no_bam  <- length(setdiff(bulk_pass_ids, binom_all_ids))

no_bulk_msg <- if (!is.na(n_bulk_pass) && n_bulk_pass == n_bulk_in) {
  paste0("Variants found in a matched bulk (non-single-cell) sample are flagged as germline",
         " and removed. In this run **no bulk sample was provided**, so all ", n_bulk_pass,
         " variants pass this stage unchanged.")
} else {
  paste0("Variants found in a matched bulk (non-single-cell) sample are flagged as germline",
         " and removed. **", n_bulk_in - n_bulk_pass, "** variants were identified as",
         " bulk germline and discarded.")
}

md <- c(
  paste0("# Variant Filter Tracking Report — Group: ", group_id),
  "",
  paste0("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "---",
  "",
  "## Overview",
  "",
  paste0("This report tracks how many variants are retained or removed at each filtering",
         " stage of the somatic SNP/indel pipeline for group **", group_id, "**, as one",
         " continuous funnel. Stages 1–4 are the upstream per-variant filters (read from the",
         " master filter table). Stages 5–8 are the VAF-split phylogeny-matrix cascade, computed",
         " on the *unfiltered* NR/NV matrices — a different universe from the Binomial stage",
         " (includes singletons, excludes no-evidence variants), so stage 4 → 5 is a universe",
         " change, not a strict subset. HQ Statistical's input is the unfiltered matrix count."),
  "",
  "---",
  "",
  "## Filter Funnel",
  "",
  "```",
  sprintf("  %-34s  %8s  %10s  %12s  %s",
          "Stage", "N Input", "N Passing", "N Removed", "% Kept (of input)"),
  sprintf("  %-34s  %8s  %10s  %12s  %s",
          strrep("-", 34), strrep("-", 8), strrep("-", 10), strrep("-", 12), strrep("-", 17)),
  sprintf("  %-34s  %8d  %10d  %12d  %s",   "MERGE_PROCESSED_VCF (reference)", n_total, n_total, 0L, "100.00 %"),
  sprintf("  %-34s  %8d  %10d  %12d  %s %%", "VEP Germline Filter",            n_vep_in,   n_vep_pass,  n_vep_in   - n_vep_pass,  fmt_pct(n_vep_pass,  n_vep_in)),
  sprintf("  %-34s  %8d  %10d  %12d  %s %%", "Bulk Germline Filter",           n_bulk_in,  n_bulk_pass, n_bulk_in  - n_bulk_pass, fmt_pct(n_bulk_pass, n_bulk_in)),
  sprintf("  %-34s  %8d  %10d  %12d  %s %%", "Binomial / Betabinomial Filter", n_binom_in, n_binom_pass,n_binom_in - n_binom_pass,fmt_pct(n_binom_pass,n_binom_in)),
  if (have_vaf) sprintf("  %-34s  %8d  %10d  %12d  %s %%", "HQ Statistical (vs unfiltered)", n_vaf_unfilt, n_vaf_hqstat, n_vaf_unfilt - n_vaf_hqstat, fmt_pct(n_vaf_hqstat, n_vaf_unfilt)) else "  HQ Statistical                      (VAF-split matrices not found)",
  if (have_vaf) sprintf("  %-34s  %8d  %10d  %12d  %s %%", "+ Pileup QC (Verdict==Pass)",   n_vaf_hqstat, n_vaf_qc,     n_vaf_hqstat - n_vaf_qc,    fmt_pct(n_vaf_qc,    n_vaf_hqstat)) else NULL,
  if (have_vaf) sprintf("  %-34s  %8d  %10d  %12d  %s %%", "+ Depth Filter",                n_vaf_qc,     n_vaf_depth,  n_vaf_qc     - n_vaf_depth, fmt_pct(n_vaf_depth, n_vaf_qc)) else NULL,
  if (have_vaf) sprintf("  %-34s  %8d  %10d  %12d  %s %%", "+ Phylogeny coverage",          n_vaf_depth,  n_vaf_phylo,  n_vaf_depth  - n_vaf_phylo, fmt_pct(n_vaf_phylo, n_vaf_depth)) else NULL,
  "```",
  "",
  if (have_vaf) paste0("**Final output:** ", n_vaf_phylo, " phylogeny-ready variants out of ",
         n_total, " merged candidates (**", fmt_pct(n_vaf_phylo, n_total), " %** of the original set).") else
        paste0("**Upstream output:** ", n_binom_pass, " variants pass the Binomial/Betabinomial filter."),
  "",
  "---",
  "",
  "## Stage-by-Stage Explanation",
  "",
  paste0("### 1. MERGE_PROCESSED_VCF — ", n_total, " variants"),
  "",
  paste0("The starting universe: all variants detected and merged across all single cells.",
         " Includes germline variants, technical artefacts, and true somatic mutations."),
  "",
  paste0("### 2. VEP Germline Filter — ", n_vep_pass, " / ", n_vep_in, " pass (", fmt_pct(n_vep_pass, n_vep_in), " %)"),
  "",
  paste0("Variants cross-referenced against population databases via VEP; removed if present",
         " in a known germline database above the configured allele-frequency threshold (or",
         " known existing variation). Removes **", n_vep_in - n_vep_pass, "** variants."),
  "",
  paste0("### 3. Bulk Germline Filter — ", n_bulk_pass, " / ", n_bulk_in, " pass (", fmt_pct(n_bulk_pass, n_bulk_in), " %)"),
  "",
  no_bulk_msg,
  "",
  paste0("### 4. Binomial / Betabinomial Filter — ", n_binom_pass, " / ", n_binom_in, " pass (", fmt_pct(n_binom_pass, n_binom_in), " %)"),
  "",
  paste0("Each variant is tested across all single cells for a germline (binomial) vs somatic",
         " mosaic (betabinomial) allele-frequency pattern; germline-consistent variants are",
         " removed. Removes **", n_binom_in - n_binom_pass, "** variants."),
  "",
  if (n_vcf_no_bam > 0) paste0("> **Gap — ", n_vcf_no_bam, " variants absent from this stage:**  \n",
         "> ", n_vcf_no_bam, " variants pass the bulk filter but never reach the",
         " binomial/betabinomial step — called in the VCF but with no read evidence in any",
         " single-cell BAM (PresentVCF_NotInBAM), so absent from the pileup count matrices.") else NULL,
  "",
  "### 5–8. VAF-split phylogeny-matrix cascade",
  "",
  if (have_vaf) paste0("Built by VAF_SPLIT_VARIANTS_HEXBIN on the *unfiltered* NR/NV matrices",
         " (shared/singleton split + anchor-aware binarization). The cascade starts from the",
         " **unfiltered** matrix (**", n_vaf_unfilt, "** evidence-bearing variants); **Stage 5",
         " (HQ Statistical)** keeps **", n_vaf_hqstat, "** (singletons + shared retained by the",
         " 2nd-pass Sequoia Rho/qval test, dropping ", n_vaf_unfilt - n_vaf_hqstat, "). The cascade then applies the",
         " **pileup-QC** (Verdict==Pass) filter (→ ", n_vaf_qc, "), the **depth filter** (→ ",
         n_vaf_depth, "), and the **≥70%-cell phylogeny-coverage** filter (→ ", n_vaf_phylo,
         " phylogeny-ready variants).") else
        "VAF-split binary matrices were not found in the work dir; stages 5–8 are reported as NA.",
  "",
  "---",
  "",
  "## Summary Table",
  "",
  "| Stage | Variants In | Variants Out | % Retained |",
  "|---|---|---|---|",
  paste0("| Merged VCF | — | ", n_total, " | 100.00 % |"),
  paste0("| VEP Germline Filter | ", n_vep_in, " | ", n_vep_pass, " | ", fmt_pct(n_vep_pass, n_vep_in), " % |"),
  paste0("| Bulk Germline Filter | ", n_bulk_in, " | ", n_bulk_pass, " | ", fmt_pct(n_bulk_pass, n_bulk_in), " % |"),
  paste0("| Binomial/Betabinomial Filter | ", n_binom_in, " | ", n_binom_pass, " | ", fmt_pct(n_binom_pass, n_binom_in), " % |"),
  if (have_vaf) paste0("| HQ Statistical (vs unfiltered) | ", n_vaf_unfilt, " | ", n_vaf_hqstat, " | ", fmt_pct(n_vaf_hqstat, n_vaf_unfilt), " % |") else NULL,
  if (have_vaf) paste0("| + Pileup QC | ", n_vaf_hqstat, " | ", n_vaf_qc, " | ", fmt_pct(n_vaf_qc, n_vaf_hqstat), " % |") else NULL,
  if (have_vaf) paste0("| + Depth | ", n_vaf_qc, " | ", n_vaf_depth, " | ", fmt_pct(n_vaf_depth, n_vaf_qc), " % |") else NULL,
  if (have_vaf) paste0("| + Phylogeny coverage | ", n_vaf_depth, " | ", n_vaf_phylo, " | ", fmt_pct(n_vaf_phylo, n_vaf_depth), " % |") else NULL,
  ""
)
md <- md[!vapply(md, is.null, logical(1))]
writeLines(unlist(md), paste0("variant_filter_report_", group_id, ".md"))
cat("Markdown report written.\n")

# ── E. Funnel plot — horizontal waterfall (ggplot, log scale; PDF + PNG) ──────
# One horizontal bar per stage (top -> bottom). The solid bar is N_Passing
# (coloured by phase); the faded grey extension to N_Input shows the variants
# removed at that stage. Each row is annotated with N and "-removed (% kept)".
plot_pdf <- paste0("variant_filter_plot_", group_id, ".pdf")
plot_png <- paste0("variant_filter_plot_", group_id, ".png")

fd <- as.data.table(tracking)
fd <- fd[!is.na(N_Passing) & !is.na(N_Input)]      # drop VAF rows when have_vaf is FALSE
nice <- c(MERGE_PROCESSED_VCF = "Merged VCF",          VEP_Germline_Filter = "VEP germline",
          Bulk_Filter = "Bulk germline",               Binom_Betabinom_Filter = "Binom / Betabinom",
          VAFSplit_HQStatistical = "HQ Statistical", VAFSplit_QCFiltered = "+ Pileup QC",
          VAFSplit_DepthFiltered = "+ Depth",          VAFSplit_ForPhylogeny = "+ Phylogeny")
fd[, label := ifelse(Stage %in% names(nice), nice[Stage], Stage)]
fd[, phase := ifelse(grepl("^VAFSplit", Stage), "VAF-split cascade", "Upstream filters")]
fd[, lab   := factor(label, levels = rev(label))]   # preserve order, first stage on top
pal <- c("Upstream filters" = "#2166AC", "VAF-split cascade" = "#1B7837")

p_funnel <- ggplot(fd, aes(y = lab)) +
  geom_col(aes(x = N_Input), fill = "grey88", width = 0.64) +
  geom_col(aes(x = N_Passing, fill = phase), width = 0.64) +
  geom_text(aes(x = N_Passing, label = comma(N_Passing)),
            hjust = 1.12, size = 3.0, fontface = "bold", colour = "white") +
  geom_text(data = fd[N_Filtered_Out > 0],
            aes(x = N_Input, label = paste0("-", comma(N_Filtered_Out), "  (",
                sprintf("%.1f%%", Pct_Retained_vs_Input), " kept)")),
            hjust = -0.08, size = 2.6, colour = "grey35") +
  scale_x_log10(labels = comma, expand = expansion(mult = c(0, 0.42))) +
  scale_fill_manual(values = pal) +
  labs(title = paste0("Variant filtering funnel  |  Group: ", group_id),
       subtitle = sprintf("Final: %s phylogeny-ready of %s merged (%.2f%%)",
                          comma(tail(fd$N_Passing, 1)), comma(n_total),
                          tail(fd$Pct_Retained_vs_Total, 1)),
       x = "Variants (log scale)", y = NULL, fill = NULL) +
  theme_ohchibi_pubr() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

ggsave(plot_pdf, p_funnel, width = 9, height = 5.2, device = cairo_pdf)
ggsave(plot_png, p_funnel, width = 9, height = 5.2, dpi = 200, bg = "white")
cat("Filter plot written:", plot_pdf, "and", plot_png, "\n")
