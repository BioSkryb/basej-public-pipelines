## VAF (NV/NR) per variant: MEAN_VAF, MEDIAN_VAF, then hexbin vs -log10(q) and vs Rho.
## Called from CUSTOM_VARIANT_FILTER_PROVENANCE (mat_nv/mat_nr files staged in work dir).
## Usage: Rscript vaf_hexbin_plots.R <master_table.tsv> <group_id> <output.pdf>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript vaf_hexbin_plots.R <master_table.tsv> <group_id> <output.pdf>")
master_tsv  <- args[1]
group_id    <- args[2]
output_pdf  <- args[3]

library(ggplot2)
library(paletteer)
library(scales)

## ── Load master table for Binom_Germline_qval_log10 and Binom_Rho ─────────────
master <- read.table(master_tsv, sep = "\t", header = TRUE, quote = "", comment.char = "")

## ── Find mat_nv and mat_nr files (staged by Nextflow) ─────────────────────────
nv_files <- list.files(".", pattern = paste0("^mat_nv_group_", group_id, "_chr_.*\\.tsv$"), full.names = TRUE)
nr_files <- list.files(".", pattern = paste0("^mat_nr_group_", group_id, "_chr_.*\\.tsv$"), full.names = TRUE)
if (length(nv_files) == 0 || length(nr_files) == 0) {
  stop("No mat_nv/mat_nr files found for group ", group_id, " in current directory.")
}
## Match by chr suffix (e.g. mat_nv_group_Patient_chr_1 <-> mat_nr_group_Patient_chr_1)
nv_bases <- sub("\\.tsv$", "", basename(nv_files))
nr_bases <- sub("\\.tsv$", "", basename(nr_files))
chr_nv   <- sub(paste0("^mat_nv_group_", group_id, "_chr_"), "", nv_bases)
chr_nr   <- sub(paste0("^mat_nr_group_", group_id, "_chr_"), "", nr_bases)
common_chr <- intersect(chr_nv, chr_nr)
if (length(common_chr) == 0) stop("No matching mat_nv/mat_nr pairs found.")
nv_files <- nv_files[match(common_chr, chr_nv)]
nr_files <- nr_files[match(common_chr, chr_nr)]

## ── Compute per-variant VAF metrics ───────────────────────────────────────────
## POOLED_VAF   = sum(NV across cells) / sum(NR across cells)
##                Treats all cells as one pooled sample; low for somatic
##                (few cells contribute NV) vs ~0.5 for germline (all cells do).
## CELL_FRAC    = fraction of covered cells (NR > 0) where NV/NR > 0.1
##                Captures how many cells carry the variant, regardless of depth.
CELL_FRAC_THRESH <- 0.1
vaf_list <- vector("list", length(nv_files))
for (i in seq_along(nv_files)) {
  nv <- read.table(nv_files[i], sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  nr <- read.table(nr_files[i], sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(identical(nv[, 1], nr[, 1]))
  variant_ids <- nv[, 1]
  nv_mat <- as.matrix(nv[, -1, drop = FALSE])
  nr_mat <- as.matrix(nr[, -1, drop = FALSE])
  covered  <- nr_mat > 0                              # cells with any reads
  vaf_mat  <- nv_mat / nr_mat
  vaf_mat[!covered] <- NA
  pooled_vaf  <- rowSums(nv_mat) / pmax(rowSums(nr_mat), 1)  # pmax avoids /0
  cell_frac   <- rowSums(vaf_mat > CELL_FRAC_THRESH, na.rm = TRUE) /
                   pmax(rowSums(covered), 1)
  vaf_list[[i]] <- data.frame(
    VariantId  = variant_ids,
    POOLED_VAF = pooled_vaf,
    CELL_FRAC  = cell_frac,
    stringsAsFactors = FALSE
  )
}
vaf_df <- do.call(rbind, vaf_list)

## ── Merge with master ─────────────────────────────────────────────────────────
cols_need <- c("VariantId", "Binom_Germline_qval_log10", "Binom_Rho", "Binom_Verdict",
               "VEP_FILTER_STATUS")
cols_need <- intersect(cols_need, colnames(master))
if (!all(c("Binom_Germline_qval_log10", "Binom_Rho") %in% cols_need))
  stop("Master table missing Binom_Germline_qval_log10 or Binom_Rho.")
merge_df <- merge(vaf_df, master[, cols_need], by = "VariantId", all.x = FALSE, all.y = FALSE)

## ── Verdict labels & palette (mirrors binom_ggplots.R) ────────────────────────
verdict_map <- c(
  "TRUE_TRUE_TRUE"   = "Somatic (Pass)",
  "TRUE_FALSE_FALSE" = "Germline: both",
  "TRUE_FALSE_TRUE"  = "Germline: Binom only",
  "TRUE_TRUE_FALSE"  = "Germline: Betabinom only"
)
pal_verdict <- setNames(
  as.character(paletteer_d("RColorBrewer::Set1")[1:4]),
  c("Somatic (Pass)", "Germline: both", "Germline: Binom only", "Germline: Betabinom only")
)
if ("Binom_Verdict" %in% colnames(merge_df)) {
  merge_df$Verdict <- verdict_map[merge_df$Binom_Verdict]
  merge_df$Verdict <- factor(merge_df$Verdict, levels = names(pal_verdict))
}

out_tsv <- paste0("vaf_binom_merged_", group_id, ".tsv")
write.table(merge_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("Written:", out_tsv, nrow(merge_df), "rows\n")

## ── Shared theme ───────────────────────────────────────────────────────────────
theme_clean <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position  = "right",
        plot.title       = element_text(face = "bold", size = 13))

fill_scale <- scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",
                                     name = "Count", trans = "log10",
                                     labels = function(x) formatC(x, format = "d", big.mark = ","))

## ── Plot 1: Pooled VAF vs log10(q-value) ──────────────────────────────────────
sub1 <- merge_df[is.finite(merge_df$POOLED_VAF) & is.finite(merge_df$Binom_Germline_qval_log10), ]
p1 <- ggplot(sub1, aes(x = POOLED_VAF, y = Binom_Germline_qval_log10)) +
  geom_bin2d(bins = 60) +
  fill_scale +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish) +
  scale_y_continuous(limits = c(-10, 0), breaks = seq(-10, 0, 2), oob = squish) +
  labs(title    = paste0("Pooled VAF vs log10(q-value)  |  Group: ", group_id),
       x        = "Pooled VAF  (sum(NV) / sum(NR) across all cells)",
       y        = expression(log[10](q-value)),
       subtitle = paste0("Pooled VAF collapses toward 0 for somatic variants (few cells contribute NV);\n",
                         "germline variants show Pooled VAF ~0.5 across all cells")) +
  theme_clean

## ── Plot 2: Cell Fraction vs log10(q-value) ───────────────────────────────────
sub2 <- merge_df[is.finite(merge_df$CELL_FRAC) & is.finite(merge_df$Binom_Germline_qval_log10), ]
p2 <- ggplot(sub2, aes(x = CELL_FRAC, y = Binom_Germline_qval_log10)) +
  geom_bin2d(bins = 60) +
  fill_scale +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish,
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(-10, 0), breaks = seq(-10, 0, 2), oob = squish) +
  labs(title    = paste0("Cell Fraction vs log10(q-value)  |  Group: ", group_id),
       x        = paste0("Cell Fraction  (cells with VAF > ", CELL_FRAC_THRESH, " / covered cells)"),
       y        = expression(log[10](q-value)),
       subtitle = paste0("Cell fraction ~0% for somatic (1 cell carries variant); ",
                         "~100% for germline (all cells carry variant)")) +
  theme_clean

## ── Plot 3: Pooled VAF vs Rho ──────────────────────────────────────────────────
sub3 <- merge_df[is.finite(merge_df$POOLED_VAF) & is.finite(merge_df$Binom_Rho), ]
p3 <- ggplot(sub3, aes(x = POOLED_VAF, y = Binom_Rho)) +
  geom_bin2d(bins = 60) +
  fill_scale +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), oob = squish) +
  labs(title    = paste0("Pooled VAF vs Rho (overdispersion)  |  Group: ", group_id),
       x        = "Pooled VAF  (sum(NV) / sum(NR) across all cells)",
       y        = "Rho",
       subtitle = "Betabinomial overdispersion; Rho > 0.1 flags high per-cell variability") +
  theme_clean

## ── Plot 4: Cell Fraction vs Pooled VAF (coloured by q-value) ─────────────────
sub4 <- merge_df[is.finite(merge_df$POOLED_VAF) & is.finite(merge_df$CELL_FRAC) &
                   is.finite(merge_df$Binom_Germline_qval_log10), ]
p4 <- ggplot(sub4, aes(x = POOLED_VAF, y = CELL_FRAC, colour = Binom_Germline_qval_log10)) +
  geom_point(alpha = 0.25, size = 0.7) +
  scale_colour_gradientn(colours = rev(paletteer::paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71", 256)),
                         name    = expression(log[10](q)),
                         limits  = c(-10, 0), oob = squish) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish,
                     labels = scales::percent_format(accuracy = 1)) +
  labs(title    = paste0("Cell Fraction vs Pooled VAF  |  Group: ", group_id),
       x        = "Pooled VAF  (sum(NV) / sum(NR) across all cells)",
       y        = paste0("Cell Fraction  (VAF > ", CELL_FRAC_THRESH, ")"),
       subtitle = paste0("Somatic: bottom-left (low pooled VAF, few cells). ",
                         "Germline: top-right (high pooled VAF, many cells). ",
                         "Colour = log10(q-value).")) +
  theme_clean

## ── Verdict-based plots (p5–p7): only produced when Binom_Verdict is present ──
p5 <- p6 <- p7 <- NULL
if ("Verdict" %in% colnames(merge_df)) {

  verdict_scales <- list(
    scale_fill_manual(values  = pal_verdict),
    scale_colour_manual(values = pal_verdict),
    guides(fill   = guide_legend(nrow = 2),
           colour = guide_legend(nrow = 2))
  )
  theme_verdict <- theme_clean +
    theme(legend.position = "bottom", legend.title = element_blank())

  ## ── Plot 5: Cell Fraction vs Pooled VAF, coloured by Verdict ────────────────
  ## Shows WHERE in the (Pooled VAF, Cell Fraction) space each filter category sits.
  sub5 <- merge_df[!is.na(merge_df$Verdict) &
                     is.finite(merge_df$POOLED_VAF) &
                     is.finite(merge_df$CELL_FRAC), ]
  p5 <- ggplot(sub5, aes(x = POOLED_VAF, y = CELL_FRAC, colour = Verdict)) +
    geom_point(alpha = 0.3, size = 0.8) +
    verdict_scales +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish,
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title    = paste0("Cell Fraction vs Pooled VAF by Verdict  |  Group: ", group_id),
         x        = "Pooled VAF  (sum(NV) / sum(NR) across all cells)",
         y        = paste0("Cell Fraction  (VAF > ", CELL_FRAC_THRESH, ")"),
         subtitle = paste0("Binom filter acts on q-value; Betabinom filter acts on Rho. ",
                           "Both should separate somatic (bottom-left) from germline (top-right).")) +
    theme_verdict

  ## ── Plot 6: Density of Pooled VAF by Verdict ────────────────────────────────
  ## Mirrors binom_ggplots.R p1 (q-value density) but uses Pooled VAF.
  ## If Pooled VAF captures the binom signal, somatic peak should be near 0
  ## and Germline: Binom only peak near 0.5.
  sub6 <- merge_df[!is.na(merge_df$Verdict) & is.finite(merge_df$POOLED_VAF), ]
  p6 <- ggplot(sub6, aes(x = POOLED_VAF, fill = Verdict, colour = Verdict)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    verdict_scales +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), oob = squish) +
    labs(title    = paste0("Pooled VAF Distribution by Verdict  |  Group: ", group_id),
         x        = "Pooled VAF  (sum(NV) / sum(NR) across all cells)",
         y        = "Density",
         subtitle = paste0("Somatic (Pass) should peak near 0 (few cells carry NV); ",
                           "Germline variants should peak near 0.5")) +
    theme_verdict

  ## ── Plot 7: Density of Cell Fraction by Verdict ─────────────────────────────
  ## Mirrors binom_ggplots.R p2 (Rho density) but uses Cell Fraction.
  ## If Cell Fraction captures the betabinom signal, high-Rho (germline) variants
  ## should show high cell fraction and somatic variants low cell fraction.
  sub7 <- merge_df[!is.na(merge_df$Verdict) & is.finite(merge_df$CELL_FRAC), ]
  p7 <- ggplot(sub7, aes(x = CELL_FRAC, fill = Verdict, colour = Verdict)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    verdict_scales +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), oob = squish,
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title    = paste0("Cell Fraction Distribution by Verdict  |  Group: ", group_id),
         x        = paste0("Cell Fraction  (cells with VAF > ", CELL_FRAC_THRESH, " / covered cells)"),
         y        = "Density",
         subtitle = paste0("Somatic (Pass) should peak near 0%; ",
                           "Germline variants should peak near 100%")) +
    theme_verdict
}

## ── VEP-based plots (p8–p12): only for VEP-assessed variants (PASS or REMOVED) ─
## "Not assessed" (NA VEP_FILTER_STATUS) are excluded — these are variants absent
## from population databases, i.e. predominantly somatic candidates.
## Focus: do VEP: Removed variants look germline by Pooled VAF, Cell Fraction, Rho, q-value?
p8 <- p9 <- p10 <- p11 <- p12 <- NULL
if ("VEP_FILTER_STATUS" %in% colnames(merge_df)) {
  vep_assessed <- merge_df[!is.na(merge_df$VEP_FILTER_STATUS), ]
  vep_assessed$VEP_Status <- ifelse(vep_assessed$VEP_FILTER_STATUS == "PASS",
                                    "VEP: Pass", "VEP: Removed")
  vep_assessed$VEP_Status <- factor(vep_assessed$VEP_Status,
                                     levels = c("VEP: Removed", "VEP: Pass"))
  pal_vep <- c("VEP: Removed" = "#E41A1C", "VEP: Pass" = "#4DAF4A")
  vep_colour_scale <- scale_colour_manual(values = pal_vep, name = "VEP status")
  vep_fill_scale   <- scale_fill_manual(values   = pal_vep, name = "VEP status")
  theme_vep <- theme_clean +
    theme(legend.position = "bottom", legend.title = element_blank())
  vep_guides <- guides(fill   = guide_legend(nrow = 1),
                       colour = guide_legend(nrow = 1))

  n_removed <- sum(vep_assessed$VEP_Status == "VEP: Removed")
  n_pass    <- sum(vep_assessed$VEP_Status == "VEP: Pass")
  vep_n_label <- paste0("VEP: Removed n=", n_removed, "  |  VEP: Pass n=", n_pass)

  ## ── Plot 8: Cell Fraction vs Pooled VAF, coloured by VEP status ─────────────
  sub8 <- vep_assessed[is.finite(vep_assessed$POOLED_VAF) & is.finite(vep_assessed$CELL_FRAC), ]
  p8 <- ggplot(sub8, aes(x = POOLED_VAF, y = CELL_FRAC, colour = VEP_Status)) +
    geom_point(alpha = 0.4, size = 0.9) +
    vep_colour_scale +
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 2))) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), oob = squish,
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title    = paste0("Cell Fraction vs Pooled VAF by VEP Status  |  Group: ", group_id),
         x        = "Pooled VAF  (sum(NV) / sum(NR) across all cells)",
         y        = paste0("Cell Fraction  (VAF > ", CELL_FRAC_THRESH, ")"),
         subtitle = paste0(vep_n_label, "\n",
                           "VEP: Removed = known germline by population AF (expect top-right = high VAF, all cells).")) +
    theme_vep

  ## ── Plot 9: Pooled VAF density by VEP status ────────────────────────────────
  sub9 <- vep_assessed[is.finite(vep_assessed$POOLED_VAF), ]
  p9 <- ggplot(sub9, aes(x = POOLED_VAF, fill = VEP_Status, colour = VEP_Status)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    vep_fill_scale + vep_colour_scale + vep_guides +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), oob = squish) +
    labs(title    = paste0("Pooled VAF Distribution by VEP Status  |  Group: ", group_id),
         x        = "Pooled VAF  (sum(NV) / sum(NR) across all cells)",
         y        = "Density",
         subtitle = paste0(vep_n_label, "\n",
                           "VEP: Removed should peak near 0.5 (heterozygous germline)")) +
    theme_vep

  ## ── Plot 10: Cell Fraction density by VEP status ────────────────────────────
  sub10 <- vep_assessed[is.finite(vep_assessed$CELL_FRAC), ]
  p10 <- ggplot(sub10, aes(x = CELL_FRAC, fill = VEP_Status, colour = VEP_Status)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    vep_fill_scale + vep_colour_scale + vep_guides +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), oob = squish,
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title    = paste0("Cell Fraction Distribution by VEP Status  |  Group: ", group_id),
         x        = paste0("Cell Fraction  (cells with VAF > ", CELL_FRAC_THRESH, " / covered cells)"),
         y        = "Density",
         subtitle = paste0(vep_n_label, "\n",
                           "VEP: Removed should peak near 100% (all cells carry germline allele)")) +
    theme_vep

  ## ── Plot 11: log10(q-value) density by VEP status ───────────────────────────
  ## Validates concordance: VEP-removed (germline by pop. AF) should also have
  ## q-value near 0 (binomial test agrees VAF ~0.5); VEP: Pass should be more spread.
  sub11 <- vep_assessed[is.finite(vep_assessed$Binom_Germline_qval_log10), ]
  p11 <- ggplot(sub11, aes(x = Binom_Germline_qval_log10, fill = VEP_Status, colour = VEP_Status)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    vep_fill_scale + vep_colour_scale + vep_guides +
    scale_x_continuous(limits = c(-10, 0), oob = squish) +
    labs(title    = paste0("Binomial q-value Distribution by VEP Status  |  Group: ", group_id),
         x        = expression(log[10](q-value)),
         y        = "Density",
         subtitle = paste0(vep_n_label, "\n",
                           "VEP: Removed should cluster near 0 (germline VAF ~0.5 confirmed by binom test)")) +
    theme_vep

  ## ── Plot 12: Rho (overdispersion) density by VEP status ─────────────────────
  ## Betabinomial Rho captures per-cell variability of the allele fraction.
  ## Germline variants (VEP: Removed) carried by all cells should show a distinct
  ## Rho profile vs VEP: Pass variants (mixed somatic + germline survivors).
  sub12 <- vep_assessed[is.finite(vep_assessed$Binom_Rho), ]
  p12 <- ggplot(sub12, aes(x = Binom_Rho, fill = VEP_Status, colour = VEP_Status)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    vep_fill_scale + vep_colour_scale + vep_guides +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), oob = squish) +
    geom_vline(xintercept = 0.1, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    annotate("text", x = 0.12, y = Inf, label = "Rho = 0.1\n(betabinom cutoff)",
             vjust = 1.3, hjust = 0, colour = "grey40", size = 3) +
    labs(title    = paste0("Rho (Overdispersion) Distribution by VEP Status  |  Group: ", group_id),
         x        = "Rho  (betabinomial overdispersion)",
         y        = "Density",
         subtitle = paste0(vep_n_label, "\n",
                           "Rho > 0.1 flags high per-cell variability; dashed line = betabinom filter threshold")) +
    theme_vep
}

pdf(output_pdf, width = 12, height = 7)
print(p1)
print(p2)
print(p3)
print(p4)
if (!is.null(p5)) print(p5)
if (!is.null(p6)) print(p6)
if (!is.null(p7)) print(p7)
if (!is.null(p8))  print(p8)
if (!is.null(p9))  print(p9)
if (!is.null(p10)) print(p10)
if (!is.null(p11)) print(p11)
if (!is.null(p12)) print(p12)
dev.off()
cat("Written:", output_pdf, "\n")
