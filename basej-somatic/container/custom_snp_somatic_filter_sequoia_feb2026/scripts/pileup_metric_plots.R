## Pileup QC metric density plots (overall distribution, no Pass/Fail split).
## Reads all res_pileup_all_group_<group>_chr*.tsv files in the current directory,
## pools all per-sample rows, and plots density distributions of key alignment
## quality metrics with meaningful x-axis scales per metric type.
## Usage: Rscript pileup_metric_plots.R <group_id> <output.pdf>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript pileup_metric_plots.R <group_id> <output.pdf>")
group_id   <- args[1]
output_pdf <- args[2]

library(ggplot2)
library(scales)
suppressPackageStartupMessages(library(data.table))

## ── Load and pool all res_pileup files ────────────────────────────────────────
pileup_files <- list.files(".", pattern = paste0("^res_pileup_all_group_", group_id, "_chr.*[.]tsv$"),
                           full.names = TRUE)
if (length(pileup_files) == 0)
  stop("No res_pileup_all_group_", group_id, "_chr*.tsv files found in current directory.")

cat("Loading", length(pileup_files), "pileup files...\n")
# fread + rbindlist: lower peak memory than do.call(rbind, lapply(read.table ...)).
# Returned as plain data.frame so downstream df[[col]] / is.finite() / etc. behave
# identically to the previous code.
pil <- as.data.frame(rbindlist(lapply(pileup_files, fread,
                                       sep = "\t", header = TRUE,
                                       data.table = FALSE,
                                       showProgress = FALSE)))
n_rows <- nrow(pil)
cat("Total sample-variant rows:", format(n_rows, big.mark = ","), "\n\n")

gp <- paste0("  |  Group: ", group_id,
             "  |  n = ", format(n_rows, big.mark = ","), " sample-variant rows")

## ── Shared theme ──────────────────────────────────────────────────────────────
theme_clean <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title       = element_text(face = "bold", size = 13),
        strip.text       = element_text(size = 11, face = "bold"),
        legend.position  = "none")

FILL_COL   <- "#2166AC"
MEDIAN_COL <- "#B2182B"

## Helper: single density with median + IQR annotations
make_density <- function(df, col, title, xlab, subtitle,
                         x_limits, x_breaks, x_labels = waiver()) {
  v   <- df[[col]][is.finite(df[[col]])]
  med <- median(v)
  q25 <- quantile(v, 0.25)
  q75 <- quantile(v, 0.75)

  ann <- paste0("Median = ", round(med, 3),
                "   IQR = [", round(q25, 3), ", ", round(q75, 3), "]")

  ggplot(data.frame(x = v), aes(x = x)) +
    geom_density(fill = FILL_COL, colour = FILL_COL, alpha = 0.4, linewidth = 0.8) +
    geom_vline(xintercept = med, colour = MEDIAN_COL,
               linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = med, y = Inf,
             label = paste0(" median\n ", round(med, 2)),
             hjust = 0, vjust = 1.3, colour = MEDIAN_COL, size = 3.5) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks,
                       labels = x_labels, oob = squish) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title    = paste0(title, gp),
         x        = xlab,
         y        = "Density",
         subtitle = paste0(subtitle, "\n", ann)) +
    theme_clean
}

## Helper: F/R pair on one page via facet
make_pair_density <- function(df, cols, strand_labels, title, xlab, subtitle,
                               x_limits, x_breaks, x_labels = waiver()) {
  long <- do.call(rbind, lapply(cols, function(col) {
    v <- df[[col]]
    data.frame(Value  = v,
               Strand = strand_labels[[col]],
               stringsAsFactors = FALSE)
  }))
  long <- long[is.finite(long$Value), ]
  long$Strand <- factor(long$Strand, levels = as.character(strand_labels))

  med_df <- aggregate(Value ~ Strand, data = long, FUN = median)

  ggplot(long, aes(x = Value)) +
    geom_density(fill = FILL_COL, colour = FILL_COL, alpha = 0.4, linewidth = 0.8) +
    geom_vline(data = med_df, aes(xintercept = Value),
               colour = MEDIAN_COL, linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ Strand, ncol = 2, scales = "free_y") +
    scale_x_continuous(limits = x_limits, breaks = x_breaks,
                       labels = x_labels, oob = squish) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title    = paste0(title, gp),
         x        = xlab,
         y        = "Density",
         subtitle = paste0(subtitle, "\nDashed red line = median. Facets: Forward vs Reverse strand.")) +
    theme_clean
}

## ─────────────────────────────────────────────────────────────────────────────
## Plot 1: Median allele score of variant-supporting reads
## AS ranges 0–150 (capped at read length × match score for 150 bp reads).
## High values = well-aligned reads; low values flag misalignment artifacts.
## ─────────────────────────────────────────────────────────────────────────────
p_as <- make_density(
  pil, "MEDIAN_AS_VARIANT_READS",
  title    = "Median allele score (AS) of variant-supporting reads",
  xlab     = "Median AS score  (0 = worst alignment, 150 = perfect match for 150 bp reads)",
  subtitle = "Low AS scores flag poorly aligned variant-supporting reads (misalignment artifacts)",
  x_limits = c(0, 150),
  x_breaks = seq(0, 150, 10)
)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 2: Proportion of clipped bases
## Proportion of bases in variant-supporting reads that are soft-clipped.
## Values near 0 = reads align cleanly; values near 1 = heavily clipped reads.
## ─────────────────────────────────────────────────────────────────────────────
p_clip <- make_density(
  pil, "PROP_BASES_CLIPPED",
  title    = "Proportion of clipped bases in variant-supporting reads",
  xlab     = "Proportion of clipped bases  (0 = no clipping, 1 = all bases clipped)",
  subtitle = "High clipping fraction indicates reads that could not be fully aligned — artifact signal",
  x_limits = c(0, 1),
  x_breaks = seq(0, 1, 0.1),
  x_labels = percent_format(accuracy = 1)
)

## ─────────────────────────────────────────────────────────────────────────────
## Plots 3–4: Proportion of fragments with bp-start in lower / upper position tail
## These capture base-position bias: if variant-supporting reads start at positions
## clustered at the extreme ends of the amplicon, it may indicate PCR or alignment
## artifacts rather than a true somatic variant.
## Most true variants have bp-starts spread uniformly → UNDER and UPPER near 0.
## ─────────────────────────────────────────────────────────────────────────────
strand_labels_prop <- c(
  "PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F" = "Forward strand",
  "PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R" = "Reverse strand"
)
p_under <- make_pair_density(
  pil,
  cols          = names(strand_labels_prop),
  strand_labels = strand_labels_prop,
  title    = "Prop. fragments with bp-start in lower position tail (BPPos UNDER)",
  xlab     = "Proportion of fragments  (0 = none start at lower tail, 1 = all start there)",
  subtitle = "High values = fragment starts cluster at the lower end of the read position range — positional bias artifact",
  x_limits = c(0, 1),
  x_breaks = seq(0, 1, 0.1),
  x_labels = percent_format(accuracy = 1)
)

strand_labels_upper <- c(
  "PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F" = "Forward strand",
  "PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R" = "Reverse strand"
)
p_upper <- make_pair_density(
  pil,
  cols          = names(strand_labels_upper),
  strand_labels = strand_labels_upper,
  title    = "Prop. fragments with bp-start in upper position tail (BPPos UPPER)",
  xlab     = "Proportion of fragments  (0 = none start at upper tail, 1 = all start there)",
  subtitle = "High values = fragment starts cluster at the upper end of the read position range — positional bias artifact",
  x_limits = c(0, 1),
  x_breaks = seq(0, 1, 0.1),
  x_labels = percent_format(accuracy = 1)
)

## ─────────────────────────────────────────────────────────────────────────────
## Plots 5–6: SD and MAD of fragment bp-start positions
## Measure the spread of read-start positions across the variant-supporting reads.
## Units = base pairs (bp). High SD/MAD = starts are evenly spread across the
## amplicon (expected for real variants). Low SD/MAD = starts cluster at a single
## position, suggesting PCR duplicates or alignment pile-up artifacts.
## Observed range: ~0–103 bp; median ~38 bp.
## ─────────────────────────────────────────────────────────────────────────────
strand_labels_sd <- c(
  "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F" = "Forward strand",
  "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R" = "Reverse strand"
)
p_sd <- make_pair_density(
  pil,
  cols          = names(strand_labels_sd),
  strand_labels = strand_labels_sd,
  title    = "Standard deviation (SD) of fragment bp-start positions",
  xlab     = "SD of bp-start positions  (base pairs)",
  subtitle = "Low SD = read starts pile up at one position (duplication/artifact); high SD = well-spread starts (expected)",
  x_limits = c(0, 110),
  x_breaks = seq(0, 100, 10)
)

strand_labels_mad <- c(
  "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F" = "Forward strand",
  "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R" = "Reverse strand"
)
p_mad <- make_pair_density(
  pil,
  cols          = names(strand_labels_mad),
  strand_labels = strand_labels_mad,
  title    = "Median absolute deviation (MAD) of fragment bp-start positions",
  xlab     = "MAD of bp-start positions  (base pairs)",
  subtitle = "Robust version of SD; less sensitive to outlier read-start positions. Low MAD = positional clustering",
  x_limits = c(0, 110),
  x_breaks = seq(0, 100, 10)
)

## ── Histogram helpers ─────────────────────────────────────────────────────────
make_histogram <- function(df, col, title, xlab, subtitle,
                           x_limits, x_breaks, binwidth,
                           x_labels = waiver()) {
  v   <- df[[col]][is.finite(df[[col]])]
  med <- median(v)
  ann <- paste0("Median = ", round(med, 3),
                "   IQR = [", round(quantile(v, 0.25), 3),
                ", ", round(quantile(v, 0.75), 3), "]")

  ggplot(data.frame(x = v), aes(x = x)) +
    geom_histogram(binwidth = binwidth, fill = FILL_COL, colour = "white",
                   linewidth = 0.2) +
    geom_vline(xintercept = med, colour = MEDIAN_COL,
               linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = med, y = Inf,
             label = paste0(" median\n ", round(med, 2)),
             hjust = 0, vjust = 1.3, colour = MEDIAN_COL, size = 3.5) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks,
                       labels = x_labels, oob = squish) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = label_comma()) +
    labs(title    = paste0(title, gp),
         x        = paste0(xlab, "  [bin width = ", binwidth, "]"),
         y        = "Number of sample-variant rows",
         subtitle = paste0(subtitle, "\n", ann)) +
    theme_clean
}

make_pair_histogram <- function(df, cols, strand_labels, title, xlab, subtitle,
                                x_limits, x_breaks, binwidth,
                                x_labels = waiver()) {
  long <- do.call(rbind, lapply(cols, function(col) {
    v <- df[[col]]
    data.frame(Value  = v,
               Strand = strand_labels[[col]],
               stringsAsFactors = FALSE)
  }))
  long   <- long[is.finite(long$Value), ]
  long$Strand <- factor(long$Strand, levels = as.character(strand_labels))
  med_df <- aggregate(Value ~ Strand, data = long, FUN = median)

  ggplot(long, aes(x = Value)) +
    geom_histogram(binwidth = binwidth, fill = FILL_COL, colour = "white",
                   linewidth = 0.2) +
    geom_vline(data = med_df, aes(xintercept = Value),
               colour = MEDIAN_COL, linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ Strand, ncol = 2, scales = "free_y") +
    scale_x_continuous(limits = x_limits, breaks = x_breaks,
                       labels = x_labels, oob = squish) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = label_comma()) +
    labs(title    = paste0(title, gp),
         x        = paste0(xlab, "  [bin width = ", binwidth, "]"),
         y        = "Number of sample-variant rows",
         subtitle = paste0(subtitle,
                           "\nDashed red line = median. Facets: Forward vs Reverse strand.")) +
    theme_clean
}

## ── Histogram versions ────────────────────────────────────────────────────────
h_as <- make_histogram(
  pil, "MEDIAN_AS_VARIANT_READS",
  title    = "Median allele score (AS) of variant-supporting reads",
  xlab     = "Median AS score  (0 = worst, 150 = perfect match for 150 bp reads)",
  subtitle = "Low AS scores flag poorly aligned variant-supporting reads (misalignment artifacts)",
  x_limits = c(0, 150), x_breaks = seq(0, 150, 10), binwidth = 1
)

h_clip <- make_histogram(
  pil, "PROP_BASES_CLIPPED",
  title    = "Proportion of clipped bases in variant-supporting reads",
  xlab     = "Proportion of clipped bases  (0 = no clipping, 1 = all bases clipped)",
  subtitle = "High clipping fraction indicates reads that could not be fully aligned — artifact signal",
  x_limits = c(0, 1), x_breaks = seq(0, 1, 0.1),
  x_labels = percent_format(accuracy = 1), binwidth = 0.02
)

h_under <- make_pair_histogram(
  pil, cols = names(strand_labels_prop), strand_labels = strand_labels_prop,
  title    = "Prop. fragments with bp-start in lower position tail (BPPos UNDER)",
  xlab     = "Proportion of fragments  (0 = none start at lower tail, 1 = all start there)",
  subtitle = "High values = fragment starts cluster at the lower end of the read position range",
  x_limits = c(0, 1), x_breaks = seq(0, 1, 0.1),
  x_labels = percent_format(accuracy = 1), binwidth = 0.02
)

h_upper <- make_pair_histogram(
  pil, cols = names(strand_labels_upper), strand_labels = strand_labels_upper,
  title    = "Prop. fragments with bp-start in upper position tail (BPPos UPPER)",
  xlab     = "Proportion of fragments  (0 = none start at upper tail, 1 = all start there)",
  subtitle = "High values = fragment starts cluster at the upper end of the read position range",
  x_limits = c(0, 1), x_breaks = seq(0, 1, 0.1),
  x_labels = percent_format(accuracy = 1), binwidth = 0.02
)

h_sd <- make_pair_histogram(
  pil, cols = names(strand_labels_sd), strand_labels = strand_labels_sd,
  title    = "Standard deviation (SD) of fragment bp-start positions",
  xlab     = "SD of bp-start positions  (base pairs)",
  subtitle = "Low SD = read starts pile up at one position (duplication/artifact); high SD = well-spread starts",
  x_limits = c(0, 110), x_breaks = seq(0, 100, 10), binwidth = 2
)

h_mad <- make_pair_histogram(
  pil, cols = names(strand_labels_mad), strand_labels = strand_labels_mad,
  title    = "Median absolute deviation (MAD) of fragment bp-start positions",
  xlab     = "MAD of bp-start positions  (base pairs)",
  subtitle = "Robust version of SD; low MAD = positional clustering of read starts",
  x_limits = c(0, 110), x_breaks = seq(0, 100, 10), binwidth = 2
)

## ── Write PDF ─────────────────────────────────────────────────────────────────
pdf(output_pdf, width = 12, height = 7)
print(p_as);    print(h_as)
print(p_clip);  print(h_clip)
print(p_under); print(h_under)
print(p_upper); print(h_upper)
print(p_sd);    print(h_sd)
print(p_mad);   print(h_mad)
dev.off()
cat("Written:", output_pdf, "(12 pages)\n")
