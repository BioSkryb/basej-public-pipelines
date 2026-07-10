#!/usr/bin/env Rscript
# explore_qc_filter_thresholds.R
# ---------------------------------------------------------------------------
# Threshold-calibration view for the artifact QC filter of
# CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES (AS / PropClipped / BPPos).
#
# For each filter it plots the underlying per-cell metric distribution, draws
# the CURRENT configured cutoff(s) on top, and (for the univariate AS /
# PropClipped gates) a pass-rate-vs-threshold sensitivity curve so the "knee"
# is visible. The universe is restricted to variant-supporting cell-variants
# (ALT rows with NV>0); a Scheme_HQStat flag (from scheme_membership) splits
# the post-2nd-statistical-filter subset out from all evaluable cell-variants,
# since the QC filter's practical job is to clean up within the HQStat set.
#
# Memory-safe: per-chr res_pileup files are read ONE AT A TIME and reduced to
# the variant-supporting subset before accumulation (the orphaned
# pileup_metric_plots.R / pileup_bppos_plots.R rbind the full ~140M-row pool
# at once and OOM).
#
# Usage:
#   Rscript explore_qc_filter_thresholds.R <group> <scheme_membership.tsv> \
#       <cutoff_as> <cutoff_prop_clipped> <cutoff_prop_bp_under> <cutoff_prop_bp_upper> \
#       <cutoff_sd_indiv> <cutoff_mad_indiv> <cutoff_sd_both> <cutoff_mad_both> \
#       <cutoff_sd_extreme> <cutoff_mad_extreme>
# Reads res_pileup_all_group_<group>_chr*.tsv from the current directory.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(egg)
})
if (file.exists("/usr/local/bin/theme_ohchibi_pubr.R")) {
  source("/usr/local/bin/theme_ohchibi_pubr.R")
} else { theme_ohchibi_pubr <- function(...) ggplot2::theme_bw() }

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 12) stop("Usage: explore_qc_filter_thresholds.R <group> <scheme_membership> <cutoff_as> <cutoff_prop_clipped> <cutoff_prop_bp_under> <cutoff_prop_bp_upper> <cutoff_sd_indiv> <cutoff_mad_indiv> <cutoff_sd_both> <cutoff_mad_both> <cutoff_sd_extreme> <cutoff_mad_extreme>")
group     <- args[1]
memb_f    <- args[2]
cut_as              <- as.numeric(args[3])
cut_prop            <- as.numeric(args[4])
cut_prop_bp_under   <- as.numeric(args[5])
cut_prop_bp_upper   <- as.numeric(args[6])
cut_sd_indiv        <- as.numeric(args[7])
cut_mad_indiv       <- as.numeric(args[8])
cut_sd_both         <- as.numeric(args[9])
cut_mad_both        <- as.numeric(args[10])
cut_sd_extreme      <- as.numeric(args[11])
cut_mad_extreme     <- as.numeric(args[12])

# ── HQStat variant set (post 2nd statistical filter) ─────────────────────────
hq_ids <- character(0)
if (file.exists(memb_f) && file.info(memb_f)$size > 0) {
  memb <- fread(memb_f, sep = "\t", header = TRUE, data.table = FALSE)
  if ("Scheme_HQStat" %in% colnames(memb))
    hq_ids <- as.character(memb$VariantId[memb$Scheme_HQStat == 1])
}
cat(sprintf("[explore_qc] HQStat variants: %d\n", length(hq_ids)))

# ── Stream per-chr pileup files, keep only variant-supporting cell-variants ──
need <- c("VariantId","ALT","NUM_FRAGMENTS_HQ_MQ_BQ_F","NUM_FRAGMENTS_HQ_MQ_BQ_R",
          "MEDIAN_AS_VARIANT_READS","PROP_BASES_CLIPPED",
          "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F","SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
          "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F","MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
          "PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F","PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R",
          "PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F","PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R",
          "AS_Filter","PropClipped_Filter","BPPos_Filter","Verdict")

files <- list.files(".", pattern = paste0("^res_pileup_all_group_", group, "_chr.*[.]tsv$"))
if (length(files) == 0) stop("No res_pileup_all_group_", group, "_chr*.tsv files found")
cat(sprintf("[explore_qc] streaming %d pileup file(s)\n", length(files)))

# Column projection: parse only the needed fields (computed once from the shared header),
# so each per-chr file is read cheaply. Same rows/values as before.
sel <- intersect(need, names(fread(files[1], sep = "\t", header = TRUE, nrows = 0)))
acc <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- fread(files[i], sep = "\t", header = TRUE, data.table = FALSE,
             showProgress = FALSE, select = sel)
  # variant-supporting cells only: ALT rows with NV>0 (the metrics describe variant reads)
  nv <- d$NUM_FRAGMENTS_HQ_MQ_BQ_F + d$NUM_FRAGMENTS_HQ_MQ_BQ_R
  d <- d[d$ALT != "REF" & !is.na(nv) & nv > 0, , drop = FALSE]
  acc[[i]] <- d
}
pil <- rbindlist(acc, fill = TRUE); rm(acc); gc()
pil <- as.data.frame(pil)
pil$HQStat <- pil$VariantId %in% hq_ids
cat(sprintf("[explore_qc] variant-supporting cell-variants: %d (HQStat: %d)\n",
            nrow(pil), sum(pil$HQStat)))

# Primary universe for calibration = HQStat cell-variants; fall back to all if empty.
uni <- pil[pil$HQStat, , drop = FALSE]
uni_label <- "HQStat cell-variants"
if (nrow(uni) == 0) { uni <- pil; uni_label <- "all evaluable cell-variants (no HQStat)" }

sub <- paste0(group, "  -  universe: ", uni_label, "  (n=", format(nrow(uni), big.mark=","), ")")

# ── helpers ──────────────────────────────────────────────────────────────────
p0 <- function(x) x[is.finite(x)]
pct <- function(x) sprintf("%.1f%%", 100*x)

# distribution panel with cutoff line(s)
dist_panel <- function(values, cutoffs, cut_labs, title, xlab, side = "right",
                       xcap = NULL, fill = "#4393C3") {
  v <- p0(values)
  if (!is.null(xcap)) v <- pmin(v, xcap)
  df <- data.frame(v = v)
  g <- ggplot(df, aes(v)) +
    geom_histogram(bins = 60, fill = fill, colour = NA, alpha = 0.85) +
    geom_vline(xintercept = cutoffs, colour = "#B2182B", linetype = "dashed", linewidth = 0.7) +
    labs(title = title, x = xlab, y = "cell-variants") +
    theme_ohchibi_pubr()
  ymax <- max(ggplot_build(g)$data[[1]]$count, na.rm = TRUE)
  lab_df <- data.frame(x = cutoffs, lab = cut_labs)
  g + geom_text(data = lab_df, aes(x = x, y = ymax, label = lab),
                inherit.aes = FALSE, angle = 90, vjust = -0.3, hjust = 1,
                colour = "#B2182B", size = 3)
}

# sensitivity: fraction of cell-variants passing as the threshold sweeps
sens_panel <- function(values, cutoff, dir = ">=", title, xlab) {
  v <- p0(values)
  grid <- seq(min(v), quantile(v, 0.995), length.out = 200)
  frac <- if (dir == ">=") sapply(grid, function(t) mean(v >= t)) else sapply(grid, function(t) mean(v < t))
  df <- data.frame(t = grid, frac = frac)
  cur <- if (dir == ">=") mean(v >= cutoff) else mean(v < cutoff)
  ggplot(df, aes(t, frac)) +
    geom_line(colour = "#2166AC", linewidth = 0.9) +
    geom_vline(xintercept = cutoff, colour = "#B2182B", linetype = "dashed") +
    annotate("point", x = cutoff, y = cur, colour = "#B2182B", size = 2.5) +
    annotate("text", x = cutoff, y = cur, label = paste0("  current: ", pct(cur)),
             hjust = 0, colour = "#B2182B", size = 3.4) +
    scale_y_continuous(labels = percent, limits = c(0,1)) +
    labs(title = title, subtitle = paste0("pass rule: metric ", dir, " threshold"),
         x = paste0(xlab, " threshold"), y = "cell-variants passing") +
    theme_ohchibi_pubr()
}

# ── AS & PropClipped ─────────────────────────────────────────────────────────
as_v   <- uni$MEDIAN_AS_VARIANT_READS
prop_v <- uni$PROP_BASES_CLIPPED
as_dist  <- dist_panel(as_v, cut_as, paste0("cutoff_as=", cut_as),
                       "AS  (MEDIAN_AS_VARIANT_READS)", "median alignment score", fill = "#4393C3")
as_sens  <- sens_panel(as_v, cut_as, ">=", "AS pass-rate sensitivity", "AS")
prop_dist<- dist_panel(prop_v, cut_prop, paste0("cutoff=", cut_prop),
                       "PropClipped  (PROP_BASES_CLIPPED)", "prop. bases soft-clipped", fill = "#92C5DE")
prop_sens<- sens_panel(prop_v, cut_prop, "<",  "PropClipped pass-rate sensitivity", "PropClipped")

# ── BPPos metric distributions (pool forward + reverse strands) ──────────────
sd_v   <- c(uni$SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F,  uni$SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R)
mad_v  <- c(uni$MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F, uni$MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R)
pu_v   <- c(uni$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F, uni$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R)
pp_v   <- c(uni$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F, uni$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R)
sd_dist  <- dist_panel(sd_v,  c(cut_sd_indiv, cut_sd_both, cut_sd_extreme),
                       c(paste0("indiv=",cut_sd_indiv), paste0("both=",cut_sd_both), paste0("extreme=",cut_sd_extreme)),
                       "BPPos  SD of fragment bp-start", "SD(bp-start)", fill = "#F4A582")
mad_dist <- dist_panel(mad_v, c(cut_mad_indiv, cut_mad_both, cut_mad_extreme),
                       c(paste0("indiv=",cut_mad_indiv), paste0("both=",cut_mad_both), paste0("extreme=",cut_mad_extreme)),
                       "BPPos  MAD of fragment bp-start", "MAD(bp-start)", fill = "#F4A582")
pu_dist  <- dist_panel(pu_v, cut_prop_bp_under, paste0("under=",cut_prop_bp_under),
                       "BPPos  prop. fragments bp-start UNDER", "prop. under", fill = "#FDDBC7")
pp_dist  <- dist_panel(pp_v, cut_prop_bp_upper, paste0("upper=",cut_prop_bp_upper),
                       "BPPos  prop. fragments bp-start UPPER", "prop. upper", fill = "#FDDBC7")

# ── per-cell filter pass rates (bar) ─────────────────────────────────────────
pass_rate <- function(x) mean(x == "Pass", na.rm = TRUE)
# True per-cell artifact-QC pass = AS & PropClipped & BPPos all Pass in that cell.
# (The res_pileup `Verdict` column is variant-level matrix membership broadcast to every
# cell, NOT the per-cell AND, so it is deliberately not used here.)
all3 <- (uni$AS_Filter == "Pass") & (uni$PropClipped_Filter == "Pass") & (uni$BPPos_Filter == "Pass")
bar_df <- data.frame(
  filter = factor(c("AS","PropClipped","BPPos","All 3 (AND)"),
                  levels = c("AS","PropClipped","BPPos","All 3 (AND)")),
  frac = c(pass_rate(uni$AS_Filter), pass_rate(uni$PropClipped_Filter),
           pass_rate(uni$BPPos_Filter), mean(all3, na.rm = TRUE)))
bar_panel <- ggplot(bar_df, aes(filter, frac)) +
  geom_col(fill = "#2166AC", width = 0.65) +
  geom_text(aes(label = pct(frac)), vjust = -0.4, fontface = "bold", size = 4) +
  scale_y_continuous(labels = percent, limits = c(0, 1.05), expand = expansion(mult = c(0,0.05))) +
  labs(title = "Per-cell filter pass rate", subtitle = sub,
       x = NULL, y = "cell-variants passing") +
  theme_ohchibi_pubr()

# ── render: multipage PDF + per-page PNG (MultiQC-ready) ─────────────────────
p_as   <- egg::ggarrange(as_dist,  as_sens,   prop_dist, prop_sens, ncol = 2, draw = FALSE)
p_bp   <- egg::ggarrange(sd_dist,  mad_dist,  pu_dist,   pp_dist,   ncol = 2, draw = FALSE)

pdf(sprintf("qc_threshold_distributions_%s.pdf", group), width = 12, height = 9)
print(p_as); print(p_bp); print(bar_panel)
invisible(dev.off())

ggsave(sprintf("qc_threshold_as_propclipped_%s.png", group), p_as, width = 12, height = 9, dpi = 150, bg = "white")
ggsave(sprintf("qc_threshold_bppos_%s.png", group),          p_bp, width = 12, height = 9, dpi = 150, bg = "white")
ggsave(sprintf("qc_threshold_passrates_%s.png", group),  bar_panel, width = 9,  height = 5, dpi = 150, bg = "white")

# ── summary TSV: current cutoff + pass rate (univariate gates) ───────────────
summ <- data.frame(
  filter    = c("AS","PropClipped","BPPos"),
  metric    = c("MEDIAN_AS_VARIANT_READS","PROP_BASES_CLIPPED","(multivariate: prop+SD+MAD x strand-case)"),
  rule      = c(paste0(">= ", cut_as), paste0("< ", cut_prop), "see rscript_2 cascade"),
  pass_rate = c(pct(mean(p0(as_v) >= cut_as)), pct(mean(p0(prop_v) < cut_prop)), pct(pass_rate(uni$BPPos_Filter))),
  universe  = uni_label,
  n         = nrow(uni),
  stringsAsFactors = FALSE)
write.table(summ, sprintf("qc_threshold_summary_%s.tsv", group),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[explore_qc] done\n"); print(summ)
