#!/usr/bin/env Rscript
# plot_comparative_signature_activities.R
#
# Faceted 3-panel stacked bar chart comparing COSMIC SBS signature activities
# across three methods run in the mutsignatures_from_annotated_vcf workflow:
#   SigProfilerAssignment (COSMIC v3.5), MuSiCaL (COSMIC v3.2), SigDyn (COSMIC v3.2).
#
# Sample order on the x-axis is derived from the SigProfilerAssignment activity
# matrix using Bray-Curtis distance + ward.D hierarchical clustering.
# A single shared palette is applied across all three panels; a single legend
# is displayed at the bottom.
#
# Usage:
#   Rscript plot_comparative_signature_activities.R \
#       <spa_activities>     merged_signature_activities_SBS.txt
#       <musical_activities> *_musical_activities.tsv
#       <sigdyn_activities>  *_sigdyn_activities.tsv
#       <cohort_id>
#       <output_prefix>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
})

# ── Inline theme (equivalent to theme_ohchibi_pubr.R) ────────────────────────
theme_ohchibi_pubr <- function(base_size = 13, base_family = "") {
  half_line <- base_size / 2
  base <- theme_bw(base_size = base_size, base_family = base_family)
  overlay <- theme(
    panel.background    = element_rect(fill = "white", colour = NA),
    panel.grid.major.x  = element_blank(),
    panel.grid.major.y  = element_line(linetype = "dotted", colour = "grey"),
    panel.grid.minor.y  = element_line(colour = "grey93", linewidth = 0.2),
    panel.grid.minor.x  = element_blank(),
    panel.border        = element_rect(fill = NA, colour = "black", linewidth = 0.3),
    axis.line           = element_blank(),
    axis.ticks          = element_line(colour = "black", linewidth = 0.5),
    axis.text.x         = element_text(colour = "grey30", face = "plain", size = 9,
                                        angle = 40, hjust = 1, vjust = 1),
    axis.text.y         = element_text(colour = "grey30", face = "plain", size = 9),
    axis.title.x        = element_text(colour = "black", face = "plain", size = 9),
    axis.title.y        = element_text(colour = "black", face = "plain", size = 9,
                                        angle = 90, vjust = 1),
    strip.background    = element_blank(),
    strip.text          = element_text(face = "plain", size = 9,
                                        margin = margin(b = 4)),
    legend.background   = element_blank(),
    legend.key          = element_blank(),
    legend.position     = "bottom",
    legend.direction    = "horizontal",
    legend.title        = element_text(face = "plain", colour = "black", size = 9),
    legend.text         = element_text(colour = "black", face = "plain", size = 9),
    legend.key.size     = grid::unit(0.70, "cm"),
    legend.key.width    = grid::unit(0.70, "cm"),
    legend.spacing.x    = grid::unit(0.4, "cm"),
    legend.margin       = margin(t = 8, b = 4),
    plot.title          = element_text(hjust = 0.5, face = "plain", size = 9),
    plot.subtitle       = element_text(hjust = 0.5, face = "plain", size = 9),
    plot.margin         = margin(half_line, half_line, half_line, half_line)
  )
  `%+replace%`(base, overlay)
}

# ── Helper functions (mirror of plot_matrix_signature_bargraphs.R) ────────────
bray_curtis_dist <- function(mat) {
  n <- nrow(mat)
  d <- matrix(0, n, n, dimnames = list(rownames(mat), rownames(mat)))
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      xi <- mat[i, ]; xj <- mat[j, ]
      s  <- sum(xi) + sum(xj)
      bc <- if (s > 0) 1 - 2 * sum(pmin(xi, xj)) / s else 0
      d[i, j] <- d[j, i] <- bc
    }
  }
  as.dist(d)
}

cluster_order <- function(mat) {
  if (nrow(mat) <= 1L) return(rownames(mat))
  hc <- hclust(bray_curtis_dist(mat), method = "ward.D")
  hc$labels[hc$order]
}

make_sig_colors <- function(sigs) {
  n   <- length(sigs)
  # The first four palettes are kept first and in the same order as before, so the
  # leading colours (assigned to the SPA-detected signatures) are unchanged. Paired
  # + Set2 are appended to supply extra distinct colours for signatures only seen in
  # MuSiCaL/SigDyn; if even more are needed, interpolate to length n.
  pal <- unique(c(
    brewer.pal(12, "Set3"),
    brewer.pal(8,  "Dark2"),
    brewer.pal(9,  "Set1"),
    brewer.pal(8,  "Accent"),
    brewer.pal(12, "Paired"),
    brewer.pal(8,  "Set2")
  ))
  if (n > length(pal)) pal <- grDevices::colorRampPalette(pal)(n)
  setNames(pal[seq_len(n)], sigs)
}

# ── Load activity file → long data frame ─────────────────────────────────────
load_activities <- function(path, method) {
  df  <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(df[, -1L, drop = FALSE])
  rownames(mat) <- df[[1L]]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  as.data.frame(mat) %>%
    tibble::rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Signature", values_to = "Activity") %>%
    filter(Activity > 0) %>%
    mutate(Method = method)
}

# ── Parse arguments ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop(paste(
    "Usage: plot_comparative_signature_activities.R",
    "<spa_activities> <musical_activities> <sigdyn_activities>",
    "<cohort_id> <output_prefix>"
  ))
}
spa_file     <- args[1L]
musical_file <- args[2L]
sigdyn_file  <- args[3L]
cohort_id    <- args[4L]
out_prefix   <- args[5L]

# ── Compute sample order from SPA (Bray-Curtis + ward.D) ─────────────────────
cat("Computing sample order from SigProfilerAssignment activities...\n")
spa_df  <- read.delim(spa_file, check.names = FALSE, stringsAsFactors = FALSE)
spa_mat <- as.matrix(spa_df[, -1L, drop = FALSE])
rownames(spa_mat) <- spa_df[[1L]]
spa_mat <- spa_mat[, colSums(spa_mat) > 0, drop = FALSE]
sample_order <- cluster_order(spa_mat)
cat(sprintf("  %d samples, order derived from Bray-Curtis clustering.\n", length(sample_order)))

# ── Build signature palette ordered by SPA activity ──────────────────────────
all_sigs   <- names(sort(colSums(spa_mat), decreasing = TRUE))
sig_colors <- make_sig_colors(all_sigs)

# ── Load all three activity matrices ─────────────────────────────────────────
cat("Loading activity matrices...\n")
method_order  <- c("SigProfilerAssignment", "MuSiCaL", "SigDyn")
method_labels <- c(
  SigProfilerAssignment = "SigProfilerAssignment  (COSMIC v3.5)",
  MuSiCaL               = "MuSiCaL  (COSMIC v3.2)",
  SigDyn                = "SigDyn  (COSMIC v3.2)"
)

all_long <- bind_rows(
  load_activities(spa_file,     "SigProfilerAssignment"),
  load_activities(musical_file, "MuSiCaL"),
  load_activities(sigdyn_file,  "SigDyn")
)

# Intersect sample set to samples present in SPA
common_samples <- intersect(sample_order, unique(all_long$Sample))
sample_order   <- sample_order[sample_order %in% common_samples]

# Factor levels & palette.
# SPA-detected signatures keep their COSMIC v3.5 colours and SPA-prevalence order
# (they come first, so make_sig_colors assigns them the same leading colours as
# before). Signatures reported only by MuSiCaL/SigDyn — not detected by SPA — are
# appended after, each receiving its OWN additional colour from the extended
# palette, instead of falling through to ggplot's default grey "missing" fill.
spa_present <- all_sigs[all_sigs %in% unique(all_long$Signature)]
extra_sigs  <- sort(setdiff(unique(all_long$Signature), all_sigs))
plot_sigs   <- c(spa_present, extra_sigs)
sig_colors  <- make_sig_colors(plot_sigs)

all_long <- all_long %>%
  filter(Sample %in% sample_order) %>%
  mutate(
    Sample    = factor(Sample,    levels = sample_order),
    Signature = factor(Signature, levels = plot_sigs),
    Method    = factor(Method,    levels = method_order)
  )

# ── Build comparative faceted plot ───────────────────────────────────────────
cat("Building faceted comparative plot...\n")

n_samples <- length(sample_order)
cols_here <- sig_colors[names(sig_colors) %in% plot_sigs]

p <- ggplot(all_long, aes(x = Sample, y = Activity, fill = Signature)) +
  geom_bar(stat = "identity", width = 1,
           position = position_stack(reverse = TRUE)) +
  scale_fill_manual(
    values = cols_here,
    name   = "Signature",
    breaks = names(cols_here),
    drop   = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  facet_grid(
    Method ~ .,
    scales   = "free_y",
    switch   = "y",
    labeller = labeller(Method = method_labels)
  ) +
  labs(
    title    = paste0("Comparative Mutational Signature Activities — ", cohort_id),
    subtitle = paste0(n_samples, " samples  |  x-axis ordered by SigProfilerAssignment Bray-Curtis clustering"),
    x        = NULL,
    y        = "Activity"
  ) +
  theme_ohchibi_pubr(base_size = 11) +
  theme(
    axis.text.x         = element_blank(),
    axis.ticks.x        = element_blank(),
    axis.ticks.length.x = grid::unit(0, "pt"),
    strip.placement     = "outside",
    strip.text.y.left   = element_text(angle = 0, hjust = 1, size = 9, face = "bold"),
    legend.position     = "bottom",
    legend.direction    = "horizontal"
  ) +
  guides(fill = guide_legend(
    title          = "Signature",
    nrow           = ceiling(length(cols_here) / 10),
    override.aes   = list(size = 4)
  ))

# ── Save outputs ──────────────────────────────────────────────────────────────
w_in <- min(max(20, 0.18 * n_samples + 12), 40)
h_in <- 18

png_path <- paste0(out_prefix, "_comparative_signature_activities.png")
ggsave(png_path, p, width = w_in, height = h_in, dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("Saved PNG  -> %s\n", png_path))

pdf_path <- paste0(out_prefix, "_comparative_signature_activities.pdf")
grDevices::pdf(pdf_path, width = w_in, height = h_in, useDingbats = FALSE, paper = "special")
tryCatch(print(p), finally = grDevices::dev.off())
if (!file.exists(pdf_path) || file.info(pdf_path)$size < 1L) {
  stop("Failed to write ", pdf_path, call. = FALSE)
}
cat(sprintf("Saved PDF  -> %s\n", pdf_path))

cat("Comparative signature plot complete.\n")

# ══════════════════════════════════════════════════════════════════════════════
# ── Concordance analysis ───────────────────────────────────────────────────────
#
# Three complementary views of cross-method robustness:
#   1. Signature detection rate heatmap  — which signatures each method detects
#      and how consistently (% samples with activity > 0)
#   2. Sample-level concordance          — per-sample cosine similarity between
#      method pairs + total-activity scatter with Pearson r and Spearman ρ
#   3. Signature-level Spearman ρ        — for each detected signature, how
#      correlated are per-sample activities between method pairs
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(patchwork))

cat("Running concordance analyses...\n")

# ── Rebuild full per-method matrices (all values, zeros filled for missing sigs) ─
rebuild_wide <- function(path) {
  df  <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(df[, -1L, drop = FALSE])
  rownames(mat) <- df[[1L]]
  mat[is.na(mat)] <- 0
  mat
}
spa_mat_raw    <- rebuild_wide(spa_file)
music_mat_raw  <- rebuild_wide(musical_file)
sigdyn_mat_raw <- rebuild_wide(sigdyn_file)

# Align to common samples × union of all signatures
common_samps <- Reduce(intersect, list(
  rownames(spa_mat_raw), rownames(music_mat_raw), rownames(sigdyn_mat_raw)
))
sigs_union <- Reduce(union, list(
  colnames(spa_mat_raw), colnames(music_mat_raw), colnames(sigdyn_mat_raw)
))

expand_mat <- function(mat, samps, sigs) {
  out <- matrix(0.0, nrow = length(samps), ncol = length(sigs),
                dimnames = list(samps, sigs))
  rs <- intersect(rownames(mat), samps)
  cs <- intersect(colnames(mat), sigs)
  if (length(rs) > 0L && length(cs) > 0L) out[rs, cs] <- mat[rs, cs]
  out
}
spa_m    <- expand_mat(spa_mat_raw,    common_samps, sigs_union)
music_m  <- expand_mat(music_mat_raw,  common_samps, sigs_union)
sigdyn_m <- expand_mat(sigdyn_mat_raw, common_samps, sigs_union)

n_com <- length(common_samps)
cat(sprintf("  %d common samples, %d union signatures.\n", n_com, length(sigs_union)))

# ── Big-picture summary ────────────────────────────────────────────────────────
det_df <- data.frame(
  Signature = sigs_union,
  SPA       = colMeans(spa_m    > 0) * 100,
  MuSiCaL   = colMeans(music_m  > 0) * 100,
  SigDyn    = colMeans(sigdyn_m > 0) * 100,
  stringsAsFactors = FALSE
)
det_df$max_det  <- apply(det_df[, 2:4], 1, max)
det_df$n_detect <- rowSums(det_df[, 2:4] > 0)

cat(sprintf(
  "  Signature overlap:\n    all 3 methods : %d\n    any 2 methods : %d\n    SPA only      : %d  |  MuSiCaL only: %d  |  SigDyn only: %d\n",
  sum(det_df$n_detect == 3L),
  sum(det_df$n_detect >= 2L),
  sum(det_df$SPA > 0 & det_df$MuSiCaL == 0 & det_df$SigDyn == 0),
  sum(det_df$SPA == 0 & det_df$MuSiCaL > 0 & det_df$SigDyn == 0),
  sum(det_df$SPA == 0 & det_df$MuSiCaL == 0 & det_df$SigDyn > 0)
))

# ── Plot 1: Signature detection rate heatmap ───────────────────────────────────
det_plot_df <- det_df[det_df$max_det >= 5, ]
det_plot_df <- det_plot_df[order(-det_plot_df$max_det, -rowMeans(det_plot_df[, 2:4])), ]

det_long <- det_plot_df %>%
  select(Signature, SPA, MuSiCaL, SigDyn) %>%
  pivot_longer(-Signature, names_to = "Method", values_to = "DetRate") %>%
  mutate(
    Method    = factor(Method, levels = c("SPA", "MuSiCaL", "SigDyn")),
    Signature = factor(Signature, levels = rev(det_plot_df$Signature)),
    lbl       = ifelse(DetRate >= 3, sprintf("%.0f%%", DetRate), ""),
    txt_col   = ifelse(DetRate > 60, "white", "grey20")
  )

p_det <- ggplot(det_long, aes(x = Method, y = Signature, fill = DetRate)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = lbl, colour = txt_col), size = 2.5, show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colors = c("white", "#deebf7", "#9ecae1", "#3182bd", "#08306b"),
    limits = c(0, 100), name = "% samples\ndetected"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title    = paste0("Signature Detection Rate — ", cohort_id),
    subtitle = sprintf(
      "%d samples  |  %d sigs detected by all 3  |  %d by ≥2 methods",
      n_com, sum(det_df$n_detect == 3L), sum(det_df$n_detect >= 2L)
    ),
    x = NULL, y = NULL
  ) +
  theme_ohchibi_pubr(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 9, face = "bold"),
    axis.text.y        = element_text(size = 8),
    legend.position    = "right",
    legend.direction   = "vertical"
  )

det_h <- min(max(5.0, nrow(det_plot_df) * 0.38 + 2.5), 40)
ggsave(paste0(out_prefix, "_concordance_detection.png"), p_det,
       width = 6, height = det_h, dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("Saved detection heatmap  -> %s_concordance_detection.png\n", out_prefix))

# ── Plot 2: Sample-level concordance ──────────────────────────────────────────
# 2a. Per-sample cosine similarity between method pairs
.cos_sim <- function(a, b) {
  d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (d == 0) NA_real_ else sum(a * b) / d
}

cs_df <- data.frame(
  Sample               = common_samps,
  `SPA vs MuSiCaL`    = vapply(seq_len(n_com), function(i) .cos_sim(spa_m[i, ], music_m[i, ]),  numeric(1L)),
  `SPA vs SigDyn`     = vapply(seq_len(n_com), function(i) .cos_sim(spa_m[i, ], sigdyn_m[i, ]), numeric(1L)),
  `MuSiCaL vs SigDyn` = vapply(seq_len(n_com), function(i) .cos_sim(music_m[i, ], sigdyn_m[i, ]), numeric(1L)),
  check.names = FALSE
)

cs_long <- cs_df %>%
  pivot_longer(-Sample, names_to = "Pair", values_to = "CS") %>%
  mutate(Pair = factor(Pair, levels = c("SPA vs MuSiCaL", "SPA vs SigDyn", "MuSiCaL vs SigDyn")))

cs_stats <- cs_long %>%
  group_by(Pair) %>%
  summarise(med = median(CS, na.rm = TRUE), mn = mean(CS, na.rm = TRUE), .groups = "drop")

pair_pal <- c("SPA vs MuSiCaL" = "#4393c3", "SPA vs SigDyn" = "#d6604d", "MuSiCaL vs SigDyn" = "#74c476")

p_cs <- ggplot(cs_long, aes(x = Pair, y = CS, fill = Pair)) +
  geom_violin(alpha = 0.6, width = 0.85, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.size = 0.4, fill = "white", alpha = 0.85,
               outlier.alpha = 0.25) +
  geom_text(data = cs_stats,
            aes(x = Pair, y = 1.12, label = sprintf("med %.2f / mean %.2f", med, mn)),
            size = 2.8, hjust = 0.5, lineheight = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(values = pair_pal, guide = "none") +
  scale_y_continuous(limits = c(NA, 1.22), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Per-sample Cosine Similarity",
    subtitle = "Signature profile similarity between method pairs",
    x = NULL, y = "Cosine Similarity"
  ) +
  theme_ohchibi_pubr(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9))

# 2b. Total activity scatter plots (one per method pair)
tot_df <- data.frame(
  SPA     = rowSums(spa_m),
  MuSiCaL = rowSums(music_m),
  SigDyn  = rowSums(sigdyn_m)
)

.scatter <- function(x_col, y_col, x_lab, y_lab) {
  x  <- tot_df[[x_col]]; y <- tot_df[[y_col]]
  pr <- cor(x, y, method = "pearson",  use = "complete.obs")
  sr <- cor(x, y, method = "spearman", use = "complete.obs")
  ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_bin2d(bins = 45) +
    scale_fill_gradientn(
      colors = c("#f7fbff", "#9ecae1", "#2171b5", "#08306b"), name = "n"
    ) +
    geom_smooth(method = "lm", se = FALSE, colour = "#d73027", linewidth = 0.7,
                formula = y ~ x) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.3,
             label = sprintf("r=%.2f  ρ=%.2f", pr, sr), size = 2.8) +
    labs(x = paste0(x_lab, "\ntotal activity"),
         y = paste0(y_lab, "\ntotal activity")) +
    theme_ohchibi_pubr(base_size = 9) +
    theme(legend.position = "none")
}

p_sc1 <- .scatter("SPA",     "MuSiCaL", "SPA",     "MuSiCaL")
p_sc2 <- .scatter("SPA",     "SigDyn",  "SPA",     "SigDyn")
p_sc3 <- .scatter("MuSiCaL", "SigDyn",  "MuSiCaL", "SigDyn")

p_sample <- (p_sc1 | p_sc2 | p_sc3) / p_cs +
  plot_annotation(
    title    = paste0("Sample-level Method Concordance — ", cohort_id),
    subtitle = paste0(n_com, " samples"),
    theme    = theme(
      plot.title    = element_text(hjust = 0.5, size = 11, face = "plain"),
      plot.subtitle = element_text(hjust = 0.5, size = 9,  face = "plain")
    )
  ) +
  plot_layout(heights = c(1, 1.3))

ggsave(paste0(out_prefix, "_concordance_sample.png"), p_sample,
       width = 14, height = 10, dpi = 300, bg = "white")
cat(sprintf("Saved sample concordance -> %s_concordance_sample.png\n", out_prefix))

# ── Plot 3: Per-signature Spearman correlation heatmap ─────────────────────────
sig_corr_df <- bind_rows(lapply(sigs_union, function(sig) {
  a <- spa_m[, sig]; b <- music_m[, sig]; cc <- sigdyn_m[, sig]
  if (max(mean(a > 0), mean(b > 0), mean(cc > 0)) * 100 < 5) return(NULL)
  data.frame(
    Signature           = sig,
    `SPA vs MuSiCaL`   = suppressWarnings(cor(a, b,  method = "spearman")),
    `SPA vs SigDyn`    = suppressWarnings(cor(a, cc, method = "spearman")),
    `MuSiCaL vs SigDyn`= suppressWarnings(cor(b, cc, method = "spearman")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}))
sig_corr_df$mean_r <- rowMeans(sig_corr_df[, 2:4], na.rm = TRUE)
sig_corr_df <- sig_corr_df[order(-sig_corr_df$mean_r), ]

sig_corr_long <- sig_corr_df %>%
  select(-mean_r) %>%
  pivot_longer(-Signature, names_to = "Pair", values_to = "SpearmanR") %>%
  mutate(
    Pair      = factor(Pair, levels = c("SPA vs MuSiCaL", "SPA vs SigDyn", "MuSiCaL vs SigDyn")),
    Signature = factor(Signature, levels = rev(sig_corr_df$Signature)),
    lbl       = ifelse(!is.na(SpearmanR), sprintf("%.2f", SpearmanR), "–"),
    txt_col   = ifelse(!is.na(SpearmanR) & abs(SpearmanR) > 0.5, "white", "grey20")
  )

p_corr <- ggplot(sig_corr_long, aes(x = Pair, y = Signature, fill = SpearmanR)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = lbl, colour = txt_col), size = 2.5, show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#d73027", mid = "white", high = "#1a9850",
    midpoint = 0, limits = c(-1, 1), name = "Spearman ρ", na.value = "grey85"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title    = paste0("Per-signature Activity Correlation — ", cohort_id),
    subtitle = "Spearman ρ of per-sample activities  |  signatures detected in ≥5% of samples",
    x = NULL, y = NULL
  ) +
  theme_ohchibi_pubr(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 9, face = "bold"),
    axis.text.y        = element_text(size = 8),
    legend.position    = "right",
    legend.direction   = "vertical"
  )

corr_h <- min(max(5.0, nrow(sig_corr_df) * 0.38 + 2.5), 40)
ggsave(paste0(out_prefix, "_concordance_signature_corr.png"), p_corr,
       width = 7, height = corr_h, dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("Saved signature corr     -> %s_concordance_signature_corr.png\n", out_prefix))

cat("Concordance analysis complete.\n")
