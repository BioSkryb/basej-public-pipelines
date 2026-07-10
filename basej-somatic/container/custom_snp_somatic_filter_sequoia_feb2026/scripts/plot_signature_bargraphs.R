suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
  library(egg)
})

# ── helpers ──────────────────────────────────────────────────────────────────

bray_curtis_dist <- function(mat) {
  n <- nrow(mat)
  d <- matrix(0, n, n, dimnames = list(rownames(mat), rownames(mat)))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      xi <- mat[i, ]; xj <- mat[j, ]
      bc <- 1 - 2 * sum(pmin(xi, xj)) / (sum(xi) + sum(xj))
      d[i, j] <- d[j, i] <- bc
    }
  }
  as.dist(d)
}

cluster_order <- function(mat) {
  hc <- hclust(bray_curtis_dist(mat), method = "ward.D")
  hc$labels[hc$order]
}

make_sig_colors <- function(sigs) {
  n   <- length(sigs)
  pal <- unique(c(
    brewer.pal(12, "Set3"),
    brewer.pal(8,  "Dark2"),
    brewer.pal(9,  "Set1"),
    brewer.pal(8,  "Accent")
  ))[seq_len(n)]
  setNames(pal, sigs)
}

# ── shared color palette across both files ────────────────────────────────────

compute_shared_colors <- function(sig_files) {
  all_long <- lapply(sig_files, function(f) {
    sig <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
    mat <- as.matrix(sig[, -1])
    rownames(mat) <- sig$Samples
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
    as.data.frame(mat) %>%
      tibble::rownames_to_column("SampleId") %>%
      pivot_longer(-SampleId, names_to = "Signature", values_to = "Mutations")
  })
  sigs_ordered <- bind_rows(all_long) %>%
    group_by(Signature) %>%
    summarise(total = sum(Mutations), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(Signature)
  make_sig_colors(sigs_ordered)
}

# ── compute sample order from a sig file (Bray-Curtis + ward.D on all samples) ─

compute_sample_order <- function(sig_file) {
  sig <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(sig[, -1])
  rownames(mat) <- sig$Samples
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (nrow(mat) == 1) rownames(mat) else cluster_order(mat)
}

# ── per-file plot builder ─────────────────────────────────────────────────────

make_plots <- function(sig_file, shared_sig_cols, title,
                       sample_order, y_max, show_legend = FALSE) {

  sig <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(sig[, -1])
  rownames(mat) <- sig$Samples
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  sample_order <- sample_order[sample_order %in% rownames(mat)]

  df <- as.data.frame(mat) %>%
    tibble::rownames_to_column("SampleId") %>%
    pivot_longer(-SampleId, names_to = "Signature", values_to = "Mutations") %>%
    filter(Mutations > 0) %>%
    mutate(SampleId = factor(SampleId, levels = sample_order))

  sigs_here <- names(shared_sig_cols)[names(shared_sig_cols) %in% unique(df$Signature)]
  cols_here  <- shared_sig_cols[sigs_here]
  df <- df %>% mutate(Signature = factor(Signature, levels = sigs_here))

  ggplot(df, aes(x = SampleId, y = Mutations, fill = Signature)) +
    geom_bar(stat = "identity", width = 1,
             position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols_here, name = "Signature",
                      breaks = sigs_here, drop = FALSE) +
    scale_y_continuous(limits = c(0, y_max),
                       expand = expansion(mult = c(0, 0.03)),
                       breaks = scales::pretty_breaks(n = 6)(c(0, y_max))) +
    labs(x = NULL, y = "Number of somatic mutations", title = title) +
    theme_classic(base_size = 13) +
    theme(
      panel.background    = element_rect(fill = "white", color = NA),
      panel.border        = element_blank(),
      plot.background     = element_rect(fill = "white", color = NA),
      axis.text.x         = element_blank(),
      axis.title.x        = element_blank(),
      axis.ticks.x        = element_blank(),
      axis.ticks.length.x = unit(0, "pt"),
      axis.text.y         = element_text(size = 12),
      axis.title.y        = element_text(size = 13),
      panel.grid.major.y  = element_line(color = "#D9D9D9", linewidth = 0.4,
                                         linetype = "dashed"),
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      strip.background    = element_rect(fill = "transparent", color = NA),
      strip.text          = element_text(face = "bold", size = 14),
      plot.title          = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position     = if (show_legend) "right" else "none",
      legend.key.size     = unit(0.55, "cm"),
      legend.text         = element_text(size = 11),
      legend.title        = element_text(size = 12, face = "bold")
    )
}

# ── coverage panel builder ────────────────────────────────────────────────────

make_coverage_plot <- function(coverage_file, sample_order) {
  cov <- read.delim(coverage_file, stringsAsFactors = FALSE, check.names = FALSE)
  # Ensure every sample in sample_order is present; fill missing with 0
  df <- data.frame(SampleID = sample_order, stringsAsFactors = FALSE)
  df <- merge(df, cov[, c("SampleID", "ProportionVariantsAssignedMutSig")],
              by = "SampleID", all.x = TRUE)
  df$ProportionVariantsAssignedMutSig[is.na(df$ProportionVariantsAssignedMutSig)] <- 0
  df$SampleID <- factor(df$SampleID, levels = sample_order)
  cov <- df

  ggplot(df, aes(x = SampleID, y = ProportionVariantsAssignedMutSig)) +
    geom_bar(stat = "identity", fill = "#4D9DE0", width = 1) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = NULL,
         y = "Proportion assigned\nmutational signature",
         title = "Proportion of input variants assigned a mutational signature") +
    theme_classic(base_size = 11) +
    theme(
      panel.background    = element_rect(fill = "white", color = NA),
      panel.border        = element_blank(),
      plot.background     = element_rect(fill = "white", color = NA),
      axis.text.x         = element_blank(),
      axis.title.x        = element_blank(),
      axis.ticks.x        = element_blank(),
      axis.ticks.length.x = unit(0, "pt"),
      axis.text.y         = element_text(size = 11),
      axis.title.y        = element_text(size = 11),
      panel.grid.major.y  = element_line(color = "#D9D9D9", linewidth = 0.4,
                                         linetype = "dashed"),
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      plot.title          = element_text(hjust = 0.5, size = 13, face = "bold")
    )
}

# ── inputs ────────────────────────────────────────────────────────────────────
# Usage: Rscript plot_signature_bargraphs.R <zero_filtered.tsv> <cosine_filtered.tsv> <mutsig_coverage.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
  stop("Usage: Rscript plot_signature_bargraphs.R <zero_filtered.tsv> <cosine_filtered.tsv> <mutsig_coverage.tsv>")

sig_files     <- list(zero = args[1], cosine = args[2])
coverage_file <- args[3]

# ── build plots ───────────────────────────────────────────────────────────────

shared_cols  <- compute_shared_colors(unname(sig_files))
shared_order <- compute_sample_order(sig_files$zero)

# shared y-axis: max per-sample total mutations across both files
shared_y_max <- max(sapply(unname(sig_files), function(f) {
  sig <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  max(rowSums(sig[, -1]))
}))

plot_cov <- make_coverage_plot(coverage_file, sample_order = shared_order)

plot_A <- make_plots(sig_files$zero,   shared_cols,
                     title = "Unfiltered",      sample_order = shared_order,
                     y_max = shared_y_max, show_legend = TRUE)
plot_B <- make_plots(sig_files$cosine, shared_cols,
                     title = "Cosine filtered", sample_order = shared_order,
                     y_max = shared_y_max, show_legend = FALSE)

combined <- egg::ggarrange(plot_cov, plot_A, plot_B,
                           ncol = 1, heights = c(0.5, 1, 1),
                           labels = c("A", "B", "C"))

ggsave("signature_bargraphs_combined.png",
       plot = combined, width = 20, height = 18, dpi = 300, bg = "white")
message("Saved: signature_bargraphs_combined.png")
message("Done.")
