suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
})

# Matrix-based SigProfilerAssignment: merged_signature_activities.txt + merged coverage TSV.
# Single stacked-bar figure; cohort assignment rate (assigned / input variants) is shown in the title.
# Writes signature_bargraphs_combined.png and signature_bargraphs_combined.pdf

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

compute_shared_colors <- function(sig_file) {
  sig <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(sig[, -1])
  rownames(mat) <- sig$Samples
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  long <- as.data.frame(mat) %>%
    tibble::rownames_to_column("SampleId") %>%
    pivot_longer(-SampleId, names_to = "Signature", values_to = "Mutations")
  sigs_ordered <- long %>%
    group_by(Signature) %>%
    summarise(total = sum(Mutations), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(Signature)
  make_sig_colors(sigs_ordered)
}

compute_sample_order <- function(sig_file) {
  sig <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(sig[, -1])
  rownames(mat) <- sig$Samples
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (nrow(mat) == 1) rownames(mat) else cluster_order(mat)
}

cohort_assignment_pct <- function(coverage_file) {
  cov <- read.delim(coverage_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("InputVariants", "NumberVariantsAssignedMutSig") %in% names(cov))) {
    return(NA_real_)
  }
  inp <- sum(as.numeric(cov$InputVariants), na.rm = TRUE)
  asg <- sum(as.numeric(cov$NumberVariantsAssignedMutSig), na.rm = TRUE)
  if (!is.finite(inp) || inp <= 0) return(NA_real_)
  100 * asg / inp
}

make_assignment_plot <- function(sig_file, shared_sig_cols, title, sample_order, y_max, show_legend) {
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
      axis.ticks.length.x = grid::unit(0, "pt"),
      axis.text.y         = element_text(size = 12),
      axis.title.y        = element_text(size = 13),
      panel.grid.major.y  = element_line(color = "#D9D9D9", linewidth = 0.4,
                                         linetype = "dashed"),
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      strip.background    = element_rect(fill = "transparent", color = NA),
      strip.text          = element_text(face = "bold", size = 14),
      plot.title          = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position     = if (show_legend) "right" else "none",
      legend.key.size     = grid::unit(0.55, "cm"),
      legend.text         = element_text(size = 11),
      legend.title        = element_text(size = 12, face = "bold")
    )
}

# Usage: Rscript plot_matrix_signature_bargraphs.R <merged_signature_activities.txt> <mutsig_coverage_merged.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_matrix_signature_bargraphs.R <merged_signature_activities.txt> <mutsig_coverage_merged.tsv>")
}

sig_file      <- args[1]
coverage_file <- args[2]

pct <- cohort_assignment_pct(coverage_file)
title_main <- if (is.finite(pct)) {
  sprintf(
    "Matrix-based COSMIC assignment (%.0f%% of mutations contributed to assignment)",
    pct
  )
} else {
  "Matrix-based COSMIC assignment (cohort assignment rate unavailable)"
}

shared_cols  <- compute_shared_colors(sig_file)
shared_order <- compute_sample_order(sig_file)

sig <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
y_max <- max(rowSums(sig[, -1]))

plot_main <- make_assignment_plot(
  sig_file, shared_cols,
  title = title_main,
  sample_order = shared_order,
  y_max = y_max, show_legend = TRUE
)

w_in <- 20
h_in <- 11

ggsave(
  "signature_bargraphs_combined.png",
  plot = plot_main, width = w_in, height = h_in, dpi = 300, bg = "white"
)
message("Saved: signature_bargraphs_combined.png")

# Vector PDF: base pdf() + print() — reliable without Cairo-linked ggsave(cairo_pdf) in minimal images.
pdf_path <- "signature_bargraphs_combined.pdf"
grDevices::pdf(pdf_path, width = w_in, height = h_in, useDingbats = FALSE, paper = "special")
tryCatch(
  print(plot_main),
  finally = grDevices::dev.off()
)
if (!file.exists(pdf_path) || file.info(pdf_path)$size < 1L) {
  stop("Failed to write ", pdf_path, call. = FALSE)
}
message("Saved: ", pdf_path)

message("Done.")
