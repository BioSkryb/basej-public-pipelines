#!/usr/bin/env Rscript
## Two-panel hexbin plot: Rho vs log10(Germline_qval) for the binomial/beta-binomial
## first pass (left) and Sequoia second pass (right), each faceted by filter verdict.
##
## Usage:
##   Rscript rscript_plot_rho_qval_hexbin.R <group_id> <binom_betabinom_tsv> \
##       <sequoia_filtering_all_txt> <output_pdf>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(pals)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript rscript_plot_rho_qval_hexbin.R ",
       "<group_id> <binom_betabinom_tsv> <sequoia_filtering_all_txt> <output_pdf>")
}
group_id     <- args[1]
binom_file   <- args[2]
sequoia_file <- args[3]
out_file     <- args[4]

# ── Theme ─────────────────────────────────────────────────────────────────────
argv0      <- commandArgs(trailingOnly = FALSE)
file_arg   <- argv0[grepl("^--file=", argv0)]
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg[1])) else "."
theme_candidates <- c("/usr/local/bin/theme_ohchibi_pubr.R",
                      file.path(script_dir, "theme_ohchibi_pubr.R"))
theme_path <- theme_candidates[file.exists(theme_candidates)][1]
if (!is.na(theme_path) && nzchar(theme_path)) {
  source(theme_path)
} else {
  theme_ohchibi_pubr <- ggplot2::theme_bw
}

# ── Shared sqrt-compressed axis for negative log10 q-values ──────────────────
# forward: -sqrt(-x)  expands 0 to -100 region, compresses extreme tail
sqrt_neg_trans <- scales::trans_new(
  name      = "sqrt_neg",
  transform = function(x) -sqrt(-x),
  inverse   = function(y) -(y^2),
  breaks    = function(x) {
    b <- c(0, -1, -5, -10, -25, -50, -100, -150, -200, -250, -300, -323)
    b[b >= min(x) & b <= max(x)]
  },
  domain = c(-Inf, 0)
)

pal_colors <- pals::kovesi.rainbow_bgyrm_35_85_c71(256)

# ── Panel builder ─────────────────────────────────────────────────────────────
make_hex_panel <- function(df, rho_col, qval_col, filter_col, panel_title) {
  df_plot <- df %>%
    rename(Rho          = !!rho_col,
           Germline_qval = !!qval_col,
           FilterStatus  = !!filter_col) %>%
    filter(!is.na(Rho), !is.na(Germline_qval)) %>%
    mutate(
      qval_floored = ifelse(Germline_qval == 0, .Machine$double.xmin, Germline_qval),
      log10_qval   = log10(qval_floored)
    )

  n <- nrow(df_plot)
  cat(panel_title, "— variants plotted:", n, "\n")

  ggplot(df_plot, aes(x = Rho, y = log10_qval)) +
    geom_hex(bins = 50) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      oob    = scales::squish
    ) +
    scale_y_continuous(
      trans  = sqrt_neg_trans,
      labels = function(x) as.character(x)
    ) +
    scale_fill_gradientn(colors = pal_colors, name = "count") +
    facet_wrap(~ FilterStatus, ncol = 2) +
    labs(
      x     = "Rho",
      y     = "log10(Germline q-value)  [sqrt-compressed]",
      title = sprintf("%s\nn = %s variants", panel_title, format(n, big.mark = ","))
    ) +
    theme_ohchibi_pubr()
}

# ── Left panel: binomial / beta-binomial first pass ──────────────────────────
df_binom <- read.table(binom_file, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE)
p_left <- make_hex_panel(
  df_binom,
  rho_col     = "Rho",
  qval_col    = "Germline_qval",
  filter_col  = "BinomialBetabinomialFilter",
  panel_title = paste0(group_id, " — Binomial/Beta-binomial")
)

# ── Right panel: Sequoia second pass ─────────────────────────────────────────
df_seq <- read.table(sequoia_file, sep = " ", header = TRUE,
                     stringsAsFactors = FALSE)
p_right <- make_hex_panel(
  df_seq,
  rho_col     = "Rho",
  qval_col    = "Germline_qval",
  filter_col  = "SecondPassFilter",
  panel_title = paste0(group_id, " — Sequoia second pass")
)

# ── Combine and save ──────────────────────────────────────────────────────────
p_out <- p_left + p_right + patchwork::plot_layout(ncol = 2)
ggsave(out_file, p_out, width = 16, height = 6)
cat("Saved:", out_file, "\n")
