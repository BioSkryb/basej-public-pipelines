suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(jsonlite)
})

# ── helpers (identical to plot_signature_bargraphs.R) ─────────────────────────

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

compute_sample_order <- function(sig_file) {
  sig <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(sig[, -1])
  rownames(mat) <- sig$Samples
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (nrow(mat) == 1) rownames(mat) else cluster_order(mat)
}

prepare_panel_data <- function(sig_file, shared_sig_cols, sample_order) {
  sig <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(sig[, -1])
  rownames(mat) <- sig$Samples
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  sample_order_here <- sample_order[sample_order %in% rownames(mat)]
  sig_names <- names(shared_sig_cols)

  lapply(sample_order_here, function(s) {
    sigs_present <- intersect(sig_names, colnames(mat))
    mut_vals     <- setNames(as.numeric(mat[s, sigs_present]), sigs_present)
    mut_vals     <- mut_vals[mut_vals > 0]
    list(
      sample    = s,
      mutations = if (length(mut_vals) > 0) as.list(mut_vals)
                  else setNames(list(), character(0))
    )
  })
}

# load cohort_signature_cosine_similarities.tsv → nested list: sample -> sig -> value
load_cosine_similarities <- function(cosine_sim_file) {
  df <- read.delim(cosine_sim_file, check.names = FALSE, stringsAsFactors = FALSE)
  # first column is Sample
  setNames(lapply(seq_len(nrow(df)), function(i) {
    vals <- as.numeric(df[i, -1])
    names(vals) <- colnames(df)[-1]
    vals <- round(vals[!is.na(vals)], 4)
    as.list(vals)
  }), df[[1]])
}

# ── inputs ────────────────────────────────────────────────────────────────────
# Usage: Rscript generate_interactive_bargraph.R <zero_filtered.tsv> <cosine_filtered.tsv> <cosine_similarities.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
  stop("Usage: Rscript generate_interactive_bargraph.R <zero_filtered.tsv> <cosine_filtered.tsv> <cosine_similarities.tsv>")

sig_files <- list(zero = args[1], cosine = args[2])
cosine_sim_file <- args[3]

# ── main ──────────────────────────────────────────────────────────────────────

message("Computing shared colours and sample order...")
shared_cols  <- compute_shared_colors(unname(sig_files))
shared_order <- compute_sample_order(sig_files$zero)

message("Preparing panel data...")
panel_A <- prepare_panel_data(sig_files$zero,   shared_cols, shared_order)
panel_B <- prepare_panel_data(sig_files$cosine, shared_cols, shared_order)

message("Loading cosine similarities...")
cosine_sim <- load_cosine_similarities(cosine_sim_file)

signatures_list <- lapply(seq_along(shared_cols), function(i) {
  list(name = names(shared_cols)[i], color = unname(shared_cols)[i])
})

all_totals <- c(
  sapply(panel_A, function(s) sum(unlist(s$mutations))),
  sapply(panel_B, function(s) sum(unlist(s$mutations)))
)
y_max <- ceiling(max(all_totals) / 5) * 5

chart_data <- list(
  panels = list(
    list(id = "A", title = "Unfiltered",      samples = panel_A),
    list(id = "B", title = "Cosine filtered",  samples = panel_B)
  ),
  signatures = signatures_list,
  yMax       = y_max,
  cosineSim  = cosine_sim
)

json_str <- toJSON(chart_data, auto_unbox = TRUE, pretty = FALSE, null = "null")
message(paste("JSON size:", nchar(json_str), "chars"))

# ── inject JSON into HTML template ───────────────────────────────────────────

template_file <- "/usr/local/bin/signature_bargraphs_interactive_template.html"
template      <- paste(readLines(template_file, warn = FALSE), collapse = "\n")
html          <- sub("{{CHART_DATA}}", json_str, template, fixed = TRUE)

out_file <- "signature_bargraphs_interactive.html"
writeLines(html, out_file, useBytes = FALSE)
message("Saved: ", out_file)
message("Done.")
