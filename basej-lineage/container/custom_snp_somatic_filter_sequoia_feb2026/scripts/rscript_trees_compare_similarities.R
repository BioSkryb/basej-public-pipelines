#!/usr/bin/env Rscript
# rscript_trees_compare_similarities.R
#
# All-vs-all comparison of phylogenetic trees produced across filtering schemes.
# Processes tree types snv / indel / both in a single run and writes:
#   - Per-type individual TSV tables and PDF heatmaps
#   - One master PDF combining all plots with SNV / INDEL / BOTH section headers
#
# Five pairwise metrics per type (when ≥ 2 branch-length trees are found):
#   1. Cophenetic distance correlation (Pearson)   — branch-length trees
#   2. Normalised Robinson-Foulds (RF) distance    — branch-length trees
#   3. Branch score distance (topology + lengths)  — branch-length trees
#   4. Baker's Gamma dendrogram similarity         — hclust on cophenetic distances
#   Metrics 1–2 are also repeated for consensus .contree files.
#   5. Tanglegrams (all-vs-all, multi-page)        — branch-length trees

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(ggplot2)
  library(reshape2)
})

option_list <- list(
  make_option(c("-i", "--donor_id"),    type = "character",
              help = "Donor / patient ID used as filename prefix"),
  make_option(c("-t", "--tree_types"),  type = "character", default = "snv,indel,both",
              help = "Comma-separated variant types to compare [default: snv,indel,both]"),
  make_option(c("-d", "--input_dir"),   type = "character", default = ".",
              help = "Directory (searched recursively) for .tree/.contree files [default: .]"),
  make_option(c("-o", "--output_dir"),  type = "character", default = ".",
              help = "Directory for output files [default: .]")
)

opt        <- parse_args(OptionParser(option_list = option_list))
donor_id   <- opt$donor_id
tree_types <- strsplit(opt$tree_types, ",")[[1]]
input_dir  <- opt$input_dir
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────

extract_label <- function(f, donor, type) {
  bn <- basename(f)
  m  <- regmatches(bn, regexec(
    sprintf("^%s_(.+)_%s_tree_with_branch_length_selectedscheme\\.tree$", donor, type), bn))[[1]]
  if (length(m) >= 2) return(m[2])
  m2 <- regmatches(bn, regexec(
    sprintf("^%s_(.+)_%s_for_MPBoot\\.fa\\.contree$", donor, type), bn))[[1]]
  if (length(m2) >= 2) return(m2[2])
  tools::file_path_sans_ext(bn)
}

read_named_trees <- function(files, donor, type) {
  result <- list()
  for (f in files) {
    t <- tryCatch(read.tree(f), error = function(e) { warning(e$message); NULL })
    if (is.null(t)) next
    result[[extract_label(f, donor, type)]] <- t
  }
  result
}

prune_to_common <- function(t1, t2) {
  tips <- intersect(t1$tip.label, t2$tip.label)
  list(t1 = drop.tip(t1, setdiff(t1$tip.label, tips)),
       t2 = drop.tip(t2, setdiff(t2$tip.label, tips)),
       tips = tips)
}

coph_cor <- function(t1, t2) {
  p <- prune_to_common(t1, t2)
  if (length(p$tips) < 3) return(NA_real_)
  d1 <- cophenetic.phylo(p$t1)[p$tips, p$tips]
  d2 <- cophenetic.phylo(p$t2)[p$tips, p$tips]
  cor(d1[lower.tri(d1)], d2[lower.tri(d2)], method = "pearson")
}

norm_rf <- function(t1, t2) {
  p      <- prune_to_common(t1, t2)
  if (length(p$tips) < 4) return(NA_real_)
  rf_raw <- as.numeric(dist.topo(p$t1, p$t2, method = "PH85"))
  max_rf <- 2 * (length(p$tips) - 3)
  if (max_rf == 0) return(0)
  rf_raw / max_rf
}

branch_score <- function(t1, t2) {
  p    <- prune_to_common(t1, t2)
  if (length(p$tips) < 4) return(NA_real_)
  bs   <- as.numeric(dist.topo(p$t1, p$t2, method = "score"))
  norm <- mean(c(sum(p$t1$edge.length), sum(p$t2$edge.length)))
  if (norm == 0) return(0)
  bs / norm
}

bakers_gamma <- function(h1, h2) {
  tips <- intersect(h1$labels, h2$labels)
  if (length(tips) < 3) return(NA_real_)
  m1 <- as.matrix(cophenetic(h1))[tips, tips]
  m2 <- as.matrix(cophenetic(h2))[tips, tips]
  cor(m1[lower.tri(m1)], m2[lower.tri(m2)], method = "spearman")
}

# All cells initialised to NA; only diagonal set to diag_val.
# Lower triangle computed, then mirrored to upper triangle.
pairwise_matrix <- function(trees, FUN, diag_val = NA_real_) {
  nm  <- names(trees)
  n   <- length(nm)
  mat <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
  diag(mat) <- diag_val
  for (i in seq_len(n - 1))
    for (j in (i + 1):n) {
      v <- tryCatch(FUN(trees[[nm[i]]], trees[[nm[j]]]), error = function(e) NA_real_)
      mat[i, j] <- v; mat[j, i] <- v
    }
  mat
}

# Per-metric interpretation captions shown at the bottom of every heatmap.
METRIC_LEGENDS <- list(
  coph_bl = paste(
    "Cophenetic Correlation (branch-length trees):",
    "Pearson correlation of all pairwise cophenetic distances between two trees (Sokal & Rohlf 1962).",
    "Range -1 to +1. Values closer to +1 indicate more similar evolutionary distances between samples;",
    "values above 0.9 are commonly interpreted as strong concordance. Diagonal = 1 (self-comparison).",
    "Sensitive to large branch-length outliers as it is Pearson-based."
  ),
  rf_bl = paste(
    "Normalised Robinson-Foulds Distance (branch-length trees):",
    "Fraction of bipartitions (splits) that differ between two trees, divided by the maximum possible",
    "RF distance [ 2*(n-3) for n taxa ]. Range 0-1. Values closer to 0 = more similar topology;",
    "0 = topologically identical, 1 = no splits shared. Values < 0.1 typically indicate very similar",
    "trees; values > 0.5 indicate substantial topological disagreement. Diagonal = 0."
  ),
  bs_bl = paste(
    "Branch Score Distance (branch-length trees):",
    "Based on Kuhner & Felsenstein (1994): sum of squared branch-length differences over all",
    "bipartitions in the union of both trees (absent bipartitions set to 0), normalised here by",
    "the mean total branch length of the two trees (implementation-specific normalisation).",
    "Lower values indicate more similar trees (topology + branch lengths); 0 = identical.",
    "Not bounded above; values are comparable only within the same dataset."
  ),
  coph_ct = paste(
    "Cophenetic Correlation (consensus trees):",
    "Same Pearson cophenetic correlation as above, applied to majority-rule consensus trees",
    "(only splits with > 50% bootstrap support are retained; rest collapsed to polytomies).",
    "Range -1 to +1; closer to +1 = more similar. Values are typically lower than for branch-length",
    "trees because consensus trees capture only well-supported splits — marginal splits near the",
    "50% bootstrap threshold are sensitive to which variants are included and can flip in/out."
  ),
  rf_ct = paste(
    "Normalised Robinson-Foulds Distance (consensus trees):",
    "RF distance applied to majority-rule consensus trees. Range 0-1; closer to 0 = more similar.",
    "Values are typically higher than branch-length RF because consensus trees retain only",
    "strongly supported splits: a split just above 50% in one scheme may drop below 50% in another,",
    "registering as a full topological disagreement even if the underlying signal barely changed.",
    "Collapsed (polytomous) nodes do not contribute splits, making the comparison asymmetric."
  ),
  bg = paste(
    "Baker's Gamma Dendrogram Similarity (Baker 1974):",
    "Goodman-Kruskal gamma rank correlation of the merge-order ranks of all taxon pairs across",
    "two hierarchical clusterings (UPGMA on cophenetic distances). Range -1 to +1.",
    "Values closer to +1 indicate more concordant dendrogram structure; near 0 = no association;",
    "negative = inverse rank ordering. Insensitive to absolute branch lengths (rank-based only).",
    "No universally fixed threshold; permutation testing is recommended for significance. Diagonal = 1."
  ),
  tanglegram = paste(
    "Tanglegram: two rooted phylogenies drawn facing each other with matching tip labels connected by lines.",
    "Crossing lines indicate topological disagreement — tips that are clustered together in one tree",
    "are separated in the other. More crossings = less similar topology.",
    "nRF = normalised Robinson-Foulds distance (0 = identical topology, 1 = no splits shared).",
    "Branch score = Kuhner-Felsenstein distance normalised by mean total branch length (lower = more similar)."
  )
)

make_heatmap <- function(mat, title, caption = NULL) {
  df           <- reshape2::melt(mat, varnames = c("Tree1", "Tree2"), value.name = "Value")
  df$Tree1     <- factor(df$Tree1, levels = rev(rownames(mat)))
  df$Tree2     <- factor(df$Tree2, levels = colnames(mat))
  df$text_col  <- ifelse(!is.na(df$Value) & df$Value > 0.5, "white", "grey15")
  ggplot(df, aes(Tree2, Tree1, fill = Value)) +
    geom_tile(color = "grey80", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.3f", Value), color = text_col),
              size = 4.5, na.rm = TRUE) +
    scale_color_identity() +
    scale_fill_gradientn(
      colours  = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      limits   = c(0, 1),
      oob      = scales::squish,
      na.value = "grey90") +
    labs(title = title, x = NULL, y = NULL, fill = NULL, caption = caption) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          plot.title   = element_text(size = 11, face = "bold"),
          plot.caption = element_text(size = 7.5, color = "grey40",
                                      hjust = 0, margin = margin(t = 8)),
          panel.grid   = element_blank())
}

write_mat_tsv <- function(mat, path) {
  df <- cbind(label = rownames(mat), as.data.frame(round(mat, 4)))
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Written: %s\n", path))
}

# Section header page (base graphics — compatible with the master PDF device)
print_section_header <- function(label) {
  old_par <- par(mar = c(0, 0, 0, 0), bg = "white")
  plot.new()
  rect(0, 0, 1, 1, col = "#2C3E50", border = NA)
  text(0.5, 0.55, label,    col = "white", cex = 3.5, font = 2, adj = c(0.5, 0.5))
  text(0.5, 0.35, donor_id, col = "#BDC3C7", cex = 1.6, adj = c(0.5, 0.5))
  par(old_par)
}

# ── Main loop: process each type ──────────────────────────────────────────────
# Store per-type results so we can render the master PDF afterwards.
type_results <- list()

for (tree_type in tree_types) {
  cat(sprintf("\n[%s:%s] Processing ...\n", donor_id, tree_type))
  type_dir <- file.path(output_dir, sprintf("tree_comparison_%s", tree_type))
  tryCatch({   # ── per-type error boundary ──────────────────────────────────
  dir.create(type_dir, showWarnings = FALSE, recursive = TRUE)
  prefix   <- file.path(type_dir, sprintf("%s_%s", donor_id, tree_type))

  # Discover files
  tree_files <- list.files(input_dir,
    pattern = sprintf("_%s_tree_with_branch_length_selectedscheme\\.tree$", tree_type),
    recursive = TRUE, full.names = TRUE)
  contree_files <- list.files(input_dir,
    pattern = sprintf("_%s_for_MPBoot\\.fa\\.contree$", tree_type),
    recursive = TRUE, full.names = TRUE)

  cat(sprintf("[%s:%s] Found %d branch-length trees, %d contrees\n",
              donor_id, tree_type, length(tree_files), length(contree_files)))

  skip_type <- function(reason) {
    cat(sprintf("[%s:%s] Skipping — %s\n", donor_id, tree_type, reason))
    writeLines(sprintf("Skipped: %s\n", reason),
               file.path(type_dir, sprintf("%s_%s_skipped.txt", donor_id, tree_type)))
    type_results[[tree_type]] <<- list(skipped = TRUE, skip_reason = reason,
                                       type_dir = type_dir)
  }

  if (length(tree_files) < 2) {
    skip_type(sprintf("found %d branch-length tree file(s) on disk", length(tree_files)))
    next
  }

  trees_bl <- read_named_trees(tree_files,    donor_id, tree_type)
  trees_ct <- read_named_trees(contree_files, donor_id, tree_type)

  if (length(trees_bl) < 2) {
    skip_type(sprintf(
      "%d file(s) found on disk but only %d loaded successfully (check warnings above)",
      length(tree_files), length(trees_bl)))
    next
  }

  # Compute matrices
  cat(sprintf("[%s:%s] Computing pairwise metrics ...\n", donor_id, tree_type))
  mat_coph_bl <- pairwise_matrix(trees_bl, coph_cor,     diag_val = 1)
  mat_rf_bl   <- pairwise_matrix(trees_bl, norm_rf,      diag_val = 0)
  mat_bs_bl   <- pairwise_matrix(trees_bl, branch_score, diag_val = 0)
  mat_coph_ct <- if (length(trees_ct) >= 2) pairwise_matrix(trees_ct, coph_cor, diag_val = 1) else NULL
  mat_rf_ct   <- if (length(trees_ct) >= 2) pairwise_matrix(trees_ct, norm_rf,  diag_val = 0) else NULL
  mat_bg      <- if (length(trees_bl) >= 3) {
    hc_list <- lapply(trees_bl, function(t) {
      hclust(as.dist(cophenetic.phylo(t)[t$tip.label, t$tip.label]), method = "average")
    })
    pairwise_matrix(hc_list, bakers_gamma, diag_val = 1)
  } else NULL

  # TSVs
  write_mat_tsv(mat_coph_bl, sprintf("%s_cophenetic_cor_branch_length.tsv",   prefix))
  write_mat_tsv(mat_rf_bl,   sprintf("%s_normalized_RF_branch_length.tsv",    prefix))
  write_mat_tsv(mat_bs_bl,   sprintf("%s_branch_score_branch_length.tsv",     prefix))
  if (!is.null(mat_coph_ct)) {
    write_mat_tsv(mat_coph_ct, sprintf("%s_cophenetic_cor_contree.tsv", prefix))
    write_mat_tsv(mat_rf_ct,   sprintf("%s_normalized_RF_contree.tsv",  prefix))
  }
  if (!is.null(mat_bg))
    write_mat_tsv(mat_bg, sprintf("%s_bakers_gamma_dendrogram.tsv", prefix))

  # Build ggplot objects
  tt <- toupper(tree_type)
  plots <- list(
    coph_bl = make_heatmap(mat_coph_bl,
      sprintf("%s %s — Cophenetic Correlation (branch-length trees)", donor_id, tt),
      caption = METRIC_LEGENDS$coph_bl),
    rf_bl   = make_heatmap(mat_rf_bl,
      sprintf("%s %s — Normalised RF Distance (branch-length trees)", donor_id, tt),
      caption = METRIC_LEGENDS$rf_bl),
    bs_bl   = make_heatmap(mat_bs_bl,
      sprintf("%s %s — Branch Score Distance (branch-length trees)", donor_id, tt),
      caption = METRIC_LEGENDS$bs_bl)
  )
  if (!is.null(mat_coph_ct)) {
    plots$coph_ct <- make_heatmap(mat_coph_ct,
      sprintf("%s %s — Cophenetic Correlation (consensus trees)", donor_id, tt),
      caption = METRIC_LEGENDS$coph_ct)
    plots$rf_ct   <- make_heatmap(mat_rf_ct,
      sprintf("%s %s — Normalised RF Distance (consensus trees)", donor_id, tt),
      caption = METRIC_LEGENDS$rf_ct)
  }
  if (!is.null(mat_bg))
    plots$bg <- make_heatmap(mat_bg,
      sprintf("%s %s — Baker's Gamma Dendrogram Similarity", donor_id, tt),
      caption = METRIC_LEGENDS$bg)

  # Write individual heatmap PDFs
  n <- nrow(mat_coph_bl)
  w <- max(6, n * 1.4); h <- max(5, n * 1.2)
  for (nm in names(plots)) {
    fname <- switch(nm,
      coph_bl = sprintf("%s_cophenetic_cor_branch_length.pdf", prefix),
      rf_bl   = sprintf("%s_normalized_RF_branch_length.pdf",  prefix),
      bs_bl   = sprintf("%s_branch_score_branch_length.pdf",   prefix),
      coph_ct = sprintf("%s_cophenetic_cor_contree.pdf",       prefix),
      rf_ct   = sprintf("%s_normalized_RF_contree.pdf",        prefix),
      bg      = sprintf("%s_bakers_gamma_dendrogram.pdf",      prefix)
    )
    ggsave(fname, plots[[nm]], width = w, height = h, limitsize = FALSE)
    cat(sprintf("  Written: %s\n", fname))
  }

  # Build rooted trees for tanglegrams (root all by same outgroup)
  all_tips      <- Reduce(intersect, lapply(trees_bl, `[[`, "tip.label"))
  outgroup      <- sort(all_tips)[1]
  trees_rooted  <- lapply(trees_bl, function(t) {
    ladderize(root(t, outgroup = outgroup, resolve.root = TRUE))
  })

  # Write individual tanglegram PDF
  tangle_path <- sprintf("%s_tanglegrams.pdf", prefix)
  nm_r <- names(trees_rooted); n_r <- length(nm_r)
  pdf(tangle_path, width = 22, height = 12)
  for (i in seq_len(n_r - 1)) {
    for (j in (i + 1):n_r) {
      tips  <- intersect(trees_rooted[[nm_r[i]]]$tip.label,
                         trees_rooted[[nm_r[j]]]$tip.label)
      rf_v  <- tryCatch(norm_rf(trees_bl[[nm_r[i]]], trees_bl[[nm_r[j]]]),
                        error = function(e) NA_real_)
      bs_v  <- tryCatch(branch_score(trees_bl[[nm_r[i]]], trees_bl[[nm_r[j]]]),
                        error = function(e) NA_real_)
      tryCatch({
        cophyloplot(
          drop.tip(trees_rooted[[nm_r[i]]], setdiff(trees_rooted[[nm_r[i]]]$tip.label, tips)),
          drop.tip(trees_rooted[[nm_r[j]]], setdiff(trees_rooted[[nm_r[j]]]$tip.label, tips)),
          assoc = cbind(tips, tips), show.tip.label = TRUE, font = 1,
          cex = 0.55, space = 35, col = "steelblue", lwd = 0.6,
          main = sprintf("%s  [%s]  vs  [%s]\nnRF = %.3f   branch score = %.3f",
                         donor_id, nm_r[i], nm_r[j],
                         ifelse(is.na(rf_v), NaN, rf_v),
                         ifelse(is.na(bs_v), NaN, bs_v)))
        mtext(METRIC_LEGENDS$tanglegram, side = 1, line = 3,
              cex = 0.65, col = "grey40", adj = 0)
      },
        error = function(e) { plot.new(); title(sprintf("Tanglegram failed: %s", e$message)) }
      )
    }
  }
  dev.off()
  cat(sprintf("  Written: %s\n", tangle_path))

  # Stdout summary
  cat(sprintf("\n=== Summary: %s [%s] ===\n", donor_id, tree_type))
  cat("Cophenetic correlation (branch-length):\n"); print(round(mat_coph_bl, 3))
  cat("\nNormalised RF distance (branch-length):\n"); print(round(mat_rf_bl, 3))
  cat("\nBranch score distance (branch-length):\n");  print(round(mat_bs_bl, 3))
  if (!is.null(mat_bg)) { cat("\nBaker's Gamma:\n"); print(round(mat_bg, 3)) }

  type_results[[tree_type]] <- list(
    skipped      = FALSE,
    type_dir     = type_dir,
    plots        = plots,
    trees_bl     = trees_bl,
    trees_rooted = trees_rooted,
    mat_rf_bl    = mat_rf_bl,
    mat_bs_bl    = mat_bs_bl,
    plot_dims    = c(w = w, h = h),
    n_pairs      = n_r * (n_r - 1) / 2
  )

  }, error = function(e) {   # ── per-type error handler ──────────────────────
    msg <- sprintf("Unexpected error processing type '%s': %s", tree_type, e$message)
    cat(sprintf("[%s:%s] ERROR — %s\n", donor_id, tree_type, e$message))
    writeLines(msg, file.path(type_dir,
               sprintf("%s_%s_error.txt", donor_id, tree_type)))
    type_results[[tree_type]] <<- list(skipped = TRUE,
                                       skip_reason = msg,
                                       type_dir    = type_dir)
  })   # ── end per-type tryCatch ──────────────────────────────────────────────
}

# ── Master PDF ────────────────────────────────────────────────────────────────
master_path <- file.path(output_dir,
                         sprintf("%s_tree_comparison_master.pdf", donor_id))
cat(sprintf("\n[%s] Writing master PDF: %s\n", donor_id, master_path))

# Use a wide page that accommodates both heatmaps and tanglegrams.
# Open device inside tryCatch so any device-open failure is caught cleanly.
tryCatch(
  pdf(master_path, width = 22, height = 14),
  error = function(e) stop(sprintf("Cannot open master PDF device: %s", e$message))
)

# Report title page
old_par <- par(mar = c(0, 0, 0, 0), bg = "white")
plot.new()
rect(0, 0, 1, 1, col = "#1A252F", border = NA)
text(0.5, 0.60, "Tree Comparison Report",    col = "white",   cex = 4,   font = 2)
text(0.5, 0.45, donor_id,                    col = "#85C1E9", cex = 2.8, font = 1)
text(0.5, 0.28, paste("Types:", paste(tree_types, collapse = " / ")),
     col = "#BDC3C7", cex = 1.8)
text(0.5, 0.18, format(Sys.time(), "%Y-%m-%d"),
     col = "#7F8C8D", cex = 1.4)
par(old_par)

for (tree_type in tree_types) {
  res <- type_results[[tree_type]]

  # Section header
  print_section_header(toupper(tree_type))

  if (isTRUE(res$skipped)) {
    old_par <- par(mar = c(2, 2, 2, 2))
    plot.new()
    reason <- if (!is.null(res$skip_reason)) res$skip_reason
              else sprintf("no trees available for type '%s'", tree_type)
    text(0.5, 0.55, sprintf("Skipped: %s", toupper(tree_type)), cex = 2,   col = "grey30", font = 2)
    text(0.5, 0.40, reason,                                     cex = 1.2, col = "grey50")
    par(old_par)
    next
  }

  # Heatmaps — each wrapped in tryCatch so one bad plot can't break the PDF
  for (nm in names(res$plots)) {
    tryCatch(
      print(res$plots[[nm]]),
      error = function(e) {
        plot.new()
        title(sprintf("Plot failed [%s %s]: %s", tree_type, nm, e$message))
      }
    )
  }

  # Tanglegrams — re-render to the active device
  nm_r <- names(res$trees_rooted); n_r <- length(nm_r)
  if (n_r >= 2) {
    plot.new()
    text(0.5, 0.5,
         sprintf("%s   %s — Tanglegrams (%d pairs)",
                 donor_id, toupper(tree_type), res$n_pairs),
         cex = 2, font = 2, col = "#2C3E50")

    for (i in seq_len(n_r - 1)) {
      for (j in (i + 1):n_r) {
        tips <- intersect(res$trees_rooted[[nm_r[i]]]$tip.label,
                          res$trees_rooted[[nm_r[j]]]$tip.label)
        rf_v <- tryCatch(norm_rf(res$trees_bl[[nm_r[i]]], res$trees_bl[[nm_r[j]]]),
                         error = function(e) NA_real_)
        bs_v <- tryCatch(branch_score(res$trees_bl[[nm_r[i]]], res$trees_bl[[nm_r[j]]]),
                         error = function(e) NA_real_)
        tryCatch({
          cophyloplot(
            drop.tip(res$trees_rooted[[nm_r[i]]],
                     setdiff(res$trees_rooted[[nm_r[i]]]$tip.label, tips)),
            drop.tip(res$trees_rooted[[nm_r[j]]],
                     setdiff(res$trees_rooted[[nm_r[j]]]$tip.label, tips)),
            assoc = cbind(tips, tips), show.tip.label = TRUE, font = 1,
            cex = 0.55, space = 35, col = "steelblue", lwd = 0.6,
            main = sprintf("[%s]  vs  [%s]   |   %s %s\nnRF = %.3f   branch score = %.3f",
                           nm_r[i], nm_r[j], donor_id, toupper(tree_type),
                           ifelse(is.na(rf_v), NaN, rf_v),
                           ifelse(is.na(bs_v), NaN, bs_v)))
          mtext(METRIC_LEGENDS$tanglegram, side = 1, line = 3,
                cex = 0.65, col = "grey40", adj = 0)
        },
          error = function(e) {
            plot.new()
            title(sprintf("Tanglegram failed [%s vs %s]: %s", nm_r[i], nm_r[j], e$message))
          }
        )
      }
    }
  }
}

tryCatch(dev.off(), error = function(e)
  cat(sprintf("[%s] WARNING: dev.off() error: %s\n", donor_id, e$message)))
cat(sprintf("[%s] Master PDF complete: %s\n", donor_id, master_path))
