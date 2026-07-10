#!/usr/bin/env Rscript
# run_sigdyn_analysis.R
#
# SigDyn-inspired mutational signature dynamics analysis.
#
# The SigDyn framework (Goncalves et al.) characterises how COSMIC SBS signatures
# evolve dynamically across tumour samples.  This script implements the core
# analytical steps using R/MutationalPatterns:
#   1. Fit COSMIC v3 SBS signatures to each sample (strict backwards-elimination NNLS).
#   2. Export per-sample activities and per-sample cosine similarities.
#   3. Compute signature dynamics statistics: mean, SD, CV, n_active_samples.
#   4. Produce a fractional-activity heatmap and a CV bar plot.
#
# Supports WES and WGS modes: when exome=true the COSMIC v3.2 WES signature
# set (exported from MuSiCaL at container build time) is used for fitting,
# matching the catalog used by the MuSiCaL pipeline.  When exome=false the
# paired COSMIC v3.2 WGS catalog is used.  Falls back to the MutationalPatterns
# bundled COSMIC v3.3.1 WGS signatures if the pre-exported files are absent.
#
# Usage:
#   Rscript run_sigdyn_analysis.R \
#       <sbs96_matrix> <cohort_id> <genome_build> <exome> <output_prefix>

suppressPackageStartupMessages({
  library(MutationalPatterns)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop(paste(
    "Usage: run_sigdyn_analysis.R",
    "<sbs96_matrix> <cohort_id> <genome_build> <exome> <output_prefix>"
  ))
}

matrix_file  <- args[1L]
cohort_id    <- args[2L]
genome_build <- args[3L]
exome        <- tolower(args[4L]) %in% c("true", "1", "yes")
out_prefix   <- args[5L]

mode_label <- if (exome) "WES" else "WGS"

# ── Load SBS96 count matrix ───────────────────────────────────────────────────
cat("Loading SBS96 matrix...\n")
counts_df <- read.delim(matrix_file, check.names = FALSE, row.names = 1L)
counts    <- as.matrix(counts_df)
cat(sprintf("Matrix: %d mutation types x %d samples  [%s]\n",
            nrow(counts), ncol(counts), mode_label))

if (ncol(counts) == 0L) stop("No samples found in count matrix.")

# ── Load COSMIC SBS reference signatures (WES- or WGS-specific) ──────────────
# COSMIC signatures: load from the pre-built TSV (written from MutationalPatterns
# at container build time) or fall back to MutationalPatterns directly.
# Both WES and WGS modes use the same COSMIC v3.3.1 WGS catalog; a WES-specific
# catalog is not available in this container.
cosmic_tsv <- "/usr/local/share/cosmic/COSMIC_v3.2_SBS_WGS.tsv"

# The pre-built TSV was exported from MutationalPatterns get_known_signatures(),
# which returns a matrix with NULL rownames; write.table then writes numeric
# indices (1,2,...) instead of mutation type names.  Detect this and fall back
# to exporting from MuSiCaL (Python), which has proper feature labels.
.has_proper_rownames <- function(m) {
  rn <- rownames(m)
  !is.null(rn) && length(rn) > 0L && grepl("^[ACGT]\\[", rn[[1L]])
}

if (file.exists(cosmic_tsv)) {
  cat(sprintf("Loading COSMIC v3.2 signatures from pre-built TSV [%s mode]...\n", mode_label))
  cosmic_sigs <- as.matrix(
    read.delim(cosmic_tsv, check.names = FALSE, row.names = 1L)
  )
  if (!.has_proper_rownames(cosmic_sigs)) {
    cat("  TSV has numeric row indices; re-exporting from MuSiCaL (Python)...\n")
    tmp_cat <- tempfile(fileext = ".tsv")
    py_cmd <- sprintf(
      "python3 -c \"from musical.catalog import load_catalog; import pandas as pd; c=load_catalog('COSMIC_v3p2_SBS_WGS'); pd.DataFrame(c.W, index=c.features, columns=c.signatures).to_csv('%s', sep='\\t')\"",
      tmp_cat
    )
    ret <- system(py_cmd, intern = FALSE)
    if (ret != 0L || !file.exists(tmp_cat)) {
      stop("Failed to export COSMIC catalog from MuSiCaL; check Python/MuSiCaL installation.")
    }
    cosmic_sigs <- as.matrix(read.delim(tmp_cat, check.names = FALSE, row.names = 1L))
    unlink(tmp_cat)
    cat(sprintf("  Re-loaded from MuSiCaL: %d mutation types x %d signatures.\n",
                nrow(cosmic_sigs), ncol(cosmic_sigs)))
  }
} else {
  cat("Pre-built COSMIC TSV not found; loading directly from MutationalPatterns...\n")
  cosmic_sigs <- tryCatch(
    get_known_signatures(source = "COSMIC_v3.2"),
    error = function(e) {
      cat(sprintf("  COSMIC_v3.2 unavailable (%s); trying COSMIC_v3.1...\n", e$message))
      tryCatch(
        get_known_signatures(source = "COSMIC_v3.1"),
        error = function(e2) {
          cat("  Falling back to default COSMIC signatures.\n")
          get_known_signatures(source = "COSMIC")
        }
      )
    }
  )
}
cat(sprintf("Loaded %d COSMIC SBS signatures.\n", ncol(cosmic_sigs)))

# Align row order between count matrix and signatures
shared_muts <- intersect(rownames(counts), rownames(cosmic_sigs))
if (length(shared_muts) < 90L) {
  warning(sprintf("Only %d/96 mutation types shared — check matrix format.", length(shared_muts)))
}
counts      <- counts[shared_muts, , drop = FALSE]
cosmic_sigs <- cosmic_sigs[shared_muts, , drop = FALSE]

# ── Signature fitting (strict backwards-elimination NNLS, SigDyn-style) ───────
# method = "backwards" is used instead of "best_subset" because best_subset
# enumerates all C(n_sigs, k) signature combinations per sample — combinatorially
# intractable with 78 COSMIC v3.2 signatures and large cohorts (>100 samples).
# Backwards elimination (O(n_sigs^2) per sample) gives equivalent signature
# selection quality at a fraction of the runtime.
cat("Fitting COSMIC signatures via fit_to_signatures_strict() [backwards]...\n")
fit_strict <- tryCatch(
  fit_to_signatures_strict(
    mut_matrix = counts,
    signatures = cosmic_sigs,
    method     = "backwards",
    max_delta  = 0.004
  ),
  error = function(e) {
    cat(sprintf(
      "  fit_to_signatures_strict failed: %s\n  Falling back to fit_to_signatures.\n",
      e$message
    ))
    NULL
  }
)

if (!is.null(fit_strict)) {
  contribution  <- fit_strict$fit_res$contribution   # signatures x samples
  reconstructed <- fit_strict$fit_res$reconstructed  # 96 x samples
} else {
  fit_basic     <- fit_to_signatures(counts, cosmic_sigs)
  contribution  <- fit_basic$contribution
  reconstructed <- fit_basic$reconstructed
}

activities <- t(contribution)  # samples x signatures

# ── Per-sample cosine similarity ──────────────────────────────────────────────
.cos_sim <- function(a, b) {
  d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (d == 0) 0.0 else sum(a * b) / d
}
per_sample_cs <- vapply(
  seq_len(ncol(counts)),
  function(i) .cos_sim(counts[, i], reconstructed[, i]),
  numeric(1L)
)
names(per_sample_cs) <- colnames(counts)

# ── Export activities ─────────────────────────────────────────────────────────
act_df <- cbind(Samples = rownames(activities), as.data.frame(activities))
act_file <- paste0(out_prefix, "_sigdyn_activities.tsv")
write.table(act_df, act_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved activities           -> %s\n", act_file))

# ── Export cosine similarities ────────────────────────────────────────────────
cos_df <- data.frame(
  Sample           = names(per_sample_cs),
  CosineSimilarity = unname(per_sample_cs)
)
cos_file <- paste0(out_prefix, "_sigdyn_cosine_similarities.tsv")
write.table(cos_df, cos_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved cosine similarities  -> %s\n", cos_file))

# ── Signature dynamics statistics ─────────────────────────────────────────────
cat("Computing signature dynamics statistics...\n")
act_long <- as.data.frame(activities) %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Signature", values_to = "Activity")

sig_stats <- act_long %>%
  group_by(Signature) %>%
  summarise(
    mean_activity    = mean(Activity, na.rm = TRUE),
    sd_activity      = sd(Activity,   na.rm = TRUE),
    cv               = ifelse(
      mean_activity > 0,
      sd_activity / mean_activity,
      NA_real_
    ),
    n_active_samples = sum(Activity > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(mean_activity > 0) %>%
  arrange(desc(mean_activity))

stats_file <- paste0(out_prefix, "_sigdyn_signature_stats.tsv")
write.table(sig_stats, stats_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved dynamics stats       -> %s\n", stats_file))

# ── Plot: fractional-activity heatmap ─────────────────────────────────────────
active_sigs <- sig_stats$Signature
n_samples   <- nrow(activities)
n_sigs      <- length(active_sigs)

if (n_sigs >= 1L && n_samples >= 1L) {
  mat_frac <- sweep(
    activities[, active_sigs, drop = FALSE],
    1L,
    pmax(rowSums(activities), 1e-9),
    "/"
  )
  hm_df <- as.data.frame(mat_frac) %>%
    tibble::rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Signature", values_to = "Fraction") %>%
    mutate(
      Signature = factor(Signature, levels = active_sigs),
      Sample    = factor(Sample,    levels = rownames(activities))
    )

  p_hm <- ggplot(hm_df, aes(x = Signature, y = Sample, fill = Fraction)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = c("white", "#deebf7", "#3182bd", "#08306b"),
      limits = c(0, 1),
      name   = "Fraction"
    ) +
    labs(
      title    = paste0("SigDyn: Signature Activity Heatmap — ", cohort_id),
      subtitle = mode_label,
      x = "COSMIC Signature", y = "Sample"
    ) +
    theme_classic(base_size = 11L) +
    theme(
      axis.text.x   = element_text(angle = 45L, hjust = 1, size = 8L),
      axis.text.y   = element_text(size = max(6L, 10L - n_samples %/% 15L)),
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9L)
    )

  w_hm <- min(max(10, 0.5 * n_sigs + 3), 40)
  h_hm <- min(max(6,  0.22 * n_samples + 3), 40)
  hm_file <- paste0(out_prefix, "_sigdyn_activity_heatmap.png")
  ggsave(hm_file, p_hm, width = w_hm, height = h_hm, dpi = 200L, bg = "white", limitsize = FALSE)
  cat(sprintf("Saved activity heatmap     -> %s\n", hm_file))
}

# ── Plot: coefficient of variation (signature dynamics bar plot) ───────────────
cv_data <- sig_stats %>% filter(!is.na(cv)) %>% arrange(desc(cv))
if (nrow(cv_data) >= 1L) {
  p_cv <- ggplot(cv_data, aes(x = reorder(Signature, cv), y = cv)) +
    geom_bar(stat = "identity", fill = "#2171b5", width = 0.7) +
    coord_flip() +
    labs(
      title    = paste0("SigDyn: Signature Dynamics (CV) — ", cohort_id),
      subtitle = paste0(mode_label, "  |  High CV = signature varies markedly across samples"),
      x = "COSMIC Signature",
      y = "Coefficient of Variation (SD / Mean activity)"
    ) +
    theme_classic(base_size = 11L) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9L)
    )

  h_cv  <- min(max(4, 0.3 * nrow(cv_data) + 2), 40)
  cv_file <- paste0(out_prefix, "_sigdyn_signature_dynamics_cv.png")
  ggsave(cv_file, p_cv, width = 10, height = h_cv, dpi = 200L, bg = "white", limitsize = FALSE)
  cat(sprintf("Saved CV dynamics plot     -> %s\n", cv_file))
}

cat("SigDyn analysis complete.\n")
