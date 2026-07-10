#!/usr/bin/env Rscript
# run_musicatk_analysis.R
#
# musicatk (Campbell lab, campbio/musicatk) de novo mutational signature
# discovery from a cohort of somatic variant calls.
#
# Workflow:
#   1. Combine per-sample mutsig variant TSV files (Project/Sample/chrom/pos/ref/alt).
#   2. Load the appropriate BSgenome (GRCh38 → hg38, GRCh37 → hg19).
#   3. build_musica()  — create musicatk MusicA object from the variant table.
#   4. build_table()   — compute SBS96 trinucleotide count table.
#   5. discover_signatures() — de novo NMF/LDA decomposition with k_denovo components.
#   6. Export per-sample exposure activities, de novo signature profiles, and plots.
#
# Supports WES and WGS modes (passed for labelling; both use the same NMF step).
#
# Usage:
#   Rscript run_musicatk_analysis.R \
#       <tsv_glob_pattern> <cohort_id> <genome_build> <exome> <k_denovo> <output_prefix>
#
# tsv_glob_pattern   Shell glob or space-separated list captured as a single string
#                    of *_mutsig_somatic_variant_table.tsv files produced by
#                    EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE[_FROM_BINARY_MATRIX].
# k_denovo           Number of de novo signatures to extract (integer >= 2).

suppressPackageStartupMessages({
  library(musicatk)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6L) {
  stop(paste(
    "Usage: run_musicatk_analysis.R",
    "<tsv_list> <cohort_id> <genome_build> <exome> <k_denovo> <output_prefix>"
  ))
}

tsv_arg      <- args[1L]   # space-separated list of TSV paths staged into work dir
cohort_id    <- args[2L]
genome_build <- args[3L]
exome        <- tolower(args[4L]) %in% c("true", "1", "yes")
k_denovo     <- max(2L, as.integer(args[5L]))
out_prefix   <- args[6L]

mode_label <- if (exome) "WES" else "WGS"
cat(sprintf("musicatk — cohort: %s  genome: %s  mode: %s  k: %d\n",
            cohort_id, genome_build, mode_label, k_denovo))

# ── Collect all per-sample variant TSVs ──────────────────────────────────────
# The Nextflow process stages them into the work directory; we find them via glob.
tsv_files <- Sys.glob("*_mutsig_somatic_variant_table.tsv")
if (length(tsv_files) == 0L) {
  # Fallback: parse space-separated argument
  tsv_files <- strsplit(trimws(tsv_arg), "\\s+")[[1L]]
  tsv_files <- tsv_files[file.exists(tsv_files)]
}
if (length(tsv_files) == 0L) stop("No *_mutsig_somatic_variant_table.tsv files found.")
cat(sprintf("Found %d per-sample variant TSV files.\n", length(tsv_files)))

# ── Read and combine variant tables ──────────────────────────────────────────
# Expected columns: Project  Sample  ID  Genome  mut_type  chrom  pos_start  pos_end  ref  alt  Type
read_tsv_safe <- function(f) {
  tryCatch(
    read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      cat(sprintf("  WARNING: could not read %s — %s\n", f, e$message))
      NULL
    }
  )
}

all_variants <- dplyr::bind_rows(lapply(tsv_files, read_tsv_safe))
all_variants <- all_variants[!is.na(all_variants$chrom), ]

# Keep only SNVs (musicatk SBS96 table requires single-nucleotide substitutions)
all_variants <- all_variants[all_variants$mut_type %in% c("SNP", "SNV"), ]
cat(sprintf("Total SNV variants across cohort: %d\n", nrow(all_variants)))

if (nrow(all_variants) == 0L) {
  stop("No SNV variants found after filtering. Cannot run musicatk.")
}

# Ensure chr prefix (UCSC style required by BSgenome.UCSC packages)
if (!grepl("^chr", all_variants$chrom[1L])) {
  all_variants$chrom <- paste0("chr", all_variants$chrom)
}

# Build musicatk-compatible data frame
# musicatk build_musica() expects: chr, start, end, ref, alt, sample_id  (case-sensitive)
mut_df <- data.frame(
  chr       = all_variants$chrom,
  start     = as.integer(all_variants$pos_start),
  end       = as.integer(all_variants$pos_end),
  ref       = toupper(as.character(all_variants$ref)),
  alt       = toupper(as.character(all_variants$alt)),
  sample_id = as.character(all_variants$Sample),
  stringsAsFactors = FALSE
)
# Remove any rows with missing or multi-allelic data
mut_df <- mut_df[
  nchar(mut_df$ref) == 1L & nchar(mut_df$alt) == 1L &
  !is.na(mut_df$start) & !is.na(mut_df$end),
]
cat(sprintf("Retained %d SNVs after quality filter.\n", nrow(mut_df)))

n_samples <- length(unique(mut_df$sample_id))
cat(sprintf("Samples in variant table: %d\n", n_samples))

if (k_denovo > n_samples) {
  k_denovo <- max(2L, n_samples - 1L)
  cat(sprintf("  k_denovo reduced to %d (fewer samples than k requested).\n", k_denovo))
}

# ── Load BSgenome reference ───────────────────────────────────────────────────
cat(sprintf("Loading BSgenome for %s...\n", genome_build))
g <- tryCatch({
  if (genome_build %in% c("GRCh38", "hg38")) {
    library(BSgenome.Hsapiens.UCSC.hg38)
    BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (genome_build %in% c("GRCh37", "hg19", "hg37")) {
    library(BSgenome.Hsapiens.UCSC.hg19)
    BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else {
    stop(paste("Unsupported genome_build:", genome_build))
  }
}, error = function(e) {
  stop(paste("Could not load BSgenome:", e$message,
             "\nEnsure BSgenome.Hsapiens.UCSC.hg38/hg19 is installed."))
})

# ── Build musicatk MusicA object ──────────────────────────────────────────────
cat("Building musicatk MusicA object (build_musica)...\n")
musica <- build_musica(x = mut_df, g = g)

# Build SBS96 trinucleotide count table
cat("Building SBS96 count table...\n")
build_table(musica, "SBS96")
cat(sprintf("SBS96 table built for %d samples.\n", ncol(musica@mut_tables$SBS96)))

# ── De novo signature discovery ───────────────────────────────────────────────
cat(sprintf("Running de novo discovery (k = %d, LDA algorithm, 10 restarts)...\n", k_denovo))
set.seed(42L)
discover_signatures(
  musica         = musica,
  modality       = "SBS96",
  num_signatures = k_denovo,
  algorithm      = "lda",
  result_name    = "de_novo",
  nstart         = 10L,
  seed           = 42L
)

# ── Extract results ───────────────────────────────────────────────────────────
exp_mat <- exposures(musica, "de_novo", "SBS96")   # K x N
sig_mat <- signatures(musica, "de_novo", "SBS96")  # 96 x K

# Activities: samples x signatures
act_df <- cbind(
  Samples = colnames(exp_mat),
  as.data.frame(t(exp_mat))
)
act_file <- paste0(out_prefix, "_musicatk_denovo_activities.tsv")
write.table(act_df, act_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved activities     -> %s\n", act_file))

# Signature profiles: mutation_types x signatures
sig_df <- cbind(MutationType = rownames(sig_mat), as.data.frame(sig_mat))
sig_file <- paste0(out_prefix, "_musicatk_denovo_signatures.tsv")
write.table(sig_df, sig_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved signatures     -> %s\n", sig_file))

# ── Plot: exposure bar chart ──────────────────────────────────────────────────
tryCatch({
  p_exp <- plot_exposures(musica, "de_novo", modality = "SBS96", plot_type = "bar") +
    labs(
      title    = paste0("musicatk: De Novo Signature Exposures — ", cohort_id),
      subtitle = mode_label
    ) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  exp_png <- paste0(out_prefix, "_musicatk_denovo_exposures.png")
  ggsave(exp_png, p_exp, width = 16, height = 8, dpi = 200L, bg = "white")
  cat(sprintf("Saved exposures plot -> %s\n", exp_png))
}, error = function(e) {
  cat(sprintf("WARNING: exposures plot failed: %s\n", e$message))
})

# ── Plot: signature profile plots ─────────────────────────────────────────────
tryCatch({
  p_sig <- plot_signatures(musica, "de_novo") +
    labs(
      title    = paste0("musicatk: De Novo Signature Profiles — ", cohort_id),
      subtitle = mode_label
    ) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  sig_png <- paste0(out_prefix, "_musicatk_denovo_signatures.png")
  ggsave(sig_png, p_sig, width = 16, height = max(4, 4 * k_denovo),
         dpi = 200L, bg = "white")
  cat(sprintf("Saved signatures plot -> %s\n", sig_png))
}, error = function(e) {
  cat(sprintf("WARNING: signatures plot failed: %s\n", e$message))
})

cat("musicatk analysis complete.\n")
