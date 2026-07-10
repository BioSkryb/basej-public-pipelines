#!/usr/bin/env Rscript
# rscript_minimal_variant_placement.R
#
# Minimal treemut placement: read a branch-length tree + NR/NV (SNV/indel/both rows should be
# pre-filtered upstream when needed), reconcile tips vs columns, assign_to_tree only.
# No VAF discretization, di2multi, second pass, or mutation-count branch lengths.
#
# Required: --input_tree, --input_nr, --input_nv, --donor_id
# Optional: --variant_type snv|indel|both (output basename; default both), --output_dir, --tree_mut_pval
#
# Outputs:
#   <output_dir>/<donor_id>_<variant_type>_placed_variants_all.tsv — same schema as SEQUOIA_PHYLOGENY_*
#     (Chr, Pos, Ref, Alt, edge_ml, Branch, p_else_where, pass_tree_mut_pval, Patient, SampleID)
#   <output_dir>/<donor_id>_<variant_type>_placed_variants_all.rds
#   <output_dir>/<donor_id>_<variant_type>_assigned_to_branches_allunfilteredvariants.txt — high-confidence
#     subset (pass_tree_mut_pval), Chr..SampleID only (matches phylogeny assigned_to_branches pattern)
#   <output_dir>/<donor_id>_<variant_type>_tree_with_branch_length_allunfilteredvariants.tree
#   <output_dir>/<donor_id>_<variant_type>_assigned_to_branches_selectedscheme.txt — rbind of
#     full placement block + optional phylo *_placed_variants_all.tsv, provenance col

if (!require("optparse", quietly = TRUE, warn.conflicts = FALSE)) {
  install.packages("optparse", repos = "http://cran.us.r-project.org")
  library("optparse", quietly = TRUE, warn.conflicts = FALSE)
}

option_list <- list(
  make_option(c("-i", "--donor_id"), default = "Patient", type = "character",
              help = "Output file prefix [%default]"),
  make_option(c("-r", "--input_nr"), default = NULL, type = "character",
              help = "NR matrix TSV (col1 = variant id, remaining cols = samples)"),
  make_option(c("-v", "--input_nv"), default = NULL, type = "character",
              help = "NV matrix TSV (same layout as NR)"),
  make_option(c("--input_tree"), default = NULL, type = "character",
              help = "Newick tree (e.g. *_tree_with_branch_length_selectedscheme.tree)"),
  make_option(c("-o", "--output_dir"), default = ".", type = "character",
              help = "Output directory [%default]"),
  make_option(c("-t", "--tree_mut_pval"), default = 0.01, type = "numeric",
              help = "p_elsewhere threshold for reporting rows [%default]"),
  make_option(c("--variant_type"), default = "both", type = "character",
              help = "Output tag: snv | indel | both [%default]"),
  make_option(c("--phylo_placed_variants_tsv"),
              default = "", type = "character",
              help = "Optional phylo *_placed_variants_all.tsv (empty = skip in master table)")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

if (is.null(opt$input_nr))   stop("--input_nr is required")
if (is.null(opt$input_nv))   stop("--input_nv is required")
if (is.null(opt$input_tree)) stop("--input_tree is required")

patient_ID    <- opt$donor_id
output_dir    <- opt$output_dir
tree_mut_pval <- opt$tree_mut_pval
variant_type  <- tolower(opt$variant_type)
if (!variant_type %in% c("snv", "indel", "both")) {
  stop("--variant_type must be snv, indel, or both")
}

options(stringsAsFactors = FALSE)

for (pkg in c("ape", "data.table", "devtools")) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
}
if (!require("treemut", quietly = TRUE, warn.conflicts = FALSE)) {
  devtools::install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut", quietly = TRUE, warn.conflicts = FALSE)
}

run_assignment <- function(tr, NV_m, NR_m) {
  tr$edge.length <- rep(1, nrow(tr$edge))
  assign_to_tree(tr, mtr = as.matrix(NV_m), dep = as.matrix(NR_m))
}

if (nzchar(output_dir) && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat(sprintf("Reading tree: %s\n", opt$input_tree))
tsz <- file.info(opt$input_tree)$size
if (is.na(tsz) || tsz == 0) stop(sprintf("Empty or missing tree: %s", opt$input_tree))

tree <- read.tree(opt$input_tree)
if (is.null(tree) || length(tree$tip.label) < 2) {
  stop("Tree is NULL or has fewer than 2 tips")
}

if ("Ancestral" %in% tree$tip.label) {
  tree <- drop.tip(tree, "Ancestral")
  cat(sprintf("Dropped Ancestral tip — %d tips\n", length(tree$tip.label)))
}
if (length(tree$tip.label) < 2) stop("Fewer than 2 tips after dropping Ancestral")

cat("Reading NR/NV...\n")
NR_use <- data.table::fread(opt$input_nr, data.table = FALSE)
rownames(NR_use) <- NR_use[, 1]
NR_use <- NR_use[, -1, drop = FALSE]

NV_use <- data.table::fread(opt$input_nv, data.table = FALSE)
rownames(NV_use) <- NV_use[, 1]
NV_use <- NV_use[, -1, drop = FALSE]

if (!identical(rownames(NR_use), rownames(NV_use))) {
  stop("NR and NV rownames differ — align matrices before calling this script")
}
if (!identical(colnames(NR_use), colnames(NV_use))) {
  stop("NR and NV column names differ")
}

drop_smp <- names(which(colSums(NR_use) == 0))
if (length(drop_smp) > 0) {
  cat(sprintf("Dropping %d zero-coverage samples\n", length(drop_smp)))
  NR_use <- NR_use[, !colnames(NR_use) %in% drop_smp, drop = FALSE]
  NV_use <- NV_use[, !colnames(NV_use) %in% drop_smp, drop = FALSE]
}

tree_tips   <- tree$tip.label
mat_samples <- colnames(NR_use)
missing_in_tree <- setdiff(mat_samples, tree_tips)
missing_in_mat  <- setdiff(tree_tips, mat_samples)

if (length(missing_in_tree) > 0) {
  cat(sprintf("WARNING: dropping %d matrix columns not in tree tips\n", length(missing_in_tree)))
  k <- mat_samples %in% tree_tips
  NR_use <- NR_use[, k, drop = FALSE]
  NV_use <- NV_use[, k, drop = FALSE]
}
if (length(missing_in_mat) > 0) {
  cat(sprintf("WARNING: pruning %d tree tips absent from matrix\n", length(missing_in_mat)))
  tree <- drop.tip(tree, missing_in_mat)
}

if (length(tree$tip.label) < 2) stop("Fewer than 2 tips after reconcile")
if (ncol(NR_use) == 0) stop("No matrix columns left after reconcile")

NR_tree <- NR_use

cat(sprintf("Variants: %d  samples: %d  tips: %d\n",
            nrow(NR_tree), ncol(NR_tree), length(tree$tip.label)))

cat("Running assign_to_tree...\n")
res <- run_assignment(tree, NV_use, NR_tree)

if (nrow(res$summary) != nrow(NR_tree)) {
  stop(sprintf("res$summary rows (%d) != NR_tree rows (%d)",
               nrow(res$summary), nrow(NR_tree)))
}

out_prefix <- file.path(output_dir, paste0(patient_ID, "_", variant_type))

## Full placement table (same columns as rscript_sequoia_build_phylogeny_only.R placed_variants_all)
placed_variants_all <- as.data.frame(matrix(ncol = 4,
  unlist(strsplit(rownames(NR_tree), split = "_")), byrow = TRUE), stringsAsFactors = FALSE)
colnames(placed_variants_all) <- c("Chr", "Pos", "Ref", "Alt")
placed_variants_all$edge_ml <- res$summary$edge_ml
placed_variants_all$Branch <- tree$edge[res$summary$edge_ml, 2]
placed_variants_all$p_else_where <- res$summary$p_else_where
placed_variants_all$pass_tree_mut_pval <- placed_variants_all$p_else_where < tree_mut_pval
placed_variants_all$Patient <- patient_ID
placed_variants_all$SampleID <- paste(patient_ID, placed_variants_all$Branch, sep = "_")

out_placed_tsv <- paste0(out_prefix, "_placed_variants_all.tsv")
write.table(placed_variants_all, out_placed_tsv, quote = FALSE, row.names = FALSE, sep = "\t")
out_placed_rds <- paste0(out_prefix, "_placed_variants_all.rds")
saveRDS(placed_variants_all, out_placed_rds)
cat(sprintf("Wrote %s and %s (%d rows)\n", out_placed_tsv, out_placed_rds, nrow(placed_variants_all)))

## High-confidence branch list only (phylogeny Mutations_per_branch / assigned_to_branches pattern)
mpb <- placed_variants_all[placed_variants_all$pass_tree_mut_pval,
  c("Chr", "Pos", "Ref", "Alt", "Branch", "Patient", "SampleID"), drop = FALSE]

out_tbl <- paste0(out_prefix, "_assigned_to_branches_allunfilteredvariants.txt")
write.table(mpb, out_tbl, quote = FALSE, row.names = FALSE, sep = "\t")
cat(sprintf("Wrote %s (%d rows)\n", out_tbl, nrow(mpb)))

out_tree <- paste0(out_prefix, "_tree_with_branch_length_allunfilteredvariants.tree")
write.tree(tree, out_tree)
cat(sprintf("Wrote %s\n", out_tree))

# ── Master table: full unfiltered placement (phylo-compatible columns) + optional phylo table ──
assigned_dt <- data.table::as.data.table(placed_variants_all)
assigned_dt[, provenance := "unfiltered_variant_placement"]
data.table::setcolorder(assigned_dt, c("provenance", setdiff(names(assigned_dt), "provenance")))

master_parts <- list(assigned_dt)
phy_path <- opt$phylo_placed_variants_tsv
if (is.character(phy_path) && nzchar(phy_path)) {
  pfi <- file.info(phy_path)
  if (!is.na(pfi$size) && pfi$size > 0) {
    phy_dt <- data.table::fread(phy_path, data.table = TRUE)
    phy_dt[, provenance := "phylogeny_filtered_variant_placement"]
    data.table::setcolorder(phy_dt, c("provenance", setdiff(names(phy_dt), "provenance")))
    master_parts[[length(master_parts) + 1L]] <- phy_dt
    cat(sprintf("Appending phylo placement table: %d rows\n", nrow(phy_dt)))
  } else {
    cat("Phylo placed_variants TSV missing or empty — master = unfiltered placement only\n")
  }
}

master_dt <- data.table::rbindlist(master_parts, fill = TRUE)
if ("provenance" %in% names(master_dt)) {
  data.table::setcolorder(master_dt, c("provenance", setdiff(names(master_dt), "provenance")))
}
out_master <- paste0(out_prefix, "_assigned_to_branches_selectedscheme.txt")
data.table::fwrite(master_dt, out_master, sep = "\t", quote = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", out_master, nrow(master_dt)))
