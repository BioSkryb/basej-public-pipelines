#!/usr/bin/env Rscript
# rscript_variant_placement.R
#
# Place variants from NR/NV matrices onto a pre-built phylogenetic tree using
# treemut::assign_to_tree().  Designed to work with any NR/NV pair (e.g. the
# full unfiltered matrices) against a tree topology established from a filtered
# set (e.g. the pileup phylogeny).
#
# Key differences from sequoia_build_phylogeny.R:
#   - No germline / beta-binomial / depth filtering  (variants are pre-filtered)
#   - No MPBoot tree building                        (tree topology is an input)
#   - No mixture-model sample selection
#   - Variant-type subsetting (snv | indel | both)   handled internally
#   - Gender is supplied explicitly via --gender     (not inferred from depth)
#
# Outputs (in --output_dir):
#   <donor_id>_<variant_type>_assigned_to_tree.Rdata
#   <donor_id>_<variant_type>_tree_with_branch_length_selectedscheme.tree
#   <donor_id>_<variant_type>_assigned_to_branches.txt

# ─── optparse ────────────────────────────────────────────────────────────────
if (!require("optparse", quietly = TRUE, warn.conflicts = FALSE)) {
  install.packages("optparse", repos = "http://cran.us.r-project.org")
  library("optparse", quietly = TRUE, warn.conflicts = FALSE)
}

option_list <- list(
  make_option(c("-i", "--donor_id"),
              default = "Patient", type = "character",
              help    = "Donor / patient ID used as output file prefix [%default]"),
  make_option(c("-r", "--input_nr"),
              default = NULL, type = "character",
              help    = "NR matrix TSV (rows=variants, cols=samples)"),
  make_option(c("-v", "--input_nv"),
              default = NULL, type = "character",
              help    = "NV matrix TSV (rows=variants, cols=samples)"),
  make_option(c("--input_tree"),
              default = NULL, type = "character",
              help    = "Path to pre-built tree (.tree or .treefile); branch lengths are overwritten"),
  make_option(c("-o", "--output_dir"),
              default = "", type = "character",
              help    = "Output directory [current dir]"),
  make_option(c("--variant_type"),
              default = "both", type = "character",
              help    = "Variant type to place: snv | indel | both [%default]"),
  make_option(c("--gender"),
              default = "female", type = "character",
              help    = "Sample gender: male | female — controls VAF thresholds on chrX/Y [%default]"),
  make_option(c("-t", "--tree_mut_pval"),
              default = 0.01, type = "numeric",
              help    = "p_elsewhere threshold for treemut mutation assignment [%default]"),
  make_option(c("--keep_ancestral"),
              default = FALSE, type = "logical",
              help    = "Retain ancestral branch during assignment [%default]"),
  make_option(c("--create_multi_tree"),
              default = TRUE, type = "logical",
              help    = "Collapse zero-length internal branches to polytomy [%default]"),
  make_option(c("--vaf_absent"),
              default = 0.1, type = "numeric",
              help    = "VAF below which a variant is called absent (autosomal) [%default]"),
  make_option(c("--vaf_present"),
              default = 0.3, type = "numeric",
              help    = "VAF above which a variant is called present (autosomal) [%default]"),
  make_option(c("--genotype_conv_prob"),
              default = FALSE, type = "logical",
              help    = "Use probabilistic genotype matrix instead of hard VAF thresholds [%default]"),
  make_option(c("--min_pval_for_true_somatic"),
              default = 0.05, type = "numeric",
              help    = "Binomial p-value threshold for somatic presence (genotype_conv_prob mode) [%default]"),
  make_option(c("--min_variant_reads_shared"),
              default = 2, type = "numeric",
              help    = "Minimum variant reads for somatic call (genotype_conv_prob mode) [%default]"),
  make_option(c("--min_vaf_shared"),
              default = 2, type = "numeric",
              help    = "Minimum VAF for somatic call (genotype_conv_prob mode) [%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# ─── Validate required inputs ────────────────────────────────────────────────
if (is.null(opt$input_nr))   stop("--input_nr is required")
if (is.null(opt$input_nv))   stop("--input_nv is required")
if (is.null(opt$input_tree)) stop("--input_tree is required")

gender <- tolower(opt$gender)
if (!gender %in% c("male", "female")) {
  stop(sprintf("--gender must be 'male' or 'female', got: '%s'", gender))
}

# ─── Load packages ───────────────────────────────────────────────────────────
options(stringsAsFactors = FALSE)

cran_pkgs <- c("ape", "data.table", "devtools")
for (pkg in cran_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
}

if (!require("treemut", quietly = TRUE, warn.conflicts = FALSE)) {
  devtools::install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut", quietly = TRUE, warn.conflicts = FALSE)
}

# ─── Helper functions ─────────────────────────────────────────────────────────

add_ancestral_outgroup <- function(tree, outgroup_name = "Ancestral") {
  tmp           <- tree$edge
  N             <- length(tree$tip.label)
  newroot       <- N + 2
  renamedroot   <- N + 3
  ancestral_tip <- N + 1
  tmp <- ifelse(tmp > N, tmp + 2, tmp)
  tree$edge        <- rbind(c(newroot, renamedroot), tmp, c(newroot, ancestral_tip))
  tree$edge.length <- c(0, tree$edge.length, 0)
  tree$tip.label   <- c(tree$tip.label, outgroup_name)
  tree$Nnode       <- tree$Nnode + 1L
  mode(tree$Nnode) <- "integer"
  mode(tree$edge)  <- "integer"
  tree
}

binom_pval_matrix <- function(NV, NR, gender) {
  NR_nz    <- NR; NR_nz[NR_nz == 0] <- 1
  pval_mat <- matrix(0, nrow = nrow(NV), ncol = ncol(NV),
                     dimnames = list(rownames(NV), colnames(NV)))
  for (i in seq_len(nrow(NV))) {
    for (j in seq_len(ncol(NV))) {
      p <- if (gender == "male" && grepl("X|Y", rownames(NV)[i])) 0.95 else 0.5
      pval_mat[i, j] <- binom.test(NV[i, j], NR_nz[i, j],
                                   p = p, alternative = "less")$p.value
    }
  }
  pval_mat
}

edge_lengths_from_res <- function(tree, res, pval_thresh) {
  el_nz            <- table(res$summary$edge_ml[res$summary$p_else_where < pval_thresh])
  el               <- rep(0, nrow(tree$edge))
  names(el)        <- seq_len(nrow(tree$edge))
  el[names(el_nz)] <- el_nz
  tree$edge.length <- as.numeric(el)
  tree
}

# ─── Assign parameters ───────────────────────────────────────────────────────
output_dir         <- opt$output_dir
patient_ID         <- opt$donor_id
variant_type       <- tolower(opt$variant_type)
tree_mut_pval      <- opt$tree_mut_pval
keep_ancestral     <- opt$keep_ancestral
create_multi_tree  <- opt$create_multi_tree
VAF_absent         <- opt$vaf_absent
VAF_present        <- opt$vaf_present
genotype_conv_prob <- opt$genotype_conv_prob
min_pval_somatic   <- opt$min_pval_for_true_somatic
min_var_reads      <- opt$min_variant_reads_shared
min_vaf            <- opt$min_vaf_shared

cat(sprintf("Gender: %s\n", gender))

if (output_dir != "" && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ─── Validate and read tree (early, so we fail fast on bad/empty input) ───────
cat(sprintf("Reading tree: %s\n", opt$input_tree))

tree_size <- file.info(opt$input_tree)$size
if (is.na(tree_size) || tree_size == 0) {
  stop(sprintf("Tree file is empty or does not exist: %s", opt$input_tree))
}

tree <- tryCatch(
  read.tree(opt$input_tree),
  error   = function(e) stop(sprintf("Failed to parse tree file '%s': %s",
                                     opt$input_tree, e$message)),
  warning = function(w) {
    cat(sprintf("WARNING while reading tree: %s\n", w$message))
    suppressWarnings(read.tree(opt$input_tree))
  }
)

if (is.null(tree)) {
  stop(sprintf("read.tree() returned NULL for file: %s", opt$input_tree))
}
if (length(tree$tip.label) < 2) {
  stop(sprintf("Tree has fewer than 2 tips (%d) — cannot perform assignment: %s",
               length(tree$tip.label), opt$input_tree))
}
cat(sprintf("Tree has %d tips before any pruning\n", length(tree$tip.label)))

# Drop Ancestral tip if present (MPBoot always adds one)
if ("Ancestral" %in% tree$tip.label) {
  tree <- drop.tip(tree, "Ancestral")
  cat(sprintf("Dropped 'Ancestral' tip — tree now has %d tips\n", length(tree$tip.label)))
}

if (length(tree$tip.label) < 2) {
  stop("Tree has fewer than 2 tips after dropping 'Ancestral' — cannot perform assignment")
}

# ─── Read NR / NV matrices ───────────────────────────────────────────────────
cat("Reading NR/NV matrices...\n")
NR <- fread(opt$input_nr, data.table = FALSE)
rownames(NR) <- NR[, 1]; NR <- NR[, -1, drop = FALSE]

NV <- fread(opt$input_nv, data.table = FALSE)
rownames(NV) <- NV[, 1]; NV <- NV[, -1, drop = FALSE]

# Drop zero-coverage samples (consistent with sequoia_build_phylogeny.R)
drop_smp <- names(which(colSums(NR) == 0))
if (length(drop_smp) > 0) {
  cat(sprintf("Dropping %d zero-coverage samples: %s\n",
              length(drop_smp), paste(drop_smp, collapse = ", ")))
  NR <- NR[, !colnames(NR) %in% drop_smp, drop = FALSE]
  NV <- NV[, !colnames(NV) %in% drop_smp, drop = FALSE]
}

Muts <- rownames(NV)
cat(sprintf("Variants: %d   Samples: %d\n", length(Muts), ncol(NV)))

# ─── Variant-type row subsetting ─────────────────────────────────────────────
Muts_coord <- matrix(ncol = 4,
                     unlist(strsplit(Muts, split = "_")),
                     byrow = TRUE)
is_indel <- nchar(Muts_coord[, 3]) > 1 | nchar(Muts_coord[, 4]) > 1

keep <- switch(variant_type,
  snv   = !is_indel,
  indel = is_indel,
  both  = rep(TRUE, length(Muts)),
  stop(sprintf("Unknown --variant_type '%s'. Use snv, indel, or both.", variant_type))
)

if (sum(keep) == 0) {
  stop(sprintf("No %s variants found in matrix — aborting.", variant_type))
}

NR   <- NR[keep, , drop = FALSE]
NV   <- NV[keep, , drop = FALSE]
Muts <- rownames(NV)
cat(sprintf("After %s subsetting: %d variants\n", variant_type, length(Muts)))

# ─── Genotype discretization → present_vars_full ─────────────────────────────
XY_chromosomal <- grepl("X|Y", Muts)
autosomal      <- !XY_chromosomal
NR_nz          <- NR; NR_nz[NR_nz == 0] <- 1

if (genotype_conv_prob) {
  cat("Building probabilistic genotype matrix...\n")
  pval_mat          <- binom_pval_matrix(NV, NR, gender)
  min_var_reads_mat <- NV >= min_var_reads
  min_pval_mat      <- pval_mat > min_pval_somatic
  min_vaf_mat       <- NV / NR_nz > min_vaf
  genotype_bin      <- min_var_reads_mat * min_pval_mat * min_vaf_mat
  genotype_bin[NV > 0 & pval_mat > 0.01 & genotype_bin != 1]   <- 0.5
  genotype_bin[NV >= 3 & pval_mat > 0.001 & genotype_bin != 1] <- 0.5
  genotype_bin[(NV == 0) & (pval_mat > 0.05)]                   <- 0.5
} else {
  genotype_bin <- as.matrix(NV / NR_nz)
  if (gender == "male") {
    genotype_bin[autosomal, ][genotype_bin[autosomal, ] < VAF_absent]              <- 0
    genotype_bin[autosomal, ][genotype_bin[autosomal, ] >= VAF_present]            <- 1
    genotype_bin[XY_chromosomal, ][genotype_bin[XY_chromosomal, ] < (2*VAF_absent)]   <- 0
    genotype_bin[XY_chromosomal, ][genotype_bin[XY_chromosomal, ] >= (2*VAF_present)] <- 1
  } else {
    genotype_bin[genotype_bin < VAF_absent]  <- 0
    genotype_bin[genotype_bin >= VAF_present] <- 1
  }
  genotype_bin[genotype_bin > 0 & genotype_bin < 1] <- 0.5
}

present_vars_full <- rowSums(genotype_bin > 0) > 0
cat(sprintf("Variants present in >=1 sample: %d / %d\n",
            sum(present_vars_full), length(present_vars_full)))

if (sum(present_vars_full) == 0) {
  stop("No variants are present in any sample after genotype discretization — cannot assign")
}

NR_tree <- NR[present_vars_full, , drop = FALSE]
NV_tree <- NV[present_vars_full, , drop = FALSE]

# ─── Reconcile sample sets: tree tips vs matrix columns ──────────────────────
tree_tips   <- tree$tip.label
mat_samples <- colnames(NR_tree)

missing_in_tree <- setdiff(mat_samples, tree_tips)
missing_in_mat  <- setdiff(tree_tips, mat_samples)

if (length(missing_in_tree) > 0) {
  cat(sprintf("WARNING: %d matrix samples absent from tree — dropping from matrices: %s\n",
              length(missing_in_tree), paste(missing_in_tree, collapse = ", ")))
  keep_cols <- mat_samples %in% tree_tips
  NR_tree   <- NR_tree[, keep_cols, drop = FALSE]
  NV_tree   <- NV_tree[, keep_cols, drop = FALSE]
}
if (length(missing_in_mat) > 0) {
  cat(sprintf("WARNING: %d tree tips absent from matrix — pruning from tree: %s\n",
              length(missing_in_mat), paste(missing_in_mat, collapse = ", ")))
  tree <- drop.tip(tree, missing_in_mat)
}

if (length(tree$tip.label) < 2) {
  stop("Fewer than 2 tree tips remain after reconciling with the matrix — cannot assign")
}
if (ncol(NR_tree) == 0) {
  stop("No matrix samples remain after reconciling with the tree — cannot assign")
}

cat(sprintf("Tree tips: %d   Matrix samples: %d\n",
            length(tree$tip.label), ncol(NR_tree)))

# ─── assign_to_tree helpers ───────────────────────────────────────────────────
run_assignment <- function(tr, NV_m, NR_m, p_err = NULL) {
  tr$edge.length <- rep(1, nrow(tr$edge))
  if (!is.null(p_err)) {
    assign_to_tree(tr, mtr = as.matrix(NV_m), dep = as.matrix(NR_m),
                   error_rate = p_err)
  } else {
    assign_to_tree(tr, mtr = as.matrix(NV_m), dep = as.matrix(NR_m))
  }
}

# ─── First assignment pass ────────────────────────────────────────────────────
cat("Running assign_to_tree...\n")

if (keep_ancestral) {
  cat("Adding ancestral outgroup for assignment...\n")
  tree   <- add_ancestral_outgroup(tree)
  NR_anc <- NR_tree; NR_anc$Ancestral <- 30
  NV_anc <- NV_tree; NV_anc$Ancestral <- 0
  p_err  <- rep(0.01, ncol(NV_anc))
  p_err[colnames(NV_anc) == "Ancestral"] <- 1e-6
  res    <- run_assignment(tree, NV_anc, NR_anc, p_err)
  NV_use <- NV_anc
  NR_use <- NR_anc
} else {
  res    <- run_assignment(tree, NV_tree, NR_tree)
  NV_use <- NV_tree
  NR_use <- NR_tree
}

tree <- edge_lengths_from_res(tree, res, tree_mut_pval)

# ─── Optional: polytomy conversion + re-assignment ───────────────────────────
if (create_multi_tree) {
  cat("Converting dichotomous tree to multifurcating (di2multi)...\n")
  if (keep_ancestral) {
    ROOT        <- tree$edge[1, 1]
    anc_tip_idx <- which(tree$tip.label == "Ancestral")
    mask        <- tree$edge[, 1] == ROOT & tree$edge[, 2] != anc_tip_idx
    cur         <- tree$edge.length[mask]
    tree$edge.length[mask] <- ifelse(cur == 0, 1, cur)
  }
  tree <- di2multi(tree)
  res  <- run_assignment(tree, NV_use, NR_use)
  tree <- edge_lengths_from_res(tree, res, tree_mut_pval)
}

# ─── Save outputs ────────────────────────────────────────────────────────────
cat("Saving outputs...\n")
out_prefix <- file.path(output_dir, paste0(patient_ID, "_", variant_type))

saveRDS(res, paste0(out_prefix, "_assigned_to_tree.Rdata"))
write.tree(tree, paste0(out_prefix, "_tree_with_branch_length_selectedscheme.tree"))

# Per-branch mutation table (same structure as sequoia_build_phylogeny.R)
mpb <- as.data.frame(matrix(ncol = 4,
         unlist(strsplit(rownames(NR_tree), split = "_")),
         byrow = TRUE))
colnames(mpb) <- c("Chr", "Pos", "Ref", "Alt")
mpb$Branch    <- tree$edge[res$summary$edge_ml, 2]
mpb           <- mpb[res$summary$p_else_where < tree_mut_pval, ]
mpb$Patient   <- patient_ID
mpb$SampleID  <- paste(patient_ID, mpb$Branch, sep = "_")
write.table(mpb, paste0(out_prefix, "_assigned_to_branches.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# ─── Summary ─────────────────────────────────────────────────────────────────
cat(sprintf("\n=== Variant placement summary: %s [%s] ===\n", patient_ID, variant_type))
cat(sprintf("  Gender              : %s\n", gender))
cat(sprintf("  Input variants      : %d\n", nrow(NR)))
cat(sprintf("  Present (>=1 sample): %d\n", sum(present_vars_full)))
cat(sprintf("  Assigned (p<%g)    : %d\n",  tree_mut_pval,
            sum(res$summary$p_else_where < tree_mut_pval)))
cat(sprintf("  Tree tips           : %d\n", length(tree$tip.label)))
cat("Done.\n")
