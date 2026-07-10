#!/usr/bin/env Rscript
# rscript_postprocess_vaf_tree_egg.R
#
# Drop-in replacement for rscript_4.postprocess_vaf_drawtree_heatmap_intnodes.R.
# Same CLI contract (8 positional args) and same output files so the calling module
# (custom_postprocess_sequoia_drawvafheat_tree) needs no change other than the script name.
#
# Produces, per variant type, an egg::ggarrange composite (2 columns):
#   col 1 = phylogenetic tree; internal-node bootstrap support shown as DOTS,
#           colour = discretised support bin, size = support value (0-100).
#   col 2 = placed-variant heatmap; fill = VAF (NV/NR) on a fixed 0-1 scale.
#           Cells where NR == 0 (no coverage) are drawn in a distinct grey with
#           their own legend entry, so they are visually separable from VAF == 0.
#
# Args (positional, same order as the script it replaces):
#   1 file_placed  *_assigned_to_branches_selectedscheme.txt  (placement master)
#   2 file_nv      NV matrix TSV (rows = CHROM_POS_REF_ALT, cols = samples)
#   3 file_nr      NR matrix TSV (same layout)
#   4 file_tree    Newick tree with internal-node support labels
#   5 file_gt      df_all_gt_chosen.tsv (sample <tab> VariantId <tab> GT)  -> digital figure
#   6 file_vep     VEP relevance TSV (accepted for compatibility; not required)
#   7 y_text_size  heatmap y-axis (sample) text size
#   8 plot_tips    TRUE/FALSE (accepted for compatibility)
#   9 file_mand    mandatory_variants_qc_status TSV (optional). The mandatory
#                  variants that FAIL QC (PassesBoth != "Yes") are REMOVED from the
#                  placed set; the heatmaps therefore retain all non-mandatory
#                  placed variants plus the mandatory variants that pass BOTH QC
#                  filters (HQStat & pileup).
#                  Absent / empty / no-fail -> show all placed variants (unchanged).
#  10 file_bin     genotype_bin matrix TSV (optional; rows = CHROM_POS_REF_ALT,
#                  cols = samples, binarized 0/1). When supplied it DRIVES the digital
#                  genotype panel (combined with VAF/NR): bin==1 & VAF<0.8 -> 0/1,
#                  bin==1 & VAF>=0.8 -> 1/1, bin==0 & NR>0 -> 0/0, NR==0 -> no depth.
#                  Absent / sentinel -> digital panel falls back to the VCF GT table (arg 5).
#
# Outputs (cwd; the module renames them per type):
#   res_composition.pdf          VAF heatmap + tree
#   res_composition_digital.pdf  digital genotype heatmap + tree
#   res_figures.RDS              list(tree, vaf, digital plots + data frames)

suppressMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(ggtree)
  library(ape)
  library(reshape2)
  library(tidytree)
  library(egg)
  library(scales)
})

# Mandated lab theme (global instruction): ggplot via theme_ohchibi_pubr.
theme_path <- "/home/ubuntu/projects/cursor_rules/theme_ohchibi_pubr.R"
if (file.exists(theme_path)) {
  source(theme_path)
} else if (file.exists("/usr/local/bin/theme_ohchibi_pubr.R")) {
  source("/usr/local/bin/theme_ohchibi_pubr.R")
}
if (!exists("theme_ohchibi_pubr")) {
  # Fallback so the script never hard-fails if the theme file is absent in the image.
  theme_ohchibi_pubr <- function(...) ggplot2::theme_bw(base_size = 13)
}

# ─────────────────────────────────────────────────────────────────────────────
# Args
# ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: rscript_postprocess_vaf_tree_egg.R <placed> <nv> <nr> <tree> <gt> [vep] [y_text_size] [plot_tips]")
}
file_placed <- args[1]
file_nv     <- args[2]
file_nr     <- args[3]
file_tree   <- args[4]
file_gt     <- args[5]
file_vep    <- if (length(args) >= 6) args[6] else ""           # accepted, unused
y_text_size <- if (length(args) >= 7) suppressWarnings(as.numeric(args[7])) else 5
if (is.na(y_text_size)) y_text_size <- 5
# plot_tips (args[8]) accepted for compatibility; not used in this layout.
file_mand   <- if (length(args) >= 9) args[9] else ""           # mandatory_variants_qc_status TSV (optional)
file_bin    <- if (length(args) >= 10) args[10] else ""         # genotype_bin matrix (optional, drives digital GT)

# Fixed visual constants
VAF_LIMS      <- c(0, 1)
NOCOV_COLOUR  <- "cyan"            # NR == 0 cells (VAF -> NA). VAF==0 with NR>0 -> white.
NOCOV_LABEL   <- "NR = 0 (no coverage)"
SUPPORT_BINS  <- c(-Inf, 50, 70, 90, Inf)
SUPPORT_LABS  <- c("<50", "50-70", "70-90", ">=90")
SUPPORT_COLS  <- c("<50" = "#D73027", "50-70" = "#FC8D59",
                   "70-90" = "#FEE08B", ">=90" = "#1A9850")

# ─────────────────────────────────────────────────────────────────────────────
# Readers (mirror the proven logic of the script being replaced)
# ─────────────────────────────────────────────────────────────────────────────
read_placement_tsv <- function(path) {
  dt <- data.table::fread(path, sep = "\t", header = TRUE, quote = "",
                          fill = TRUE, na.strings = c("", "NA"),
                          data.table = FALSE, stringsAsFactors = FALSE)
  req <- c("Chr", "Pos", "Ref", "Alt", "Branch")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("placement TSV missing columns: ", paste(miss, collapse = ", "))
  dt <- dt[stats::complete.cases(dt[, req, drop = FALSE]), , drop = FALSE]
  for (c in req) dt[[c]] <- as.character(dt[[c]])
  dt$VariantId <- paste(dt$Chr, dt$Pos, dt$Ref, dt$Alt, sep = "_")
  dt$Branch    <- suppressWarnings(as.integer(dt$Branch))
  dt
}

# Mandatory-variant QC status reader. Returns the VariantIds that FAIL the QC
# filters (PassesBoth != "Yes") — these are removed from the placed set so the
# heatmaps retain all non-mandatory variants plus the mandatory variants that
# pass. Returns NULL when the feature is inactive (file absent / sentinel /
# header-only). A present file with no failures returns character(0) so the
# caller can distinguish "nothing to remove" from "feature off".
read_mandatory_fail <- function(path) {
  if (is.null(path) || !nzchar(path) || path == "/dev/null" || !file.exists(path)) return(NULL)
  sz <- file.info(path)$size
  if (is.na(sz) || sz < 2L) return(NULL)
  dt <- tryCatch(data.table::fread(path, sep = "\t", header = TRUE, quote = "",
                                   data.table = FALSE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)            # header-only -> feature off
  if (!all(c("VariantId", "PassesBoth") %in% names(dt))) {
    warning("mandatory QC TSV missing VariantId/PassesBoth columns; ignoring")
    return(NULL)
  }
  unique(as.character(dt$VariantId[toupper(trimws(dt$PassesBoth)) != "YES"]))
}

read_nr_nv_matrix <- function(f) {
  dt <- data.table::fread(f, header = TRUE, sep = "\t", check.names = FALSE,
                          data.table = FALSE, stringsAsFactors = FALSE)
  if (nrow(dt) == 0 || ncol(dt) == 0) return(matrix(numeric(0), 0, 0))
  if (!is.numeric(dt[[1]])) {                 # first col is the VariantId
    rn <- as.character(dt[[1]]); dt <- dt[, -1, drop = FALSE]; rownames(dt) <- rn
  }
  m <- as.matrix(dt); storage.mode(m) <- "numeric"; m
}

# ─────────────────────────────────────────────────────────────────────────────
# Load
# ─────────────────────────────────────────────────────────────────────────────
cat("Reading placement, NR/NV, tree...\n")
placed <- read_placement_tsv(file_placed)

m_nv <- read_nr_nv_matrix(file_nv)
m_nr <- read_nr_nv_matrix(file_nr)
common_v <- intersect(rownames(m_nv), rownames(m_nr))
common_s <- intersect(colnames(m_nv), colnames(m_nr))
if (!length(common_v) || !length(common_s)) stop("NR/NV matrices share no variants/samples")
m_nv <- m_nv[common_v, common_s, drop = FALSE]
m_nr <- m_nr[common_v, common_s, drop = FALSE]

# Optional binarized genotype matrix (drives the digital genotype panel when present).
m_bin <- matrix(numeric(0), 0, 0)
if (!is.null(file_bin) && nzchar(file_bin) && file_bin != "/dev/null" &&
    file.exists(file_bin) && !is.na(file.info(file_bin)$size) && file.info(file_bin)$size > 1L) {
  m_bin <- tryCatch(read_nr_nv_matrix(file_bin), error = function(e) matrix(numeric(0), 0, 0))
  cat(sprintf("genotype_bin matrix: %d variants x %d samples (drives digital genotype panel).\n",
              nrow(m_bin), ncol(m_bin)))
} else {
  cat("genotype_bin matrix: not provided — digital genotype panel falls back to VCF GT table.\n")
}

tree <- ape::read.tree(file = file_tree)
if (is.null(tree) || length(tree$tip.label) < 2) stop("tree NULL or < 2 tips")

# Long VAF table for placed variants only (variants that have a branch + are in the matrices)
df_vaf <- merge(reshape2::melt(m_nv), reshape2::melt(m_nr), by = c("Var1", "Var2"))
names(df_vaf) <- c("VariantId", "biosampleName", "NV", "NR")
df_vaf$VariantId     <- as.character(df_vaf$VariantId)
df_vaf$biosampleName <- as.character(df_vaf$biosampleName)
df_vaf <- df_vaf[df_vaf$VariantId %in% placed$VariantId, , drop = FALSE]
df_vaf$Branch <- placed$Branch[match(df_vaf$VariantId, placed$VariantId)]
df_vaf$VAF    <- ifelse(df_vaf$NR > 0, df_vaf$NV / df_vaf$NR, NA_real_)   # NR==0 -> NA
df_vaf$VAF    <- pmin(pmax(df_vaf$VAF, 0), 1)

if (!nrow(df_vaf)) stop("No placed variants overlap the NR/NV matrices — nothing to plot")

# ── Optional mandatory-variant QC restriction ────────────────────────────────
# If a mandatory_variants_qc_status TSV is supplied, DROP the mandatory variants
# that FAIL the QC filters (PassesBoth != "Yes"). The retained set is therefore
# all non-mandatory placed variants plus the mandatory variants that pass. Falls
# back to the full placed set when the file is absent/empty or lists no failing
# variant that is actually placed.
mandatory_active <- FALSE
mandatory_fail   <- read_mandatory_fail(file_mand)
if (is.null(mandatory_fail)) {
  cat("Mandatory QC filter: no usable file (absent/empty) — showing all placed variants.\n")
} else {
  to_remove <- intersect(mandatory_fail, unique(df_vaf$VariantId))
  cat(sprintf("Mandatory QC filter: %d mandatory variants fail QC; %d are placed & in NR/NV (to remove).\n",
              length(mandatory_fail), length(to_remove)))
  if (length(to_remove) > 0) {
    df_vaf <- df_vaf[!(df_vaf$VariantId %in% to_remove), , drop = FALSE]
    mandatory_active <- TRUE
    cat(sprintf("Removed %d failing mandatory variants; %d placed variants retained.\n",
                length(to_remove), length(unique(df_vaf$VariantId))))
    if (!nrow(df_vaf)) stop("All placed variants removed by mandatory QC filter — nothing to plot")
  } else {
    cat("No failing mandatory variant is placed — showing all placed variants.\n")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Tree panel: internal-node support as colour-binned, size-scaled dots.
# Equal branch lengths (cladogram) + aligned tips so all tips reach the SAME right
# edge, matching the original script's geom_tiplab(align=TRUE) look.
# ─────────────────────────────────────────────────────────────────────────────
# Facet (internal-node block) order = ladderized internal-node edge order (as in original).
ntip <- length(tree$tip.label)
tree_lad <- ape::ladderize(tree, right = TRUE)
order_nodes <- tree_lad$edge[tree_lad$edge[, 2] > ntip, 2]

# Equal branch lengths (cladogram) + ladderize=FALSE so the tip (row) order matches the
# original script; aligned tips reach the same right edge.
tree$edge.length <- rep(1, nrow(tree$edge))
p_tree_base <- ggtree(tree, ladderize = FALSE, size = 0.4, colour = "grey30")
td <- p_tree_base$data

internal <- td[!td$isTip, , drop = FALSE]
internal$support <- suppressWarnings(as.numeric(internal$label))
internal <- internal[!is.na(internal$support), , drop = FALSE]
internal$support_bin <- cut(internal$support, breaks = SUPPORT_BINS,
                            labels = SUPPORT_LABS, right = FALSE)

p_tree <- p_tree_base +
  geom_tiplab(size = 0, align = TRUE, linetype = "dotted", linesize = 0.2) +
  geom_point(data = internal,
             aes(x = x, y = y, fill = support_bin, size = support),
             shape = 21, colour = "grey20", stroke = 0.25) +
  scale_fill_manual(values = SUPPORT_COLS, drop = FALSE,
                    name = "Bootstrap support (bin)") +
  scale_size_continuous(range = c(1, 6), limits = c(0, 100),
                        breaks = c(25, 50, 75, 100), name = "Support value") +
  guides(fill = guide_legend(override.aes = list(size = 8), order = 1),
         size = guide_legend(order = 2)) +
  coord_cartesian(clip = "off") +
  ggtitle("Phylogeny (node = bootstrap support)") +
  theme_void() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size = 18), legend.text = element_text(size = 18),
        legend.key.size = grid::unit(1.0, "cm"),
        plot.title = element_text(hjust = 0.5, size = 12))

# Tip order (bottom -> top) so the heatmap rows align with tree tips under egg
tip_order <- td %>% dplyr::filter(isTip) %>% dplyr::arrange(y) %>% dplyr::pull(label)
tip_order <- intersect(tip_order, unique(df_vaf$biosampleName))

# ─────────────────────────────────────────────────────────────────────────────
# Heatmap helpers (rows = samples aligned to tips; cols = variants ordered by branch)
# ─────────────────────────────────────────────────────────────────────────────
# Original-style faceted blocks: keep variants placed on INTERNAL nodes only, one
# facet panel per node (per-variant columns inside, width proportional to #variants),
# facets ordered by the ladderized internal-node order, rows = samples (tip order).
df_int <- df_vaf[df_vaf$Branch %in% order_nodes, , drop = FALSE]
if (!nrow(df_int)) stop("No variants placed on internal nodes — nothing to plot")

# Annotate the tree with internal-node IDs, but ONLY for nodes that have a facet
# in the variant heatmap (i.e. nodes with >=1 placed variant). The label is the
# bare node number so it matches the facet header "<node> (n=...)"; drawn just
# left of the support dot so the two don't overlap. Plain text, sized to match the
# heatmap sample labels (element_text pt -> geom_text mm via ggplot2::.pt).
facet_node_ids <- unique(order_nodes[order_nodes %in% df_int$Branch])
node_id_df     <- td[td$node %in% facet_node_ids & !td$isTip, , drop = FALSE]
if (nrow(node_id_df)) {
  p_tree <- p_tree +
    geom_text(data = node_id_df, aes(x = x, y = y, label = node),
              inherit.aes = FALSE, size = (y_text_size * 2) / ggplot2::.pt,
              colour = "grey15", hjust = 1.3, vjust = -0.4)
}

n_per_node <- df_int %>% dplyr::distinct(VariantId, Branch) %>% dplyr::count(Branch, name = "n_var")
node_label <- function(b) sprintf("%s (n=%d)", b,
                                  n_per_node$n_var[match(b, n_per_node$Branch)])
facet_levels <- node_label(order_nodes[order_nodes %in% df_int$Branch])

# Within-facet variant order: by node order then VariantId
var_levels <- df_int %>%
  dplyr::arrange(match(Branch, order_nodes), VariantId) %>%
  dplyr::pull(VariantId) %>% unique()

prep_facet_df <- function(df) {
  df$biosampleName <- factor(df$biosampleName, levels = tip_order)
  df$LabelFacet    <- factor(node_label(df$Branch), levels = facet_levels)
  df$VariantId     <- factor(df$VariantId, levels = var_levels)
  df[!is.na(df$biosampleName) & !is.na(df$LabelFacet), , drop = FALSE]
}

# Invisible layer that injects an "NR = 0" key into a colour legend (no ggnewscale needed).
# Uses REAL first-level facet/x/y so it never introduces NA categories.
nocov_legend_layer <- function() {
  dummy <- data.frame(grp = NOCOV_LABEL, stringsAsFactors = FALSE)
  dummy$VariantId    <- factor(var_levels[1],   levels = var_levels)
  dummy$biosampleName<- factor(tip_order[1],    levels = tip_order)
  dummy$LabelFacet   <- factor(facet_levels[1], levels = facet_levels)
  list(
    geom_point(data = dummy, aes(x = VariantId, y = biosampleName, colour = grp),
               inherit.aes = FALSE, alpha = 0, na.rm = TRUE),
    scale_colour_manual(values = setNames(NOCOV_COLOUR, NOCOV_LABEL), name = NULL),
    guides(colour = guide_legend(order = 3,
             override.aes = list(alpha = 1, size = 5, shape = 15),
             label.theme = element_text(angle = 0, size = 18)))   # keep NR=0 label horizontal
  )
}

facet_heatmap_theme <- function() {
  # Heatmap labels (sample names + facet/node strips) doubled per request.
  theme_ohchibi_pubr() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = y_text_size * 2, hjust = 1),
          panel.spacing = grid::unit(1, "pt"),
          panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
          strip.text.x = element_text(size = 12, angle = 90, hjust = 0, vjust = 0.5),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18, angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")
}

build_vaf_heatmap <- function() {
  df <- prep_facet_df(df_int)
  # VAF==0 (covered, no alt) stays 0 (viridis low); NR==0 is NA -> NOCOV_COLOUR (distinct).
  ggplot(df, aes(x = VariantId, y = biosampleName)) +
    geom_raster(aes(fill = VAF)) +
    facet_grid(. ~ LabelFacet, scales = "free_x") +
    # VAF==0 (NR>0) -> white (first colour); VAF>0 -> ocean.amp ramp; NR==0 -> NA -> cyan.
    scale_fill_gradientn(colours = c("white", paletteer::paletteer_c("pals::ocean.amp", 255)),
                         limits = VAF_LIMS, oob = scales::squish,
                         na.value = NOCOV_COLOUR, name = "VAF (NV/NR)") +
    nocov_legend_layer() +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    facet_heatmap_theme() +
    labs(x = "Variants placed per internal node (facet = node (n))", y = NULL) +
    ggtitle(paste0("Placed-variant VAF, faceted by internal node",
                   if (mandatory_active) " (failing mandatory variants removed)" else ""))
}

# ─────────────────────────────────────────────────────────────────────────────
# Digital genotype heatmap (fill = genotype 0/0.5/1 from the GT table)
# ─────────────────────────────────────────────────────────────────────────────
read_gt_numeric <- function(path) {
  empty <- data.frame(biosampleName = character(0), VariantId = character(0), geno = character(0))
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(empty)
  if (is.na(file.info(path)$size) || file.info(path)$size < 2L) return(empty)
  dt <- tryCatch(data.table::fread(path, sep = "\t", header = FALSE, quote = "",
                                   data.table = FALSE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(dt) || ncol(dt) < 3) return(empty)
  dt <- dt[, 1:3]; names(dt) <- c("biosampleName", "VariantId", "GT")
  alleles <- function(g) {
    a <- unlist(strsplit(gsub("\\|", "/", as.character(g)), "/"))
    a <- suppressWarnings(as.integer(a)); a <- a[!is.na(a)]
    if (!length(a)) return(NA_character_)
    if (all(a == 0)) "0" else if (all(a > 0)) "1" else "0.5"
  }
  dt$geno <- vapply(dt$GT, alleles, character(1))
  dt[, c("biosampleName", "VariantId", "geno")]
}

# Discrete genotype classes + palette (same 4 colours as the phylogeny_from_annotated_vcf
# postprocess: rscript_4.postprocess_vaf_drawtree_heatmap_intnodes.R).
GENO_LEVELS <- c("0/0", "0/1", "1/1", "No depth")
GENO_COLS   <- c("0/0" = "#1082A2", "0/1" = "#A0CC2C", "1/1" = "#12284C", "No depth" = "#DD14D3")

# Genotype CLASS per variant x sample, one of GENO_LEVELS.
# When the genotype_bin matrix is supplied it is authoritative and combined with
# VAF/NR per the rules:
#   bin==1 & VAF <  0.8        -> 0/1
#   bin==1 & VAF >= 0.8        -> 1/1
#   bin==0 (or absent) & NR>0  -> 0/0
#   NR == 0 (no read depth)    -> No depth   (overrides the above)
geno_class_from_bin <- function(df) {
  rk <- as.character(df$VariantId); ck <- as.character(df$biosampleName)
  bin <- rep(NA_real_, nrow(df))
  ok  <- rk %in% rownames(m_bin) & ck %in% colnames(m_bin)
  if (any(ok)) bin[ok] <- m_bin[cbind(rk[ok], ck[ok])]
  vaf <- df$VAF; nr <- df$NR
  present <- !is.na(bin) & bin >= 1
  g <- rep("0/0", nrow(df))                              # bin 0 / missing -> 0/0
  g[present & !is.na(vaf) & vaf >= 0.8] <- "1/1"
  g[present & (is.na(vaf) | vaf < 0.8)] <- "0/1"
  g[!is.na(nr) & nr == 0] <- "No depth"                 # no depth wins
  factor(g, levels = GENO_LEVELS)
}

# Fallback: derive the same classes from the VCF GT table (0/0, 0/1, 1/1; NR==0 -> No depth).
geno_class_from_gt <- function(df, file_gt) {
  gt <- read_gt_numeric(file_gt)
  geno <- gt$geno[match(paste(df$biosampleName, df$VariantId),
                        paste(gt$biosampleName, gt$VariantId))]
  g <- c("0" = "0/0", "0.5" = "0/1", "1" = "1/1")[as.character(geno)]
  g[is.na(g)] <- "0/0"                                   # uncalled -> ref
  g[!is.na(df$NR) & df$NR == 0] <- "No depth"            # no coverage
  factor(unname(g), levels = GENO_LEVELS)
}

# Digital genotype, original-style faceted per internal node: per-variant columns,
# fill = discrete genotype class with the fixed 4-colour palette (0/0, 0/1, 1/1, No depth).
build_digital_heatmap <- function(file_gt) {
  df <- prep_facet_df(df_int)
  if (nrow(m_bin)) {
    df$geno_cls <- geno_class_from_bin(df)
    cat("Digital genotype panel: derived from genotype_bin + VAF/NR (discrete classes).\n")
  } else {
    df$geno_cls <- geno_class_from_gt(df, file_gt)
    cat("Digital genotype panel: derived from VCF GT table (discrete classes; no genotype_bin).\n")
  }
  ggplot(df, aes(x = VariantId, y = biosampleName)) +
    geom_raster(aes(fill = geno_cls)) +
    facet_grid(. ~ LabelFacet, scales = "free_x") +
    scale_fill_manual(values = GENO_COLS, name = "Genotype", drop = FALSE, na.value = "#DDDDDD") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    facet_heatmap_theme() +
    labs(x = "Variants placed per internal node (facet = node (n))", y = NULL) +
    ggtitle(paste0("Placed-variant digital genotype, faceted by internal node",
                   if (mandatory_active) " (failing mandatory variants removed)" else ""))
}

# ─────────────────────────────────────────────────────────────────────────────
# Assemble + save
# ─────────────────────────────────────────────────────────────────────────────
assemble <- function(p_tree, p_heat) {
  egg::ggarrange(p_tree, p_heat, nrow = 1, ncol = 2,
                 widths = c(0.3, 1), draw = FALSE)
}

save_pdf <- function(g, file, width = 22, height = 18) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (inherits(g, c("gtable", "egg", "grob", "gTree"))) grid::grid.draw(g) else print(g)
}

cat("Building VAF figure...\n")
p_vaf  <- build_vaf_heatmap()
comp_vaf <- assemble(p_tree, p_vaf)
save_pdf(comp_vaf, "res_composition.pdf")

cat("Building digital genotype figure...\n")
p_dig  <- build_digital_heatmap(file_gt)
comp_dig <- assemble(p_tree, p_dig)
save_pdf(comp_dig, "res_composition_digital.pdf")

saveRDS(list(
  p_tree = p_tree, p_vaf = p_vaf, p_digital = p_dig,
  composition_vaf = comp_vaf, composition_digital = comp_dig,
  df_int = df_int, tip_order = tip_order, facet_levels = facet_levels
), file = "res_figures.RDS")

cat(sprintf("Done: %d variants on %d internal-node facets x %d samples; %d nodes with support.\n",
            length(var_levels), length(facet_levels), length(tip_order), nrow(internal)))
