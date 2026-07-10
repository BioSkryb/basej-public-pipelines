#!/usr/bin/env Rscript
## Build master variant × filter-status table.
## Called from CUSTOM_VARIANT_FILTER_PROVENANCE. Usage: Rscript master_table.R <group_id>
## Expects inputs in current directory (Nextflow work dir); writes outputs there.
##
## Pileup aggregation + SEQUOIA columns are now precomputed by variant_prefilter_table.R
## (module VARIANT_PREFILTER_TABLE, which runs before VAF_SPLIT_VARIANTS_HEXBIN); this
## script consumes that prefilter table instead of recomputing them. It also folds in
## VAF_SPLIT_VARIANTS_HEXBIN membership as VAFSplit_* columns so the master table — and
## everything provenance derives from it — reflects the VAF-split cascade.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript master_table.R <group_id>")
group_id <- args[1]

suppressPackageStartupMessages(library(data.table))

bulk_file      <- list.files(".", pattern = "bulk_filter_provenance_.*\\.tsv$")[1]
binom_file     <- list.files(".", pattern = "_binomial_betabinomial\\.tsv$")[1]
vep_file       <- list.files(".", pattern = "vep_filter_provenance_.*\\.tsv$")[1]
prefilter_file <- paste0("variant_prefilter_table_", group_id, ".tsv")
if (!file.exists(prefilter_file)) stop("Prefilter table not found: ", prefilter_file)

cat("Files found:\n")
cat("  Bulk:      ", bulk_file,      "\n")
cat("  Binom:     ", binom_file,     "\n")
cat("  VEP:       ", vep_file,       "\n")
cat("  Prefilter: ", prefilter_file, "\n\n")

# Read TSVs with fread (fast, low-memory) and keep the original parsing contract:
# tab-delimited, "NA" -> NA, quotes treated literally (matches the old read.delim
# usage; downstream vcf_annotation_table.R also reads the master with quote="").
read_tsv <- function(f) fread(f, sep = "\t", header = TRUE, na.strings = "NA",
                              quote = "", data.table = TRUE, showProgress = FALSE)

# ── 1. Bulk provenance ─────────────────────────────────────────────────────────
bulk <- read_tsv(bulk_file)
setnames(bulk, 1, "VariantId")
setnames(bulk, 2:ncol(bulk), paste0("Bulk_", names(bulk)[-1]))

# ── 2. Binomial / betabinomial (first column already named VariantId) ──────────
binom <- read_tsv(binom_file)
setnames(binom, 2:ncol(binom), paste0("Binom_", names(binom)[-1]))

# ── 3. VEP provenance ──────────────────────────────────────────────────────────
vep <- read_tsv(vep_file)
vep[, VariantId := paste(CHROM, POS, REF, ALT, sep = "_")]
orig_vep_cols <- setdiff(names(vep), "VariantId")
setnames(vep, orig_vep_cols, paste0("VEP_", orig_vep_cols))
vep_sel <- vep[, c("VariantId", paste0("VEP_", orig_vep_cols)), with = FALSE]

# ── 4. Prefilter table (Pileup_* aggregation + SecondRunSequoia_* columns) ─────
prefilter <- read_tsv(prefilter_file)

# ── 5. Full outer join: bulk → binom → vep → prefilter (keyed data.table merge) ─
# Keyed merges are far faster / lower-memory than base merge at millions of rows.
# Column order (VariantId, Bulk_*, Binom_*, VEP_*, Pileup_*, SecondRunSequoia_*) matches
# the previous base-R Reduce(merge); the final base order() reproduces the prior row order.
master <- Reduce(function(a, b) merge(a, b, by = "VariantId", all = TRUE),
                 list(bulk, binom, vep_sel, prefilter))
master <- as.data.frame(master, stringsAsFactors = FALSE)
master <- master[order(master$VariantId), ]

# ── 6. Fold in VAF_SPLIT_VARIANTS_HEXBIN membership (VAFSplit_* columns) ────────
# Reads the VAF-split outputs (staged in the work dir) and adds, per variant:
#   VAFSplit_Class         singleton | shared | none
#   VAFSplit_SharedRetained Yes (shared & passed Rho/qval) | No (shared, not retained) | . (not shared)
#   VAFSplit_HQStat        Yes/No  (HQ statistical set = retained-shared ∪ singletons)
#   VAFSplit_DepthPass     Yes/No  (survived the depth filter)
#   VAFSplit_QCPass        Yes/No  (Pileup_Verdict==Pass cascade stage)
#   VAFSplit_ForPhylogeny  Yes/No  (>=70%-cell coverage)
read_ids  <- function(f) if (file.exists(f)) unique(readLines(f)) else character(0)
read_bvid <- function(f) if (file.exists(f))
  unique(fread(f, sep = "\t", header = TRUE, select = 1L, data.table = FALSE)[[1]]) else character(0)

g <- group_id
singleton_ids <- read_ids(paste0(g, "_singleton_variants.txt"))
shared_ids    <- read_ids(paste0(g, "_shared_variants.txt"))
retained_ids  <- read_ids(paste0(g, "_shared_variants_retained.txt"))
hqstat_ids    <- read_bvid(paste0(g, "_binary_matrix_HQRoundStatisticalFiltered.tsv"))
qc_ids        <- read_bvid(paste0(g, "_binary_matrix_HQRoundStatisticalFilteredPlusQCFiltered.tsv"))
depth_ids     <- read_bvid(paste0(g, "_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered.tsv"))
phylo_ids     <- read_bvid(paste0(g, "_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv"))

vid <- master$VariantId
master$VAFSplit_Class <- ifelse(vid %in% singleton_ids, "singleton",
                          ifelse(vid %in% shared_ids, "shared", "none"))
master$VAFSplit_SharedRetained <- ifelse(vid %in% retained_ids, "Yes",
                                   ifelse(vid %in% shared_ids, "No", "."))
yn <- function(set) ifelse(vid %in% set, "Yes", "No")
master$VAFSplit_HQStat       <- yn(hqstat_ids)
master$VAFSplit_DepthPass    <- yn(depth_ids)
master$VAFSplit_QCPass       <- yn(qc_ids)
master$VAFSplit_ForPhylogeny <- yn(phylo_ids)

write.table(master, paste0("variant_master_filter_table_", group_id, ".tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Master table written:", nrow(master), "variants x", ncol(master), "columns\n")
cat("VAF-split membership added (HQStat/QC/Depth/ForPhylogeny:",
    length(hqstat_ids), length(qc_ids), length(depth_ids), length(phylo_ids), ")\n")

# Sections 7-9 (filter tracking table, markdown report, funnel plot) live in
# variant_filter_funnel.R (module VARIANT_FILTER_FUNNEL).
