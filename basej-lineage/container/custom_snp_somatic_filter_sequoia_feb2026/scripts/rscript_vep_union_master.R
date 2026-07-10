#!/usr/bin/env Rscript
## Build master VEP table from one or more per-sample CSQ raw TSVs (one best CSQ row per VARIANT_ID).
## Cohort mode may pass a single representative TSV when all VCFs share the same variant set.
suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("optparse required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")
  library(data.table)
})

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--raw_tsvs", type = "character",
                        help = "Space-separated *_vep_csq_raw.tsv paths (one representative file is enough when cohort VCFs share the same sites)"),
  optparse::make_option("--info_header", type = "character",
                        help = "Path to file with ##INFO CSQ (and optional VEP_Existing_variation) lines"),
  optparse::make_option("--out", type = "character", help = "Output master TSV path")
)))

paths <- strsplit(opt$raw_tsvs, " ", fixed = TRUE)[[1]]
paths <- paths[nzchar(paths)]
if (!length(paths)) stop("No --raw_tsvs")
if (!file.exists(opt$info_header)) stop("Missing --info_header: ", opt$info_header)

parse_csq_format_from_header <- function(header_path) {
  lines <- readLines(header_path, warn = FALSE)
  csq_line <- lines[grepl("^##INFO=<ID=CSQ", lines)][1]
  if (is.na(csq_line)) stop("No ##INFO=<ID=CSQ line in ", header_path)
  m <- regmatches(csq_line, regexpr("Format:[^>]+", csq_line))
  if (!length(m)) stop("Could not parse CSQ Format from header")
  fmt <- sub("^Format:[[:space:]]*", "", m)
  fmt <- gsub("\"$", "", fmt)
  strsplit(fmt, "|", fixed = TRUE)[[1]]
}

csq_cols <- parse_csq_format_from_header(opt$info_header)

explode_csq_row <- function(chrom, pos, ref, alt, sample, csq_str, vep_ex) {
  vid <- paste(chrom, pos, ref, alt, sep = "_")
  if (is.na(csq_str) || !nzchar(csq_str)) {
    empty <- as.data.table(as.list(setNames(rep(NA_character_, length(csq_cols)), csq_cols)))
    empty[, `:=`(
      sample_name = sample, CHROM = chrom, POS = as.integer(pos), REF = ref, ALT = alt,
      VARIANT_ID = vid, CSQ_raw = NA_character_, VEP_Existing_variation = vep_ex
    )]
    return(empty)
  }
  parts <- strsplit(csq_str, ",", fixed = TRUE)[[1]]
  n <- length(csq_cols)
  out <- lapply(parts, function(p) {
    v <- strsplit(p, "|", fixed = TRUE)[[1]]
    if (length(v) < n) v <- c(v, rep(NA_character_, n - length(v)))
    if (length(v) > n) v <- v[seq_len(n)]
    v
  })
  mat <- as.data.table(do.call(rbind, out))
  setnames(mat, csq_cols)
  mat[, `:=`(
    sample_name = sample,
    CHROM = chrom,
    POS = as.integer(pos),
    REF = ref,
    ALT = alt,
    VARIANT_ID = vid,
    CSQ_raw = parts,
    VEP_Existing_variation = vep_ex
  )]
  mat[]
}

impact_rank <- function(x) {
  x <- toupper(as.character(x))
  ifelse(x == "HIGH", 4L,
         ifelse(x == "MODERATE", 3L,
                ifelse(x == "LOW", 2L,
                       ifelse(x == "MODIFIER", 1L, 0L))))
}

clin_rank <- function(clin) {
  if (is.na(clin) || !nzchar(clin)) return(0L)
  s <- tolower(clin)
  sc <- 0L
  if (grepl("pathogenic", s) && !grepl("likely_pathogenic|benign", s)) sc <- sc + 8L
  if (grepl("likely_pathogenic", s)) sc <- sc + 7L
  if (grepl("uncertain|vus|conflicting", s)) sc <- sc + 2L
  if (grepl("likely_benign", s)) sc <- sc - 2L
  if (grepl("benign", s) && !grepl("likely_benign|pathogenic", s)) sc <- sc - 4L
  sc
}

read_one_raw <- function(fp) {
  dt <- fread(fp, sep = "\t", header = TRUE, quote = "", na.strings = c("", "NA"))
  req <- c("sample_name", "CHROM", "POS", "REF", "ALT", "CSQ", "VEP_Existing_variation")
  if (!all(req %in% names(dt))) {
    stop("Bad columns in ", fp, ": need ", paste(req, collapse = ", "))
  }
  rows <- vector("list", nrow(dt))
  for (i in seq_len(nrow(dt))) {
    rows[[i]] <- explode_csq_row(
      dt$CHROM[i], dt$POS[i], dt$REF[i], dt$ALT[i],
      dt$sample_name[i], dt$CSQ[i], dt$VEP_Existing_variation[i]
    )
  }
  rbindlist(rows, fill = TRUE)
}

all_dt <- rbindlist(lapply(paths, read_one_raw), fill = TRUE)

## Per VARIANT_ID (+ exploded CSQ row): rank and keep single best annotation
if (!"IMPACT" %in% names(all_dt)) all_dt[, IMPACT := NA_character_]
if (!"Consequence" %in% names(all_dt)) all_dt[, Consequence := NA_character_]
if (!"CLIN_SIG" %in% names(all_dt)) all_dt[, CLIN_SIG := NA_character_]

all_dt[, `:=`(
  rank_impact = impact_rank(IMPACT),
  rank_clin = sapply(CLIN_SIG, clin_rank)
)]
all_dt[, rank_score := rank_impact * 10L + rank_clin]
setorder(all_dt, VARIANT_ID, -rank_score, -rank_impact)
best <- all_dt[, .SD[1], by = VARIANT_ID]
best[, c("rank_impact", "rank_clin", "rank_score") := NULL]

## Stable column order: ids first, then CSQ fields, then provenance
lead <- c("VARIANT_ID", "CHROM", "POS", "REF", "ALT", "sample_name",
          "VEP_Existing_variation", "CSQ_raw")
rest <- setdiff(names(best), lead)
setcolorder(best, c(intersect(lead, names(best)), sort(rest)))

fwrite(best, opt$out, sep = "\t", quote = FALSE)
message("[vep_master_csq] Wrote ", nrow(best), " rows (one per VARIANT_ID) -> ", opt$out)
