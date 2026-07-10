#!/usr/bin/env Rscript
## Add clinical relevance columns + mandatory Verdict (5 tiers) to master VEP union TSV.
suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("optparse required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")
  library(data.table)
})

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option(c("-i", "--input"), type = "character", dest = "input",
                          help = "Master union VEP TSV"),
  optparse::make_option("--out", type = "character", help = "Output TSV with extra columns")
)))

if (is.null(opt$input) || !nzchar(opt$input)) stop("--input is required")
dt <- fread(opt$input, sep = "\t", header = TRUE, quote = "", na.strings = c("", "NA"))

safe_lower <- function(x) tolower(as.character(x))
ch <- function(v) as.character(v)

## SYMBOL (HGNC-style) + Gene (Ensembl ID) → unique gene tokens; intergenic "A-B" → two genes
parse_gene_tokens <- function(sym, gene, cons) {
  s <- trimws(ch(sym))
  g <- trimws(ch(gene))
  c_low <- tolower(ch(cons))
  tok <- character(0)
  if (!is.na(s) && nzchar(s) && !s %in% c("-", ".", "NA")) {
    if (grepl("intergenic", c_low) && grepl("-", s, fixed = TRUE) &&
        !grepl(",", s, fixed = TRUE)) {
      parts <- trimws(strsplit(s, "-", fixed = TRUE)[[1]])
      tok <- c(tok, parts[nzchar(parts) & !parts %in% c("-", ".")])
    } else {
      parts <- trimws(strsplit(s, ",|;|&|\\|", perl = TRUE)[[1]])
      tok <- c(tok, parts[nzchar(parts) & !parts %in% c("-", ".")])
    }
  }
  if (!is.na(g) && nzchar(g) && !g %in% c("-", ".", "NA")) {
    parts <- trimws(strsplit(g, ",|;|&", perl = TRUE)[[1]])
    tok <- c(tok, parts[nzchar(parts) & !parts %in% c("-", ".")])
  }
  unique(tok[nzchar(tok)])
}

sym_col <- if ("SYMBOL" %in% names(dt)) dt[["SYMBOL"]] else rep(NA_character_, nrow(dt))
gene_col <- if ("Gene" %in% names(dt)) dt[["Gene"]] else rep(NA_character_, nrow(dt))
cons_col <- if ("Consequence" %in% names(dt)) dt[["Consequence"]] else rep(NA_character_, nrow(dt))

genes_per_row <- mapply(parse_gene_tokens, sym_col, gene_col, cons_col, SIMPLIFY = FALSE)
dt[, Genes_affected := vapply(genes_per_row, function(z) {
  if (!length(z)) "NoGene" else paste(z, collapse = ",")
}, character(1))]
dt[, NumberGenesAffected := vapply(genes_per_row, length, integer(1))]

ex_csq <- if ("Existing_variation" %in% names(dt)) ch(dt[["Existing_variation"]]) else rep(NA_character_, nrow(dt))
ex_vep <- if ("VEP_Existing_variation" %in% names(dt)) ch(dt[["VEP_Existing_variation"]]) else rep(NA_character_, nrow(dt))
dt[, Existing_variation_combined := vapply(seq_len(nrow(dt)), function(i) {
  a <- trimws(ex_csq[i])
  b <- trimws(ex_vep[i])
  ok_a <- !is.na(a) && nzchar(a)
  ok_b <- !is.na(b) && nzchar(b)
  if (!ok_a && !ok_b) return("NoRecord")
  if (ok_a && ok_b) paste(a, b, sep = "|")
  else if (ok_a) a
  else b
}, character(1))]

has_impact <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else NA_character_
has_cons <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else NA_character_
has_clin <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else NA_character_
has_maxaf <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else NA_character_

impact <- safe_lower(has_impact(dt, "IMPACT"))
cons <- safe_lower(has_cons(dt, "Consequence"))
clin <- safe_lower(has_clin(dt, "CLIN_SIG"))
max_af <- suppressWarnings(as.numeric(has_maxaf(dt, "MAX_AF")))

## Helper flags (order: likely_pathogenic before pathogenic substring match)
clin_lp <- grepl("likely_pathogenic", clin)
clin_path <- grepl("pathogenic", clin) & !clin_lp & !grepl("benign", clin)
clin_b <- grepl("(^|,)\\s*benign|likely_benign", clin) | grepl("benign/likely_benign", clin)
clin_conf <- grepl("conflicting|uncertain|vus", clin)

high_trunc <- (impact == "high") & grepl(
  "frameshift|stop_gained|stop_lost|start_lost|splice_donor|splice_acceptor|nonsense",
  cons
)
missense_mod <- grepl("missense|inframe", cons)
coding <- grepl(
  "missense|frameshift|stop_gained|splice|intron|synonymous|utr|upstream|downstream|intergenic|non_coding|inframe",
  cons
)

syn_only <- grepl("^synonymous_variant$|^synonymous_variant,|synonymous_variant", cons) &
  !grepl("missense|frameshift", cons)
intergenic <- grepl("intergenic_variant", cons) & !grepl("missense|frameshift|splice", cons)

af_low <- is.na(max_af) | max_af <= 0.001

## Summaries
dt[, Consequence_class := fifelse(
  high_trunc, "truncating_splice",
  fifelse(missense_mod, "missense_inframe",
          fifelse(syn_only, "synonymous",
                  fifelse(intergenic, "intergenic_distal", "other"))))
]

dt[, CLIN_SIG_summary := ifelse(
  clin_path | clin_lp,
  ifelse(clin_lp, "likely_pathogenic", "pathogenic"),
  ifelse(clin_b, "benign_family",
         ifelse(clin_conf, "conflicting_or_vus", "unreported_or_other")))
]

dt[, MAX_AF_numeric := max_af]

dt[, Coding_region_flag := ifelse(coding, "coding_or_near_gene", "other")]

ev <- rep(NA_character_, nrow(dt))
ev[clin_path | clin_lp] <- "clinical_pathogenic"
ev[is.na(ev) & high_trunc] <- "high_impact_truncating"
ev[is.na(ev) & (impact %in% c("high", "moderate")) & missense_mod & !clin_b] <- "damaging_coding"
ev[is.na(ev) & clin_conf] <- "conflicting_classification"
ev[is.na(ev) & syn_only & af_low] <- "synonymous_rare"
ev[is.na(ev) & intergenic] <- "non_coding_intergenic"
ev[is.na(ev)] <- "see_consequence"
dt[, Pathogenicity_evidence := ev]

## Genes_affected / NumberGenesAffected / Existing_variation_combined already added above

## Verdict: 5 succinct categories — assign lowest priority first, then override upward
n <- nrow(dt)
ver <- rep("Other", n)
is_unc <- is.na(cons) & is.na(clin) & is.na(impact)
ver[is_unc] <- "Uncertain"

## Non-coding / intronic / UTR without splice — usually lowest clinical interest
ver[grepl("downstream|upstream|utr|intron|intergenic|non_coding", cons) &
      !grepl("splice_donor|splice_acceptor|frameshift|stop_gained", cons) &
      !clin_path & !clin_lp] <- "Other"

ver[clin_b | (syn_only & (is.na(max_af) | max_af <= 0.05))] <- "Benign"
ver[clin_conf & !clin_path & !clin_lp] <- "Uncertain"

del_mask <- !clin_b & (
  clin_lp | high_trunc |
    (impact == "high" & (missense_mod | grepl("splice", cons))) |
    (impact == "moderate" & missense_mod & af_low)
)
ver[del_mask] <- "Deleterious"

path_mask <- clin_path | (clin_lp & (high_trunc | impact == "high"))
ver[path_mask] <- "Pathogenic"

dt[, Verdict := ver]

fwrite(dt, opt$out, sep = "\t", quote = FALSE)
message("[vep_relevance] Wrote ", nrow(dt), " rows -> ", opt$out)
