# Shared QC scoring for the basej-* QC composition plots.
#
# WHY THIS EXISTS
# ---------------
# Every QC composition script used to score cells by COUNTING how many fixed thresholds a
# cell passed (`Mat` + `rowSums`). That has two defects:
#
#   * failures are interchangeable - missing PreSeq by 5% costs the same one point as CNV
#     MAPD being 3x over threshold, so the resulting 1-5 number is not an ordering;
#   * the count is displayed as an ordered tier, so a cell that fails a gate outright still
#     renders in the second band. On real data 53 of 84 cells scoring 4 of 5 were failing
#     PreSeq, including 7 cells with PreSeq of exactly zero.
#
# Here a tier requires EVERY one of its conditions to hold (conjunctive AND), so a cell can
# never rank above its weakest metric. Tier specs are nested by construction - tier 5's
# thresholds are tighter than tier 4's on the same metrics - so the ordering is meaningful.
#
# Tier 4 is deliberately set to each assay's established gate values, so `tier >= 4` stays
# identical to the old "all gates pass" (CompositeScore == 5) and published counts reconcile.
#
# DISPLAY BANDS
#   5, 4 -> Pass        4 = meets the established gates, 5 = tighter still
#   3    -> Borderline  relaxed thresholds; usable if capacity allows
#   2, 1 -> Fail        2 = too little data to judge, 1 = fails outright
#
# Used by dna_qc_plot.R, methyl_qc_plot.R, rna_qc_plot.R, wgs_qc_plot.R, wes_qc_plot.R
# and basej-hifiqc/bin/hifiqc_composition_plot.R.

QC_TIER_LABEL <- c(
  "5" = "5 Excellent quality",
  "4" = "4 Good quality",
  "3" = "3 Borderline",
  "2" = "2 Inconclusive",
  "1" = "1 Not recommended"
)

QC_TIER_BAND <- c(
  "5 Excellent quality" = "Pass",
  "4 Good quality"      = "Pass",
  "3 Borderline"        = "Borderline",
  "2 Inconclusive"      = "Fail",
  "1 Not recommended"   = "Fail"
)

QC_BAND_ORDER <- c("Pass", "Borderline", "Fail")
QC_TIER_ORDER <- unname(QC_TIER_LABEL)

# BioSkryb brand palette; the ploidy heatmap keeps its own diverging ramp.
QC_BAND_COLORS <- c(Pass = "#12284C", Borderline = "#A0CC2C", Fail = "#F45D34")
QC_TIER_COLORS <- c(
  "5 Excellent quality" = "#12284C",
  "4 Good quality"      = "#1082A2",
  "3 Borderline"        = "#A0CC2C",
  "2 Inconclusive"      = "#777776",
  "1 Not recommended"   = "#F45D34"
)


#' Depth floor below which a cell's metrics cannot be trusted.
#'
#' Coverage-derived metrics degrade with sequencing depth, not just cell quality: on the
#' validation cohort the PreSeq pass rate ran 0% / 41% / 87% / 98% / 97% across the
#' <0.3M / 0.3-0.6M / 0.6-1.0M / 1.0-1.5M / >1.5M read bands. A cell judged below that is
#' being scored on its depth. Such cells become tier 2 (Inconclusive) - unknown, not bad.
#'
#' Two modes:
#'   absolute   - `min_value`, when the assay has an established number (Illumina low-pass DNA)
#'   relative   - `median_frac` x the run median, when it does not (HiFi, RNA)
#'
#' The relative mode is deliberately median-relative rather than a fixed bottom percentile: a
#' percentile condemns a constant fraction of cells even when the whole run is uniformly good.
#'
#' @return logical vector, TRUE where the cell is below the floor
qc_depth_floor <- function(values, min_value = NULL, median_frac = NULL) {
  v <- qc_numeric(values)
  if (!is.null(min_value)) {
    return(is.na(v) | v < min_value)
  }
  if (!is.null(median_frac)) {
    med <- stats::median(v, na.rm = TRUE)
    if (!is.finite(med)) return(rep(TRUE, length(v)))
    return(is.na(v) | v < median_frac * med)
  }
  rep(FALSE, length(v))
}


#' Coerce a metrics column to numeric, tolerating the formats the pipelines actually write.
#'
#' Several metrics tables are built for MultiQC display as much as for machine reading, so read
#' counts arrive with thousands separators ("1,234,567" - see basej-methylqc's TotalReads).
#' Plain as.numeric() turns those into NA, which under `required_cols` would quietly mark every
#' cell's depth as unmeasured and push a whole plate to Inconclusive. Percent signs and stray
#' whitespace get the same treatment. Every numeric read in this file goes through here so the
#' trap cannot be reintroduced one call site at a time.
qc_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  v <- as.character(x)
  v <- gsub("[,%[:space:]]", "", v)
  v[v %in% c("", "NA", "NaN", "NULL", ".", "-")] <- NA
  suppressWarnings(as.numeric(v))
}


#' Apply a tier spec's condition for one column.
#'
#' gt/lt are strict, ge/le inclusive. Both forms exist because the original scripts were not
#' consistent - the DNA gates were written `>` while the WGS/WES read and coverage gates were
#' written `>=`. For an integer metric such as total_reads the boundary is reachable, so
#' collapsing the two would silently move cells sitting exactly on a gate.
#'
#' A condition may carry more than one comparison, all of which must hold, which is how a
#' two-sided band is written: `c(gt = 40, lt = 90)` for "between 40 and 90". Methylation needs
#' this - too little AND too much conversion both indicate a failed library - and expressing it
#' as a range beats the alternative of duplicating the column under a second name.
qc_cmp <- function(v, cond) {
  ok <- rep(TRUE, length(v))
  for (i in seq_along(cond)) {
    op <- names(cond)[i]
    thr <- cond[[i]]
    ok <- ok & switch(op,
      gt = v >  thr,
      ge = v >= thr,
      lt = v <  thr,
      le = v <= thr,
      stop(sprintf("qc_cmp: unknown comparison '%s' (expected gt/ge/lt/le)", op))
    )
  }
  ok
}


#' Do all of a tier's conditions hold?
#'
#' @param spec named list; each element is c(gt=x), c(ge=x), c(lt=x) or c(le=x) for that column
#' @return logical vector; NA in any required column yields FALSE (unmeasured is not acceptable)
qc_tier_meets <- function(df, spec) {
  ok <- rep(TRUE, nrow(df))
  for (col in names(spec)) {
    if (!col %in% colnames(df)) {
      stop(sprintf("qc_tier_meets: column '%s' required by the tier spec is absent", col))
    }
    v <- qc_numeric(df[[col]])
    passed <- qc_cmp(v, spec[[col]])
    passed[is.na(passed)] <- FALSE
    ok <- ok & passed
  }
  ok
}


#' Assign the 1-5 quality tier.
#'
#' Tiers are tested 5 -> 4 -> 3, first match wins. Quality is tested BEFORE the depth floor on
#' purpose: low depth depresses coverage metrics, so a cell that clears a tier on few reads has
#' cleared it against the odds and keeps its tier. Only cells reaching no tier AND sitting below
#' the floor become Inconclusive, where the honest statement is "cannot tell".
#'
#' @param tiers named list "5"/"4"/"3", each a spec for qc_tier_meets()
#' @param below_floor logical vector from qc_depth_floor(), or NULL to disable tier 2
#' @param required_cols columns that must be non-NA for any tier above 1 (e.g. the CNV metrics)
#' @return data.frame(SampleId, Tier, TierLabel, Band, BlockingMetric)
qc_assign_tiers <- function(df, tiers, below_floor = NULL, required_cols = character(),
                            sample_col = "SampleId") {
  n <- nrow(df)
  if (is.null(below_floor)) below_floor <- rep(FALSE, n)

  # An unmeasured required metric can only ever be tier 1. CNV is not optional for the QC
  # pipelines whose verdict depends on it, so a missing value must not silently pass.
  measured <- rep(TRUE, n)
  for (col in required_cols) {
    if (!col %in% colnames(df)) {
      stop(sprintf("qc_assign_tiers: required column '%s' is absent", col))
    }
    measured <- measured & !is.na(qc_numeric(df[[col]]))
  }

  tier <- rep(1L, n)
  for (k in c("3", "4", "5")) {
    if (!is.null(tiers[[k]])) {
      hit <- measured & qc_tier_meets(df, tiers[[k]])
      tier[hit] <- as.integer(k)
    }
  }
  # Tier 2 means "cannot judge", and there are two ways to land there: too little data to
  # trust the metrics, or a required metric that was never produced for this cell (e.g. Ginkgo
  # did not emit CNV for it). Neither is evidence the cell is bad, so neither may be reported
  # as tier 1.
  tier[tier == 1L & (below_floor | !measured)] <- 2L

  # Which condition blocked promotion to the next tier up - so no score is unexplained.
  blocking <- rep("", n)
  next_up <- c("1" = "3", "2" = NA, "3" = "4", "4" = "5", "5" = NA)
  for (i in seq_len(n)) {
    if (tier[i] == 2L) {
      blocking[i] <- if (!measured[i]) "metric_unmeasured" else "insufficient_data"
      next
    }
    if (!measured[i])  { blocking[i] <- "metric_unmeasured"; next }
    tgt <- next_up[[as.character(tier[i])]]
    if (is.na(tgt) || is.null(tiers[[tgt]])) next
    spec <- tiers[[tgt]]
    fails <- character()
    for (col in names(spec)) {
      v <- qc_numeric(df[[col]][i])
      passed <- qc_cmp(v, spec[[col]])
      if (is.na(passed) || !passed) fails <- c(fails, col)
    }
    blocking[i] <- paste(fails, collapse = "+")
  }

  label <- unname(QC_TIER_LABEL[as.character(tier)])
  data.frame(
    SampleId = df[[sample_col]],
    Tier = tier,
    TierLabel = factor(label, levels = QC_TIER_ORDER),
    Band = factor(unname(QC_TIER_BAND[label]), levels = QC_BAND_ORDER),
    BlockingMetric = blocking,
    stringsAsFactors = FALSE
  )
}


#' Fail the run when a metric the verdict depends on is missing for every sample.
#'
#' Per-cell absence is handled by `required_cols` (that cell becomes Inconclusive). Absence
#' across the whole run means an upstream step did not produce the metric, and silently
#' scoring everything as Inconclusive would hide a broken run.
qc_require_metric <- function(df, cols, what = "CNV") {
  for (col in cols) {
    if (!col %in% colnames(df)) {
      stop(sprintf("%s metric '%s' is absent - cannot score. Did the upstream step run?", what, col))
    }
    if (all(is.na(qc_numeric(df[[col]])))) {
      stop(sprintf("%s metric '%s' is NA for every sample - cannot score. Did the upstream step run?",
                   what, col))
    }
  }
  invisible(TRUE)
}


#' Per-sample tier table for the pipeline's *_ConsensusScores.txt output.
qc_scores_table <- function(tiers_df, extra = NULL) {
  out <- data.frame(
    SampleId = tiers_df$SampleId,
    QC_Tier = tiers_df$Tier,
    QC_Label = as.character(tiers_df$TierLabel),
    QC_Band = as.character(tiers_df$Band),
    BlockingMetric = tiers_df$BlockingMetric,
    stringsAsFactors = FALSE
  )
  if (!is.null(extra)) out <- cbind(out, extra)
  out[order(-out$QC_Tier, out$SampleId), ]
}
