#!/usr/bin/env python3
"""Shared QC tier scoring for the Rust-rewrite (-rs) pipelines.

A Python port of containers/qc_plots/scripts/qc_scoring.R, so the -rs pipelines produce
the same qc_score / qc_status as the R composition scripts do for their reference
pipelines. The -rs pipelines write Parquet directly from the Rust single pass and never
invoke R, so the logic has to exist here too.

Kept a deliberate line-by-line mirror of the R rather than a re-imagining: the two are
cross-checked against each other in tests/test_qc_scoring.py, and that check is only
meaningful if they are independently readable statements of the same rules. When you
change one, change the other and re-run both tests.

WHY TIERS AND NOT A COUNT
-------------------------
The reference pipelines used to score a cell by counting how many fixed thresholds it
passed. That makes failures interchangeable - missing the read target by 5% costs the same
single point as CNV MAPD being 3x over - and the count was then displayed as an ordered
band, so a cell failing a gate outright still rendered above cells that passed everything.
On 1705 real cells, 53 of the 84 scoring 4 of 5 were failing PreSeq, 7 of them with PreSeq
of exactly zero.

Here a tier requires EVERY one of its conditions (conjunctive AND), so a cell can never
rank above its weakest metric. Tier specs are nested by construction, so the ordering
means something.

  5, 4 -> PASS        4 = meets the established gates, 5 = tighter still
  3    -> Borderline  relaxed thresholds; usable if capacity allows
  2, 1 -> FAIL        2 = too little data to judge, 1 = fails outright

Tier 4 uses each assay's established gate values, so `tier >= 4` is identical to the old
"all gates pass" and published counts reconcile.

AI-assisted implementation.
"""
from __future__ import annotations

import math
import re
from typing import Iterable, Mapping, Sequence

TIER_LABEL = {
    5: "5 Excellent quality",
    4: "4 Good quality",
    3: "3 Borderline",
    2: "2 Inconclusive",
    1: "1 Not recommended",
}

TIER_BAND = {
    "5 Excellent quality": "PASS",
    "4 Good quality": "PASS",
    "3 Borderline": "Borderline",
    "2 Inconclusive": "FAIL",
    "1 Not recommended": "FAIL",
}

TIER_ORDER = [TIER_LABEL[t] for t in (5, 4, 3, 2, 1)]
BAND_ORDER = ["PASS", "Borderline", "FAIL"]

_NUM_STRIP = re.compile(r"[,%\s]")
_NULLISH = {"", "NA", "NAN", "NONE", "NULL", ".", "-"}


def qc_numeric(value) -> float | None:
    """Coerce a metrics value to float, tolerating the formats the pipelines write.

    Mirrors qc_numeric() in qc_scoring.R. Read counts reach the R scripts formatted for
    MultiQC display ("1,234,567"), and plain float() raises on those. Returning None for a
    genuinely absent value matters: "not measured" and "measured and bad" lead to different
    tiers, so they must not collapse into each other.
    """
    if value is None:
        return None
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return None if (isinstance(value, float) and math.isnan(value)) else float(value)
    s = _NUM_STRIP.sub("", str(value))
    if s.upper() in _NULLISH:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def qc_cmp(value, cond: Mapping[str, float]) -> bool:
    """Apply a tier spec's condition for one metric.

    gt/lt are strict, ge/le inclusive. Both exist because the reference scripts were not
    consistent - the DNA gates were written `>`, the WGS/WES read and coverage gates `>=` -
    and for an integer metric such as total_reads the boundary is reachable, so collapsing
    them would silently move samples sitting exactly on a gate.

    A condition may carry several comparisons, all of which must hold, which is how a
    two-sided band is written: {"gt": 40, "lt": 90}.

    An unmeasured value never satisfies a condition.
    """
    v = qc_numeric(value)
    if v is None:
        return False
    for op, thr in cond.items():
        if op == "gt":
            ok = v > thr
        elif op == "ge":
            ok = v >= thr
        elif op == "lt":
            ok = v < thr
        elif op == "le":
            ok = v <= thr
        else:
            raise ValueError(f"qc_cmp: unknown comparison {op!r} (expected gt/ge/lt/le)")
        if not ok:
            return False
    return True


def tier_meets(row: Mapping, spec: Mapping[str, Mapping[str, float]]) -> bool:
    """Do all of a tier's conditions hold for this row?"""
    return all(qc_cmp(row.get(col), cond) for col, cond in spec.items())


def _blocking_metrics(row: Mapping, spec: Mapping[str, Mapping[str, float]]) -> list[str]:
    return [col for col, cond in spec.items() if not qc_cmp(row.get(col), cond)]


def depth_floor(value, *, min_value: float | None = None,
                median: float | None = None, median_frac: float | None = None) -> bool:
    """Is this sample below the depth at which its coverage metrics can be trusted?

    Coverage-derived metrics degrade with sequencing depth, not just with cell quality: on
    the validation cohort the PreSeq pass rate ran 0% / 41% / 87% / 98% / 97% across the
    <0.3M / 0.3-0.6M / 0.6-1.0M / 1.0-1.5M / >1.5M read bands. A sample judged below that is
    being scored on its depth, so it becomes tier 2 (Inconclusive) - unknown, not bad.

    Absolute mode (min_value) where the assay has an established target; median-relative
    mode (median + median_frac) where it does not. Median-relative rather than a bottom
    percentile: a percentile condemns a constant fraction of samples even when the whole run
    is uniformly good.
    """
    v = qc_numeric(value)
    if min_value is not None:
        return v is None or v < min_value
    if median_frac is not None:
        if median is None or not math.isfinite(median):
            return True
        return v is None or v < median_frac * median
    return False


def assign_tier(row: Mapping, tiers: Mapping[str, Mapping], *,
                below_floor: bool = False,
                required: Sequence[str] = ()) -> dict:
    """Assign the 1-5 tier for one sample.

    Tiers are tested 5 -> 4 -> 3 and the best match wins. Quality is tested BEFORE the depth
    floor on purpose: low depth depresses coverage metrics, so a sample that clears a tier on
    few reads has cleared it against the odds and keeps its tier. Only samples reaching no
    tier AND sitting below the floor become Inconclusive, where the honest statement is
    "cannot tell".

    `required` names metrics that must be present for any tier above 1 (the CNV metrics, for
    the pipelines whose verdict depends on them). A sample missing one is Inconclusive, not
    failed - absence is not evidence of poor quality.
    """
    measured = all(qc_numeric(row.get(c)) is not None for c in required)

    tier = 1
    if measured:
        for k in (3, 4, 5):
            spec = tiers.get(str(k))
            if spec and tier_meets(row, spec):
                tier = k

    # Tier 2 means "cannot judge", reached two ways: too little data to trust the metrics, or
    # a required metric that was never produced. Neither is evidence the sample is bad, so
    # neither may be reported as tier 1.
    if tier == 1 and (below_floor or not measured):
        tier = 2

    if tier == 2:
        blocking = "metric_unmeasured" if not measured else "insufficient_data"
    elif not measured:
        blocking = "metric_unmeasured"
    else:
        nxt = {1: "3", 3: "4", 4: "5"}.get(tier)
        blocking = "+".join(_blocking_metrics(row, tiers[nxt])) if nxt and tiers.get(nxt) else ""

    label = TIER_LABEL[tier]
    return {"qc_score": tier, "qc_label": label, "qc_status": TIER_BAND[label],
            "blocking_metric": blocking}


def require_metric(rows: Iterable[Mapping], cols: Sequence[str], what: str = "CNV") -> None:
    """Raise when a metric the verdict depends on is missing for EVERY sample.

    Per-sample absence is handled by `required` (that sample becomes Inconclusive). Absence
    across a whole run means an upstream step did not produce the metric, and silently
    scoring everything Inconclusive would hide a broken run.
    """
    rows = list(rows)
    if not rows:
        return
    for col in cols:
        if all(qc_numeric(r.get(col)) is None for r in rows):
            raise SystemExit(
                f"ERROR: {what} metric '{col}' is NA for every sample - cannot score. "
                "Did the upstream step run?")


# --------------------------------------------------------------------------- tier specs
# Mirrors qc_scoring_specs.R. These are BioSkryb-internal values and depend on cell type and
# biological context: the right cutoff for a tumour biopsy is not the right cutoff for a cell
# line. Tier 4 mirrors each assay's established gates so `tier >= 4` matches the historical
# "all gates pass".

def spec_dna(cutoff_preseq=3.5e9, cutoff_chimeras=0.20, cutoff_chrm=0.20,
             cutoff_cnv_mapd=0.25, cutoff_cnv_sk=0.25,
             tight_cnv_mapd=0.20, tight_cnv_sk=0.20,
             borderline_preseq=2.5e9, borderline_cnv=0.35):
    """Single-cell low-pass DNA (basej-dnaqc / basej-dnaqc)."""
    return {
        "5": {"preseq_count": {"gt": cutoff_preseq},
              "pct_chimeras": {"lt": cutoff_chimeras},
              "pct_chrm": {"lt": cutoff_chrm},
              "ginkgo_cnv_mapd": {"lt": tight_cnv_mapd},
              "ginkgo_cnv_skew": {"lt": tight_cnv_sk}},
        "4": {"preseq_count": {"gt": cutoff_preseq},
              "pct_chimeras": {"lt": cutoff_chimeras},
              "pct_chrm": {"lt": cutoff_chrm},
              "ginkgo_cnv_mapd": {"lt": cutoff_cnv_mapd},
              "ginkgo_cnv_skew": {"lt": cutoff_cnv_sk}},
        "3": {"preseq_count": {"gt": borderline_preseq},
              "ginkgo_cnv_mapd": {"lt": borderline_cnv},
              "ginkgo_cnv_skew": {"lt": borderline_cnv}},
    }


def spec_wgs(cutoff_num_reads=50e6, cutoff_pct_dup=0.25, cutoff_pct_chim=0.15,
             cutoff_1x=0.9, cutoff_5x=0.7,
             tight_pct_dup=0.15, tight_1x=0.95, tight_5x=0.85,
             borderline_reads=25e6, borderline_1x=0.8, borderline_5x=0.5):
    """Full-depth WGS (basej-wgs / basej-wgsqc, wgs mode).

    Unlike low-pass DNA, the read count IS a gate here: the assay has an established target
    and coverage breadth is read against it. Tier 5 tightens coverage breadth and
    duplication; chimeras stay at the gate value because on real data the chimera gate is
    almost never the binding constraint.
    """
    return {
        "5": {"total_reads": {"ge": cutoff_num_reads},
              "pct_duplication": {"lt": tight_pct_dup},
              "pct_chimeras": {"lt": cutoff_pct_chim},
              "pct_1x": {"ge": tight_1x},
              "pct_5x": {"ge": tight_5x}},
        "4": {"total_reads": {"ge": cutoff_num_reads},
              "pct_duplication": {"lt": cutoff_pct_dup},
              "pct_chimeras": {"lt": cutoff_pct_chim},
              "pct_1x": {"ge": cutoff_1x},
              "pct_5x": {"ge": cutoff_5x}},
        "3": {"total_reads": {"ge": borderline_reads},
              "pct_1x": {"ge": borderline_1x},
              "pct_5x": {"ge": borderline_5x}},
    }


def spec_wes(cutoff_num_reads=5e5, cutoff_10x=0.75, cutoff_zero_cov=0.05, cutoff_fold_80=5,
             tight_10x=0.90, tight_zero_cov=0.02, tight_fold_80=3,
             borderline_reads=2.5e5, borderline_10x=0.5, borderline_fold_80=8):
    """Exome (basej-wgs / basej-wgsqc, exome mode).

    Four gates rather than five, so the old maximum score was 4 - `tier >= 4` therefore still
    means "all gates pass". fold_80_base_penalty is capture uniformity, the metric that most
    often separates a usable exome from a good one.
    """
    return {
        "5": {"total_reads": {"ge": cutoff_num_reads},
              "pct_target_bases_10x": {"ge": tight_10x},
              "zero_cvg_targets_pct": {"lt": tight_zero_cov},
              "fold_80_base_penalty": {"lt": tight_fold_80}},
        "4": {"total_reads": {"ge": cutoff_num_reads},
              "pct_target_bases_10x": {"ge": cutoff_10x},
              "zero_cvg_targets_pct": {"lt": cutoff_zero_cov},
              "fold_80_base_penalty": {"lt": cutoff_fold_80}},
        "3": {"total_reads": {"ge": borderline_reads},
              "pct_target_bases_10x": {"ge": borderline_10x},
              "fold_80_base_penalty": {"lt": borderline_fold_80}},
    }
