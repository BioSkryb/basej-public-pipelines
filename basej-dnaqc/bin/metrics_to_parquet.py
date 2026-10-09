#!/usr/bin/env python3
"""
Aggregate dnaqc_summary Parquet writer for basej-dnaqc (rewrites.bio).

Collect-based rewrite (mirrors basej-methylqc-rs's methylqc_metrics_to_parquet.py):
discovers samples from the *.rustqc.json staged in the work dir and merges, per
sample, the bskryb-qc rustqc JSON (align/insert/coverage/chrM/quality-yield/gcbias/
preseq/lorenz-gini), dupblaster stats (dedup_*), fastp JSON (fastp_*), read counts,
the actual preseq tool ({sample}_preseq.txt override), and Ginkgo metrics/SegCopy/
cnvSummarizer output (ginkgo_*) into the SAME dnaqc_summary schema (column names +
types) as basej-dnaqc's QC_PLOTS. Keeping the schema identical means the
Athena/Iceberg table is unchanged and the two pipelines can be compared directly.

Emits:
  - dnaqc_summary/workspace=*/workflow_id=*/biosample=*/output.parquet (per biosample)
  - {sample}_selected_metrics_mqc.txt   (MultiQC custom-content table)
  - dnaqc_all_metrics.tsv               (flat summary of all samples)

QC scoring is the shared conjunctive 1-5 tier scheme (qc_scoring.py, a validated port of
containers/qc_plots/scripts/qc_scoring.R), so qc_score/qc_status match what basej-dnaqc
produces via dna_qc_plot.R for the same metrics.

Ginkgo ploidy is loaded from the initial Ginkgo metrics/SegCopy outputs, while
MAPD_CNV and SKEW_CNV are loaded from the R cnvSummarizer output. These metrics
feed the shared DNA QC scoring thresholds when Ginkgo runs.

NOTE: the Lorenz/Gini columns (gini_coefficient_index, roc_lorenz_curve,
total_covered_positions_of_genome, total_sequenced_bases,
total_investigated_genomic_positions) and the ginkgo_* columns are additive vs the
original dnaqc_summary and require Iceberg/Athena schema evolution before they are
queryable.

Attribution: reimplements metric parsing from Sentieon (AlignmentStat /
InsertSizeMetricAlgo / GCBias / CoverageMetrics / QualityYield), Picard
MarkDuplicates (via dupblaster, Fulcrum Genomics, MIT), fastp, and preseq.
AI-assisted; validated versions and known gaps are documented in the
container/module READMEs. Content rephrased for compliance.
"""
import argparse
import glob
import json
import os
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import qc_scoring

# Below this many reads the coverage-derived metrics measure depth rather than the cell, so it
# is reported Inconclusive instead of failed. Matches dna_qc_plot.R's cutoff_reads_floor.
DEPTH_FLOOR_READS = 1e6

# Modular chimera classification: bskryb-qc emits the chimera_cube; named mechanism
# classes are defined in chimera_classes.yaml and applied via chimera_classify.py
# (co-located in bin/). Falls back gracefully if unavailable.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from chimera_classify import classify_cube, load_ruleset
    _CHIM_RULESET = load_ruleset()
except Exception:
    classify_cube = None
    _CHIM_RULESET = None

CHIM_CLASSES = ["proper", "interchromosomal", "local_inverted_hairpin", "inverted_nonlocal",
                "everted_tandem", "large_insert", "split_read_only", "other"]


def chimera_class_fields(m):
    """Per-class chimera fractions (chim_class_<name>_pct) from the rustqc chimera_cube."""
    if classify_cube is None:
        return {}
    cube = m.get("chimera_cube")
    if not cube:
        return {}
    try:
        res = classify_cube(cube, _CHIM_RULESET)
        return {f"chim_class_{k}_pct": v["pct"] for k, v in res.items()}
    except Exception:
        return {}

# Type lists copied from basej-dnaqc so Parquet types stay identical (subset),
# plus the additive ginkgo_* columns (mirrors basej-methylqc-rs).
#
# qc_label and blocking_metric are additive too: the tier name, and the condition that held
# the cell back so no score is unexplained. Additive columns are safe for Iceberg, but the
# dnaqc_summary table schema must still be evolved before they are queryable - see the
# "Metrics Parquet schema changes" section of .kiro/steering/pipeline-dev.md.
_str_cols = ['biosample', 'dataset_id', 'pipeline', 'pipeline_version',
             'molecule_type', 'workflow_id', 'workspace', 'user',
             'ginkgo_cnv_mapd', 'ginkgo_cnv_skew', 'ginkgo_ploidy', 'qc_status',
             'qc_label', 'blocking_metric']
_double_cols = ['align_pf_mismatch_rate', 'align_pf_hq_error_rate', 'align_pf_indel_rate',
                'align_mean_read_length', 'align_pct_reads_aligned_in_pairs',
                'align_strand_balance', 'align_pct_chimeras',
                'align_chim_diff_contig_pct', 'align_chim_bad_orient_pct',
                'align_chim_large_insert_pct', 'align_chim_split_sa_pct',
                'chim_class_proper_pct', 'chim_class_interchromosomal_pct',
                'chim_class_local_inverted_hairpin_pct', 'chim_class_inverted_nonlocal_pct',
                'chim_class_everted_tandem_pct', 'chim_class_large_insert_pct',
                'chim_class_split_read_only_pct', 'chim_class_other_pct',
                'pct_aligned', 'pct_pf', 'pct_error', 'pct_chimeras', 'pct_chrm',
                'sentieon_q20_rate', 'sentieon_q30_rate',
                'gc_at_dropout', 'gc_gc_dropout', 'gc_gc_nc_0_19', 'gc_gc_nc_20_39',
                'gc_gc_nc_40_59', 'gc_gc_nc_60_79', 'gc_gc_nc_80_100',
                'insert_median', 'insert_size', 'insert_mean', 'insert_std',
                'insert_median_absolute_deviation', 'insert_min', 'insert_max',  # double, as basej-dnaqc
                'gini_coefficient_index', 'roc_lorenz_curve', 'ginkgo_average_ploidy']
_double_cols += ['fastp_q20_rate', 'fastp_q30_rate', 'fastp_gc_content',
                 'fastp_q30_rate_after', 'fastp_duplication_rate',
                 'dedup_pct_duplication', 'pct_duplication',
                 'pct_optical_duplicates', 'pct_pcr_duplicates']
_bigint_cols = ['align_total_reads', 'align_pf_reads', 'align_pf_reads_aligned',
                'align_pf_hq_aligned_reads', 'align_pf_hq_aligned_bases',
                'align_pf_hq_aligned_q20_bases', 'align_reads_aligned_in_pairs',
                'align_chim_diff_contig', 'align_chim_bad_orient',
                'align_chim_large_insert', 'align_chim_split_sa',
                'insert_read_pairs', 'cov_total_bases', 'cov_chrm_bases',
                'sentieon_total_bases', 'gc_total_clusters', 'gc_aligned_reads',
                'preseq_count',
                'total_covered_positions_of_genome', 'total_sequenced_bases',
                'total_investigated_genomic_positions',
                'total_reads', 'final_reads',
                'fastp_total_reads', 'fastp_total_bases', 'fastp_read1_mean_length',
                'fastp_read2_mean_length', 'fastp_reads_after_filter',
                'fastp_adapter_trimmed_reads', 'fastp_adapter_trimmed_bases',
                'dedup_read_pairs_examined', 'dedup_read_pair_duplicates',
                'dedup_unpaired_reads_examined', 'dedup_unpaired_read_duplicates',
                'dedup_read_pair_optical_duplicates', 'dedup_estimated_library_size',
                'qc_score']
_bool_cols = ['subsampled']


def parse_fastp(path):
    """Extract fastp_* fields from a fastp JSON (before/after filtering + adapter)."""
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        f = json.load(fh)
    before = f.get("summary", {}).get("before_filtering", {})
    after = f.get("summary", {}).get("after_filtering", {})
    adapter = f.get("adapter_cutting", {})
    return {
        "fastp_total_reads": before.get("total_reads", 0),
        "fastp_total_bases": before.get("total_bases", 0),
        "fastp_q20_rate": before.get("q20_rate", 0.0),
        "fastp_q30_rate": before.get("q30_rate", 0.0),
        "fastp_gc_content": before.get("gc_content", 0.0),
        "fastp_read1_mean_length": before.get("read1_mean_length", 0),
        "fastp_read2_mean_length": before.get("read2_mean_length", 0),
        "fastp_reads_after_filter": after.get("total_reads", 0),
        "fastp_q30_rate_after": after.get("q30_rate", 0.0),
        "fastp_adapter_trimmed_reads": adapter.get("adapter_trimmed_reads", 0),
        "fastp_adapter_trimmed_bases": adapter.get("adapter_trimmed_bases", 0),
        "fastp_duplication_rate": f.get("duplication", {}).get("rate", 0.0),
    }


def parse_dupblaster(path):
    """Parse dupblaster's `--stats` per-library TSV (Fulcrum Genomics dupblaster,
    MIT; Rust successor to samblaster / Picard MarkDuplicates by the original
    Picard MarkDuplicates author).

    Columns used: mapped_pairs, duplicate_pairs, frac_duplicates,
    estimated_library_size. Rows are one-per-library; for the single-library QC
    input we sum pair counts across data rows (ignoring an empty "Unknown Library"
    row) and take a read-weighted duplication fraction.

    NOTE: dupblaster marks duplicates but does NOT parse flow-cell coordinates, so
    it cannot split optical vs PCR duplicates. The optical subset is recomputed
    downstream by bskryb-qc (Picard/Sentieon OpticalDuplicateFinder, pixel distance
    100) and merged into dedup_read_pair_optical_duplicates / pct_optical_duplicates /
    pct_pcr_duplicates by the caller. These placeholders (0) are the fallback when no
    rustqc optical block is present."""
    out = {
        "dedup_read_pairs_examined": 0,
        "dedup_read_pair_duplicates": 0,
        "dedup_read_pair_optical_duplicates": 0,  # dupblaster: no optical detection
        "dedup_estimated_library_size": 0,
        "dedup_pct_duplication": 0.0,
        "pct_duplication": 0.0,
        "pct_optical_duplicates": 0.0,            # dupblaster: no optical detection
        "pct_pcr_duplicates": 0.0,
    }
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    if len(rows) < 2:
        return out
    header = rows[0]
    idx = {name: i for i, name in enumerate(header)}

    def col(row, name, cast=float, default=0):
        i = idx.get(name)
        if i is None or i >= len(row) or row[i] == "":
            return default
        try:
            return cast(row[i])
        except ValueError:
            return default

    # Sum across all library data rows (single-library input has exactly one).
    mapped_pairs = 0
    dup_pairs = 0
    est_lib = 0
    frac_num = 0.0   # read-weighted numerator for overall frac_duplicates
    frac_den = 0.0
    for row in rows[1:]:
        mp = col(row, "mapped_pairs", int)
        dp = col(row, "duplicate_pairs", int)
        fd = col(row, "frac_duplicates", float)
        mapped_pairs += mp
        dup_pairs += dp
        est_lib += col(row, "estimated_library_size", int)
        # weight each library's fraction by its paired-read count (2 reads/pair)
        frac_num += fd * (2 * mp)
        frac_den += (2 * mp)

    frac_dup = (frac_num / frac_den) if frac_den else 0.0

    # Single-end input (Ultima CRAM path): no pairs, every mapped read is an "orphan".
    # Picard reports these as UNPAIRED_READS_EXAMINED / UNPAIRED_READ_DUPLICATES and
    # PERCENT_DUPLICATION = dups / examined. Paired runs never take this branch, so their
    # numbers are unchanged.
    if mapped_pairs == 0:
        orphans = sum(col(row, "mapped_orphans", int) for row in rows[1:])
        orphan_dups = sum(col(row, "duplicate_orphans", int) for row in rows[1:])
        if orphans > 0:
            frac_dup = orphan_dups / orphans
            out["dedup_unpaired_reads_examined"] = orphans
            out["dedup_unpaired_read_duplicates"] = orphan_dups

    out["dedup_read_pairs_examined"] = mapped_pairs
    out["dedup_read_pair_duplicates"] = dup_pairs
    out["dedup_estimated_library_size"] = est_lib
    out["dedup_pct_duplication"] = frac_dup
    out["pct_duplication"] = frac_dup
    out["pct_pcr_duplicates"] = frac_dup   # all detected dups are PCR (no optical split)
    return out


def parse_read_counts(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    if len(lines) >= 2:
        total, final = int(lines[0]), int(lines[1])
        return {"total_reads": total, "final_reads": final, "subsampled": total != final}
    return {}


def parse_preseq(sample):
    """Library complexity from the actual preseq tool ({sample}_preseq.txt, one line =
    gc_extrap last point). Preferred over bskryb-qc's port, which reflects the minibwa
    read set rather than the Sentieon rmdup baseline. Returns int or None."""
    for f in glob.glob(f"{sample}_preseq.txt") or glob.glob(f"{sample}*_preseq.txt"):
        try:
            with open(f) as fh:
                lines = [l.strip() for l in fh if l.strip()]
            if lines:
                return int(float(lines[0]))
        except (ValueError, OSError):
            pass
    return None


def alignment_category(m):
    """AlignmentStat row used for the align_* columns: PAIR for paired-end data, falling
    back to UNPAIRED for single-end (Ultima) input - the same PAIR -> UNPAIRED preference
    basej-dnaqc's QC_PLOTS applies to Sentieon/Picard AlignmentStat. bskryb-qc >= 0.1.4
    emits UNPAIRED only when single-end reads are present."""
    stats = m.get("alignment_stat", {})
    pair = stats.get("PAIR", {})
    if pair.get("TOTAL_READS", 0):
        return pair
    unpaired = stats.get("UNPAIRED", {})
    return unpaired if unpaired.get("TOTAL_READS", 0) else pair


def rustqc_fields(m):
    """Map the bskryb-qc rustqc JSON dict to the align_/insert_/gc_/cov_/sentieon_/
    lorenz/preseq columns (same key->column mapping as basej-dnaqc's QC_PLOTS)."""
    pair = alignment_category(m)
    # PCT_CHIMERAS component breakdown (bskryb-qc CHIMERAS_BREAKDOWN); a read may
    # satisfy several categories so counts can sum above the union total.
    chim = pair.get("CHIMERAS_BREAKDOWN", {})
    ins = m.get("insert_size", {})
    # Single-end input (Ultima CRAM) has no insert-size distribution: report null, not 0
    # (same as basej-wgsqc and the Picard path of basej-dnaqc).
    single_end = not ins.get("READ_PAIRS")
    contig = m.get("contig", {})
    qy = m.get("quality_yield", {})
    gc = m.get("gcbias", {})
    ps = m.get("preseq", {})
    lz = m.get("lorenz", {})

    total = pair.get("TOTAL_READS", 0)
    aligned = pair.get("PF_READS_ALIGNED", 0)

    return {
        # alignment (from Rust single pass)
        "align_total_reads": total,
        "align_pf_reads": total,               # PF == total (no vendor-fail reads)
        "align_pf_reads_aligned": aligned,
        "align_pf_hq_aligned_reads": pair.get("PF_HQ_ALIGNED_READS", 0),
        "align_pf_hq_aligned_bases": pair.get("PF_HQ_ALIGNED_BASES", 0),
        "align_pf_hq_aligned_q20_bases": pair.get("PF_HQ_ALIGNED_Q20_BASES", 0),
        "align_pf_mismatch_rate": pair.get("PF_MISMATCH_RATE", 0.0),
        "align_pf_hq_error_rate": pair.get("PF_HQ_ERROR_RATE", 0.0),
        "align_pf_indel_rate": pair.get("PF_INDEL_RATE", 0.0),
        "align_mean_read_length": pair.get("MEAN_READ_LENGTH", 0.0),
        "align_reads_aligned_in_pairs": pair.get("READS_ALIGNED_IN_PAIRS", 0),
        "align_pct_reads_aligned_in_pairs": pair.get("PCT_READS_ALIGNED_IN_PAIRS", 0.0),
        "align_strand_balance": pair.get("STRAND_BALANCE", 0.0),
        "align_pct_chimeras": pair.get("PCT_CHIMERAS", 0.0),
        # PCT_CHIMERAS component breakdown (counts + fraction of pairs)
        "align_chim_diff_contig": chim.get("different_contig", 0),
        "align_chim_diff_contig_pct": chim.get("different_contig_pct", 0.0),
        "align_chim_bad_orient": chim.get("bad_orientation", 0),
        "align_chim_bad_orient_pct": chim.get("bad_orientation_pct", 0.0),
        "align_chim_large_insert": chim.get("large_insert_gt_100kb", 0),
        "align_chim_large_insert_pct": chim.get("large_insert_gt_100kb_pct", 0.0),
        "align_chim_split_sa": chim.get("split_read_sa_tag", 0),
        "align_chim_split_sa_pct": chim.get("split_read_sa_tag_pct", 0.0),
        # derived convenience columns (match basej-dnaqc names)
        "pct_aligned": (aligned / total) if total else 0.0,
        "pct_pf": 1.0 if total else 0.0,
        "pct_error": pair.get("PF_MISMATCH_RATE", 0.0),
        "pct_chimeras": pair.get("PCT_CHIMERAS", 0.0),
        # insert size
        "insert_median": None if single_end else ins.get("MEDIAN_INSERT_SIZE", 0.0),
        "insert_size": None if single_end else ins.get("MEDIAN_INSERT_SIZE", 0.0),
        "insert_mean": None if single_end else ins.get("MEAN_INSERT_SIZE", 0.0),
        "insert_std": None if single_end else ins.get("STANDARD_DEVIATION", 0.0),
        "insert_read_pairs": ins.get("READ_PAIRS", 0),
        "insert_median_absolute_deviation": None if single_end else ins.get("MEDIAN_ABSOLUTE_DEVIATION"),
        "insert_min": None if single_end else ins.get("MIN_INSERT_SIZE"),
        "insert_max": None if single_end else ins.get("MAX_INSERT_SIZE"),
        # coverage / mitochondrial fraction (CoverageMetrics-equivalent)
        "cov_total_bases": contig.get("cov_total_bases", 0),
        "cov_chrm_bases": contig.get("cov_chrm_bases", 0),
        "pct_chrm": contig.get("pct_chrm", 0.0),
        # QualityYield (Q20/Q30 over all primary reads)
        "sentieon_total_bases": qy.get("TOTAL_BASES", 0),
        "sentieon_q20_rate": qy.get("PCT_Q20", 0.0),
        "sentieon_q30_rate": qy.get("PCT_Q30", 0.0),
        # GC bias (present only when --reference was provided to bskryb-qc)
        "gc_total_clusters": gc.get("TOTAL_CLUSTERS"),
        "gc_aligned_reads": gc.get("ALIGNED_READS"),
        "gc_at_dropout": gc.get("AT_DROPOUT"),
        "gc_gc_dropout": gc.get("GC_DROPOUT"),
        "gc_gc_nc_0_19": gc.get("GC_NC_0_19"),
        "gc_gc_nc_20_39": gc.get("GC_NC_20_39"),
        "gc_gc_nc_40_59": gc.get("GC_NC_40_59"),
        "gc_gc_nc_60_79": gc.get("GC_NC_60_79"),
        "gc_gc_nc_80_100": gc.get("GC_NC_80_100"),
        # preseq library complexity (gc_extrap bootstrap-median reimplementation)
        "preseq_count": ps.get("preseq_count"),
        # Lorenz-curve coverage evenness / Gini (bam-lorenz-coverage reproduction)
        "gini_coefficient_index": lz.get("gini_coefficient_index"),
        "roc_lorenz_curve": lz.get("roc_lorenz_curve"),
        "total_covered_positions_of_genome": lz.get("total_covered_positions_of_genome"),
        "total_sequenced_bases": lz.get("total_sequenced_bases"),
        "total_investigated_genomic_positions": lz.get("total_investigated_genomic_positions"),
    }


def normalize_sample(sample):
    """Strip one terminal Ginkgo `_sorted` suffix from a sample ID."""
    sample = str(sample)
    return sample[:-len("_sorted")] if sample.endswith("_sorted") else sample


def load_ginkgo(metrics_path, segcopy_path):
    """Return ploidy dictionaries from the initial Ginkgo metrics and SegCopy."""
    ploidy, avg_ploidy = {}, {}
    if metrics_path and os.path.exists(metrics_path):
        try:
            df = pd.read_csv(metrics_path, sep="\t")
            if "SampleId" in df.columns:
                df = df.rename(columns={"SampleId": "biosample"})
            df["biosample"] = df["biosample"].map(normalize_sample)
            for _, row in df.iterrows():
                if pd.notna(row.get("GenomePloidy")):
                    ploidy[row["biosample"]] = str(row.get("GenomePloidy"))
        except Exception as e:
            print(f"Warning: could not load Ginkgo metrics: {e}")
    if segcopy_path and os.path.exists(segcopy_path):
        try:
            seg = pd.read_csv(segcopy_path, sep="\t")
            for col in seg.columns:
                if col in ("CHR", "START", "END"):
                    continue
                avg_ploidy[normalize_sample(col)] = float(seg[col].mean())
        except Exception as e:
            print(f"Warning: could not load SegCopy: {e}")
    return ploidy, avg_ploidy


def load_cnv_summary(cnv_summary_path):
    """Return per-sample MAPD_CNV_Log2 and SKEW_CNV dictionaries."""
    if not cnv_summary_path:
        return {}, {}
    if not os.path.exists(cnv_summary_path):
        raise SystemExit(
            f"ERROR: CNV summary '{cnv_summary_path}' does not exist; "
            "check the --cnv-summary path."
        )
    try:
        df = pd.read_csv(cnv_summary_path, sep="\t")
    except pd.errors.EmptyDataError:
        raise SystemExit(
            f"ERROR: CNV summary '{cnv_summary_path}' is empty; regenerate it with sample rows."
        )
    except Exception as e:
        raise SystemExit(
            f"ERROR: could not read CNV summary '{cnv_summary_path}': {e}; "
            "check its permissions and TSV format."
        ) from e
    if df.empty:
        raise SystemExit(
            f"ERROR: CNV summary '{cnv_summary_path}' is empty; regenerate it with sample rows."
        )
    required = {"SampleId", "MAPD_CNV_Log2", "SKEW_CNV"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(
            f"ERROR: CNV summary '{cnv_summary_path}' is missing required column(s): "
            f"{', '.join(missing)}."
        )

    mapd, skew = {}, {}
    for _, row in df.iterrows():
        sample = normalize_sample(row["SampleId"])
        if sample in mapd:
            raise SystemExit(
                f"ERROR: CNV summary '{cnv_summary_path}' has duplicate sample ID "
                f"'{sample}' after normalization; make sample IDs unique."
            )
        mapd_value = pd.to_numeric(row["MAPD_CNV_Log2"], errors="coerce")
        skew_value = pd.to_numeric(row["SKEW_CNV"], errors="coerce")
        mapd[sample] = float(mapd_value) if pd.notna(mapd_value) else None
        skew[sample] = float(skew_value) if pd.notna(skew_value) else None
    return mapd, skew


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset-id", required=True)
    ap.add_argument("--workspace", required=True)
    ap.add_argument("--workflow-id", required=True)
    ap.add_argument("--pipeline-version", required=True)
    ap.add_argument("--user", required=True)
    ap.add_argument("--ginkgo-metrics", default=None)
    ap.add_argument("--segcopy", default=None)
    ap.add_argument("--cnv-summary", default=None)
    args = ap.parse_args()

    # Discover samples from the rustqc JSONs staged in the work dir.
    samples = sorted(os.path.basename(p)[:-len(".rustqc.json")]
                     for p in glob.glob("*.rustqc.json"))
    if not samples:
        raise SystemExit("ERROR: no *.rustqc.json files found")
    print(f"Found {len(samples)} samples: {samples}")

    ploidy_by_sample, avg_ploidy_by_sample = load_ginkgo(args.ginkgo_metrics, args.segcopy)
    mapd_by_sample, skew_by_sample = load_cnv_summary(args.cnv_summary)
    if args.cnv_summary:
        missing = [sample for sample in samples if sample not in mapd_by_sample]
        if missing:
            raise SystemExit(
                "ERROR: CNV summary is missing expected sample(s): "
                f"{', '.join(missing)}; regenerate it for every Rust QC sample."
            )
        cnv_rows = [
            {"ginkgo_cnv_mapd": mapd_by_sample[sample],
             "ginkgo_cnv_skew": skew_by_sample[sample]}
            for sample in samples
        ]
        qc_scoring.require_metric(
            cnv_rows, ("ginkgo_cnv_mapd", "ginkgo_cnv_skew"), "CNV summary"
        )

    all_summaries = []
    for sample in samples:
        summary = {
            "biosample": sample,
            "dataset_id": args.dataset_id,
            "pipeline": "bj-dnaqc",
            "pipeline_version": args.pipeline_version,
            "molecule_type": "dna",
            "user": args.user,
            "workflow_id": args.workflow_id,
            "workspace": args.workspace,
        }

        with open(f"{sample}.rustqc.json") as fh:
            m = json.load(fh)
        summary.update(rustqc_fields(m))
        summary.update(chimera_class_fields(m))

        # Merge in fastp, dupblaster (dedup), and read-count fields
        summary.update(parse_fastp(f"{sample}_fastp.json"))
        summary.update(parse_dupblaster(f"{sample}.dupblaster.tsv"))
        summary.update(parse_read_counts(f"{sample}_read_counts.txt"))

        # No FASTQ on the Ultima CRAM path, so no fastp JSON: report Q20/Q30 from the
        # rustqc QualityYield instead (same fallback as basej-dnaqc's QC_PLOTS).
        if "fastp_q30_rate" not in summary:
            summary["fastp_q20_rate"] = summary.get("sentieon_q20_rate")
            summary["fastp_q30_rate"] = summary.get("sentieon_q30_rate")

        # Override bskryb-qc preseq_count with the actual preseq tool when staged.
        _preseq = parse_preseq(sample)
        if _preseq is not None:
            summary["preseq_count"] = _preseq

        # Optical-duplicate split. dupblaster marks duplicates but does not parse
        # flow-cell coordinates, so bskryb-qc recomputes the optical subset from read
        # names (Picard/Sentieon OpticalDuplicateFinder, pixel distance 100). Split the
        # total duplication into optical vs PCR to match Sentieon Dedup's
        # READ_PAIR_OPTICAL_DUPLICATES.
        optical = int(m.get("optical", {}).get("read_pair_optical_duplicates", 0))
        read_pairs_examined = summary.get("dedup_read_pairs_examined", 0) or 0
        frac_dup = summary.get("dedup_pct_duplication", 0.0) or 0.0
        if read_pairs_examined > 0:
            pct_optical = optical / read_pairs_examined
            summary["dedup_read_pair_optical_duplicates"] = optical
            summary["pct_optical_duplicates"] = pct_optical
            # remaining duplication is PCR (clamp at 0 for safety)
            summary["pct_pcr_duplicates"] = max(frac_dup - pct_optical, 0.0)

        # Ginkgo: ploidy from the initial metrics/SegCopy and evenness metrics
        # from the R cnvSummarizer output.
        summary["ginkgo_cnv_mapd"] = mapd_by_sample.get(sample)
        summary["ginkgo_cnv_skew"] = skew_by_sample.get(sample)
        summary["ginkgo_ploidy"] = ploidy_by_sample.get(sample)
        summary["ginkgo_average_ploidy"] = avg_ploidy_by_sample.get(sample)

        # ---- QC tier ----
        #
        # Same rules basej-dnaqc applies via dna_qc_plot.R, from qc_scoring.py, which is
        # validated against the R library (15000/15000 on randomised inputs, see
        # tests/test_qc_scoring.py). Applied here because this pipeline writes Parquet
        # straight from the Rust single pass and never invokes R for final aggregation.
        verdict = qc_scoring.assign_tier(
            summary,
            qc_scoring.spec_dna(),
            below_floor=qc_scoring.depth_floor(summary.get("total_reads"),
                                               min_value=DEPTH_FLOOR_READS),
            required=("ginkgo_cnv_mapd", "ginkgo_cnv_skew"),
        )
        summary["qc_score"] = verdict["qc_score"]
        summary["qc_status"] = verdict["qc_status"]
        summary["qc_label"] = verdict["qc_label"]
        summary["blocking_metric"] = verdict["blocking_metric"]

        # A zero-read cell is an outright FAIL with no score, matching basej-dnaqc. The tier
        # scoring would call it Inconclusive (below the depth floor), which is defensible, but
        # both pipelines feed the same Athena table so they must agree.
        if int(qc_scoring.qc_numeric(summary.get("total_reads")) or 0) == 0:
            summary["qc_status"] = "FAIL"
            summary["qc_score"] = None
            summary["qc_label"] = None
            summary["blocking_metric"] = "no_reads"

        df = pd.DataFrame([summary])
        for c in _double_cols:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors='coerce').astype('float64')
        for c in _bigint_cols:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors='coerce').astype('Int64')
        for c in _bool_cols:
            if c in df.columns:
                df[c] = df[c].astype('boolean')
        for c in _str_cols:
            if c in df.columns:
                df[c] = df[c].astype('string')

        out_dir = (f"dnaqc_summary/workspace={args.workspace}"
                   f"/workflow_id={args.workflow_id}/biosample={sample}")
        os.makedirs(out_dir, exist_ok=True)
        pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                       os.path.join(out_dir, "output.parquet"))
        all_summaries.append(summary)

        # MultiQC custom-content table (subset of columns). The chimera breakdown
        # pct columns sit next to PCT_CHIMERAS (matches basej-methylqc-rs).
        # Higher-level chimera topology classification (chim_class_*) is appended at
        # the END of the table so the existing column order stays stable. One column
        # per named mechanism class from chimera_classes.yaml (CHIM_CLASSES order).
        chim_class_titles = {
            "proper": "ChimCls_Proper_pct",
            "interchromosomal": "ChimCls_Interchrom_pct",
            "local_inverted_hairpin": "ChimCls_LocalInvHairpin_pct",
            "inverted_nonlocal": "ChimCls_InvNonlocal_pct",
            "everted_tandem": "ChimCls_EvertedTandem_pct",
            "large_insert": "ChimCls_LargeInsertCls_pct",
            "split_read_only": "ChimCls_SplitReadOnly_pct",
            "other": "ChimCls_Other_pct",
        }
        chim_class_cols = [chim_class_titles[k] for k in CHIM_CLASSES]
        cols = ["sample_name", "QC_Status", "Score", "preseq_count",
                "PCT_CHIMERAS", "Chim_DiffContig_pct", "Chim_BadOrient_pct",
                "Chim_LargeInsert_pct", "Chim_SplitSA_pct", "chrM", "MAPD_CNV",
                "SKEW_CNV", "TotalReads", "FinalReads", "insert_size", "pct_duplication",
                "fastp_q30_rate", "pct_aligned", "sentieon_q30_rate"] + chim_class_cols
        # Raw numeric values; display precision is controlled by the per-column
        # `headers.format` block below so small rate metrics are not rounded to
        # 0.00 by MultiQC's default {:,.1f} formatting.
        def _num(v):
            return "" if v is None else v
        vals = [sample, summary.get("qc_status") or "NA",
                _num(summary.get("qc_score"))]
        vals.append(_num(summary.get("preseq_count")))
        vals.append(_num(summary.get('align_pct_chimeras')))
        vals.append(_num(summary.get('align_chim_diff_contig_pct')))
        vals.append(_num(summary.get('align_chim_bad_orient_pct')))
        vals.append(_num(summary.get('align_chim_large_insert_pct')))
        vals.append(_num(summary.get('align_chim_split_sa_pct')))
        vals.append(_num(summary.get('pct_chrm')))
        vals.append(_num(summary.get('ginkgo_cnv_mapd')))
        vals.append(_num(summary.get('ginkgo_cnv_skew')))
        vals.append(_num(summary.get('total_reads')))
        vals.append(_num(summary.get('final_reads')))
        vals.append(_num(summary.get('insert_size')))
        vals.append(_num(summary.get('dedup_pct_duplication')))
        vals.append(_num(summary.get('fastp_q30_rate')))
        vals.append(_num(summary.get('pct_aligned')))
        vals.append(_num(summary.get('sentieon_q30_rate')))
        for k in CHIM_CLASSES:
            vals.append(_num(summary.get(f"chim_class_{k}_pct")))
        vals = [str(v) for v in vals]

        _fmt = {
            "sample_name": None, "QC_Status": None, "Score": "{:.2f}",
            "preseq_count": "{:,.0f}", "PCT_CHIMERAS": "{:.4f}",
            "Chim_DiffContig_pct": "{:.4f}", "Chim_BadOrient_pct": "{:.4f}",
            "Chim_LargeInsert_pct": "{:.4f}", "Chim_SplitSA_pct": "{:.4f}",
            "chrM": "{:.4f}", "MAPD_CNV": "{:.4f}", "SKEW_CNV": "{:.4f}",
            "TotalReads": "{:,.0f}", "FinalReads": "{:,.0f}",
            "insert_size": "{:.0f}", "pct_duplication": "{:.4f}",
            "fastp_q30_rate": "{:.4f}", "pct_aligned": "{:.4f}",
            "sentieon_q30_rate": "{:.4f}",
        }
        _fmt.update({c: "{:.4f}" for c in chim_class_cols})
        with open(f"{sample}_selected_metrics_mqc.txt", "w") as fh:
            fh.write("# id: 'dnaqc_summary'\n")
            fh.write("# plot_type: 'table'\n")
            fh.write("# section_name: 'QC Summary'\n")
            fh.write("# description: 'Per-sample alignment, chimera breakdown, coverage, and library complexity QC metrics.'\n")
            fh.write("# pconfig:\n")
            fh.write("#   id: 'dnaqc_summary_table'\n")
            fh.write("# headers:\n")
            for i, c in enumerate(cols):
                fh.write(f"#   {c}:\n")
                fh.write(f"#     title: '{c}'\n")
                fh.write(f"#     placement: {(i + 1) * 10}\n")
                fmt = _fmt.get(c)
                if fmt:
                    fh.write(f"#     format: '{fmt}'\n")
            fh.write("\t".join(cols) + "\n")
            fh.write("\t".join(vals) + "\n")
        print(f"  {sample}: dnaqc_summary parquet written")

    all_df = pd.DataFrame(all_summaries)
    all_df.to_csv("dnaqc_all_metrics.tsv", sep="\t", index=False)
    # Per-biosample QC verdict index (same columns as basej-dnaqc's QC_PLOTS output).
    all_df[["biosample", "qc_status", "pipeline", "pipeline_version"]].rename(
        columns={"biosample": "biosampleName"}
    ).to_csv("per_biosample_status.csv", index=False)
    print(f"Wrote {len(samples)} parquets to dnaqc_summary/")


if __name__ == "__main__":
    main()
