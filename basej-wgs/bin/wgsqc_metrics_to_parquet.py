#!/usr/bin/env python3
"""Convert bskryb-qc and dupblaster metrics to canonical basej-wgs outputs.

The per-biosample Parquet uses the exact 97-field order and Arrow types emitted by
basej-wgs in WGS mode. Target-capture fields remain typed nulls for WGS. Diagnostic
scoring details are logged but are intentionally not added to the shared Iceberg table.
"""
import argparse
import json
import os

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import het_sensitivity
import qc_scoring

DEPTH_FLOOR_FRAC_OF_TARGET = 0.1

S = pa.string()
I = pa.int64()
F = pa.float64()

# Physical contract observed in the canonical basej-wgs connector Parquet. Keep names,
# order, and types stable: the lakehouse table is shared by Sentieon and Rust producers.
WGSQC_SCHEMA = pa.schema([
    pa.field("biosample", S),
    pa.field("dataset_id", S),
    pa.field("pipeline", S),
    pa.field("pipeline_version", S),
    pa.field("molecule_type", S),
    pa.field("mode", S),
    pa.field("genome", S),
    pa.field("workspace", S),
    pa.field("workflow_id", S),
    pa.field("user", S),
    pa.field("qc_status", S),
    pa.field("align_total_reads", I),
    pa.field("pf_reads_aligned", I),
    pa.field("pf_aligned_bases", I),
    pa.field("pf_hq_aligned_reads", I),
    pa.field("pf_hq_aligned_bases", I),
    pa.field("pf_hq_aligned_q20_bases", I),
    pa.field("pf_mismatch_rate", F),
    pa.field("pf_hq_error_rate", F),
    pa.field("pf_indel_rate", F),
    pa.field("mean_read_length", F),
    pa.field("pct_reads_aligned_in_pairs", F),
    pa.field("pct_chimeras", F),
    pa.field("pct_adapter", F),
    pa.field("strand_balance", F),
    pa.field("sentieon_total_bases", I),
    pa.field("sentieon_q20_rate", F),
    pa.field("sentieon_q30_rate", F),
    pa.field("pct_duplication", F),
    pa.field("estimated_library_size", I),
    pa.field("read_pairs_examined", I),
    pa.field("read_pair_duplicates", I),
    pa.field("read_pair_optical_duplicates", I),
    pa.field("insert_median", F),
    pa.field("insert_mad", F),
    pa.field("insert_min", I),
    pa.field("insert_max", I),
    pa.field("insert_mean", F),
    pa.field("insert_std", F),
    pa.field("insert_read_pairs", I),
    pa.field("at_dropout", F),
    pa.field("gc_dropout", F),
    pa.field("gc_nc_0_19", F),
    pa.field("gc_nc_20_39", F),
    pa.field("gc_nc_40_59", F),
    pa.field("gc_nc_60_79", F),
    pa.field("gc_nc_80_100", F),
    pa.field("mean_coverage", F),
    pa.field("sd_coverage", F),
    pa.field("median_coverage", F),
    pa.field("mad_coverage", F),
    pa.field("pct_1x", F),
    pa.field("pct_5x", F),
    pa.field("pct_10x", F),
    pa.field("pct_15x", F),
    pa.field("pct_20x", F),
    pa.field("pct_25x", F),
    pa.field("pct_30x", F),
    pa.field("pct_40x", F),
    pa.field("pct_50x", F),
    pa.field("pct_60x", F),
    pa.field("pct_70x", F),
    pa.field("pct_80x", F),
    pa.field("pct_90x", F),
    pa.field("pct_100x", F),
    pa.field("pct_exc_mapq", F),
    pa.field("pct_exc_dupe", F),
    pa.field("pct_exc_unpaired", F),
    pa.field("pct_exc_baseq", F),
    pa.field("pct_exc_overlap", F),
    pa.field("pct_exc_capped", F),
    pa.field("pct_exc_total", F),
    pa.field("het_snp_sensitivity", F),
    pa.field("het_snp_q", I),
    pa.field("mean_target_coverage", F),
    pa.field("fold_enrichment", F),
    pa.field("pct_selected_bases", F),
    pa.field("fold_80_base_penalty", F),
    pa.field("zero_cvg_targets_pct", F),
    pa.field("on_bait_bases", I),
    pa.field("near_bait_bases", I),
    pa.field("off_bait_bases", I),
    pa.field("pct_target_bases_1x", F),
    pa.field("pct_target_bases_2x", F),
    pa.field("pct_target_bases_10x", F),
    pa.field("pct_target_bases_20x", F),
    pa.field("pct_target_bases_30x", F),
    pa.field("pct_target_bases_40x", F),
    pa.field("pct_target_bases_50x", F),
    pa.field("pct_target_bases_100x", F),
    pa.field("pct_target_bases_250x", F),
    pa.field("pct_target_bases_500x", F),
    pa.field("pct_target_bases_1000x", F),
    pa.field("total_reads", I),
    pa.field("total_reads_source", S),
    pa.field("final_reads", I),
    pa.field("qc_score", I),
])

REQUIRED_COVERAGE_KEYS = [
    "MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE", "MAD_COVERAGE",
    "PCT_1X", "PCT_5X", "PCT_10X", "PCT_15X", "PCT_20X", "PCT_25X",
    "PCT_30X", "PCT_40X", "PCT_50X", "PCT_60X", "PCT_70X", "PCT_80X",
    "PCT_90X", "PCT_100X", "PCT_EXC_MAPQ", "PCT_EXC_DUPE",
    "PCT_EXC_UNPAIRED", "PCT_EXC_BASEQ", "PCT_EXC_OVERLAP",
    "PCT_EXC_CAPPED", "PCT_EXC_TOTAL",
]

TARGET_NULL_FIELDS = [
    "mean_target_coverage", "fold_enrichment", "pct_selected_bases",
    "fold_80_base_penalty", "zero_cvg_targets_pct", "on_bait_bases",
    "near_bait_bases", "off_bait_bases", "pct_target_bases_1x",
    "pct_target_bases_2x", "pct_target_bases_10x", "pct_target_bases_20x",
    "pct_target_bases_30x", "pct_target_bases_40x", "pct_target_bases_50x",
    "pct_target_bases_100x", "pct_target_bases_250x",
    "pct_target_bases_500x", "pct_target_bases_1000x",
]


def parse_dupblaster(path):
    """Parse dupblaster's stats TSV into canonical dedup fields."""
    if not path:
        return {}
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if len(rows) < 2:
        return {}
    index = {name: i for i, name in enumerate(rows[0])}

    def value(row, name, cast=float, default=0):
        i = index.get(name)
        if i is None or i >= len(row) or row[i] == "":
            return default
        try:
            return cast(row[i])
        except ValueError:
            return default

    mapped_pairs = duplicate_pairs = estimated_library_size = 0
    fraction_numerator = fraction_denominator = 0.0
    for row in rows[1:]:
        pairs = value(row, "mapped_pairs", int)
        mapped_pairs += pairs
        duplicate_pairs += value(row, "duplicate_pairs", int)
        estimated_library_size += value(row, "estimated_library_size", int)
        fraction_numerator += value(row, "frac_duplicates", float) * (2 * pairs)
        fraction_denominator += 2 * pairs

    return {
        "read_pairs_examined": mapped_pairs,
        "read_pair_duplicates": duplicate_pairs,
        "read_pair_optical_duplicates": 0,
        "estimated_library_size": estimated_library_size,
        "pct_duplication": (
            fraction_numerator / fraction_denominator if fraction_denominator else 0.0
        ),
    }


HS_FLOAT_FIELDS = {
    "mean_target_coverage": "MEAN_TARGET_COVERAGE",
    "fold_enrichment": "FOLD_ENRICHMENT",
    "pct_selected_bases": "PCT_SELECTED_BASES",
    "zero_cvg_targets_pct": "ZERO_CVG_TARGETS_PCT",
    "pct_target_bases_1x": "PCT_TARGET_BASES_1X",
    "pct_target_bases_2x": "PCT_TARGET_BASES_2X",
    "pct_target_bases_10x": "PCT_TARGET_BASES_10X",
    "pct_target_bases_20x": "PCT_TARGET_BASES_20X",
    "pct_target_bases_30x": "PCT_TARGET_BASES_30X",
    "pct_target_bases_40x": "PCT_TARGET_BASES_40X",
    "pct_target_bases_50x": "PCT_TARGET_BASES_50X",
    "pct_target_bases_100x": "PCT_TARGET_BASES_100X",
    "pct_target_bases_250x": "PCT_TARGET_BASES_250X",
    "pct_target_bases_500x": "PCT_TARGET_BASES_500X",
    "pct_target_bases_1000x": "PCT_TARGET_BASES_1000X",
}
HS_INT_FIELDS = {
    "on_bait_bases": "ON_BAIT_BASES",
    "near_bait_bases": "NEAR_BAIT_BASES",
    "off_bait_bases": "OFF_BAIT_BASES",
}


def parse_hsmetrics(path):
    """Target-capture fields from Picard CollectHsMetrics (exome mode).

    Same field mapping as basej-wgs's WGS_QC_METRICS_TO_PARQUET. Picard writes '?' for an
    undefined value (e.g. FOLD_80_BASE_PENALTY with zero target coverage), mapped to null.
    Missing/empty file -> {} so the caller can fail loudly."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    with open(path) as fh:
        lines = [line.rstrip("\n") for line in fh]
    for i, line in enumerate(lines):
        if line.startswith("BAIT_SET") or ("MEAN_TARGET_COVERAGE" in line and not line.startswith("#")):
            if i + 1 >= len(lines) or not lines[i + 1].strip():
                return {}
            row = dict(zip(line.split("\t"), lines[i + 1].split("\t")))

            def num(key, cast):
                raw = row.get(key)
                if raw in (None, "", "?"):
                    return None
                try:
                    return cast(float(raw)) if cast is int else cast(raw)
                except ValueError:
                    return None

            out = {name: num(key, float) for name, key in HS_FLOAT_FIELDS.items()}
            out.update({name: num(key, int) for name, key in HS_INT_FIELDS.items()})
            out["fold_80_base_penalty"] = num("FOLD_80_BASE_PENALTY", float)
            return out
    return {}


def parse_flagstat(path):
    """Duplication from the existing 0x400 flags via `samtools flagstat` (CRAM path).

    Pre-aligned Ultima CRAMs arrive with duplicates already marked by the vendor demux, so
    no dupblaster stats exist. Same derivation as basej-wgs's flagstat fallback:
    PERCENT_DUPLICATION = primary duplicates / primary mapped (QC-pass counts)."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    primary = primary_dup = primary_mapped = 0
    with open(path) as fh:
        for line in fh:
            text = line.strip()
            tokens = text.split()
            if len(tokens) < 4:
                continue
            try:
                count = int(tokens[0])
            except ValueError:
                continue
            descriptor = text.split(None, 3)[-1]
            if descriptor.startswith("primary mapped"):
                primary_mapped = count
            elif descriptor.startswith("primary duplicates"):
                primary_dup = count
            elif descriptor == "primary":
                primary = count
    denominator = primary_mapped if primary_mapped > 0 else primary
    return {
        "read_pairs_examined": denominator,
        "read_pair_duplicates": primary_dup,
        "read_pair_optical_duplicates": 0,
        "pct_duplication": (primary_dup / denominator) if denominator else 0.0,
    }


def alignment_category(metrics):
    """PAIR for paired-end data; UNPAIRED (bskryb-qc >= 0.1.4) for single-end input."""
    stats = metrics["alignment_stat"]
    pair = stats["PAIR"]
    if pair.get("TOTAL_READS", 0):
        return pair
    unpaired = stats.get("UNPAIRED", {})
    return unpaired if unpaired.get("TOTAL_READS", 0) else pair


def parse_read_counts(path, align_total_reads):
    """Preserve raw, final, and alignment-derived read-count semantics."""
    if path:
        with open(path) as fh:
            lines = [line.strip() for line in fh if line.strip()]
        if lines:
            return {
                "total_reads": int(lines[0]),
                "total_reads_source": "raw_input",
                "final_reads": int(lines[1]) if len(lines) > 1 else None,
            }
    return {
        "total_reads": int(align_total_reads),
        "total_reads_source": "alignment",
        "final_reads": None,
    }


def require_coverage(coverage, sample):
    missing = [key for key in REQUIRED_COVERAGE_KEYS if coverage.get(key) is None]
    if missing:
        raise SystemExit(
            f"ERROR: {sample}: missing required WGS coverage metrics {missing}. "
            "The run is not scorable; verify --genome-territory and the bskryb-qc image."
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True)
    parser.add_argument("--dupblaster-stats", default=None)
    parser.add_argument("--flagstat", default=None,
                        help="samtools flagstat of a pre-aligned CRAM (used when no dupblaster stats)")
    parser.add_argument("--read-counts", default=None)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", default="wgs", choices=["wgs", "exome"])
    parser.add_argument("--hsmetrics", default=None,
                        help="Picard CollectHsMetrics output (required for --mode exome)")
    parser.add_argument("--genome", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--workflow-id", required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--user", required=True)
    args = parser.parse_args()

    with open(args.json) as fh:
        metrics = json.load(fh)

    pair = alignment_category(metrics)
    insert = metrics.get("insert_size", {})
    if not insert.get("READ_PAIRS"):
        # Single-end input (Ultima CRAM) has no insert-size distribution; report the
        # insert_* fields as null rather than 0, as basej-wgs does for its CRAM path.
        insert = {}
    quality = metrics.get("quality_yield", {})
    gc = metrics.get("gcbias", {})
    coverage = metrics.get("coverage", {})
    exome = args.mode == "exome"
    prefix = "wes" if exome else "wgs"
    if exome:
        # WGS-style coverage is not part of exome QC (basej-wgs skips WgsMetricsAlgo in
        # exome mode); the pct_Nx/mean coverage block stays null and target coverage
        # comes from CollectHsMetrics instead.
        coverage = {}
    else:
        require_coverage(coverage, args.sample)

    align_total_reads = int(pair["TOTAL_READS"])
    summary = {
        "biosample": args.sample,
        "dataset_id": args.dataset_id,
        "pipeline": "basej-wgs",
        "pipeline_version": args.pipeline_version,
        "molecule_type": "dna",
        "mode": args.mode,
        "genome": args.genome,
        "workspace": args.workspace,
        "workflow_id": args.workflow_id,
        "user": args.user,
        "qc_status": None,
        "align_total_reads": align_total_reads,
        "pf_reads_aligned": pair["PF_READS_ALIGNED"],
        "pf_aligned_bases": pair.get("PF_ALIGNED_BASES", 0),
        "pf_hq_aligned_reads": pair["PF_HQ_ALIGNED_READS"],
        "pf_hq_aligned_bases": pair["PF_HQ_ALIGNED_BASES"],
        "pf_hq_aligned_q20_bases": pair["PF_HQ_ALIGNED_Q20_BASES"],
        "pf_mismatch_rate": pair["PF_MISMATCH_RATE"],
        "pf_hq_error_rate": pair["PF_HQ_ERROR_RATE"],
        "pf_indel_rate": pair["PF_INDEL_RATE"],
        "mean_read_length": pair["MEAN_READ_LENGTH"],
        "pct_reads_aligned_in_pairs": pair["PCT_READS_ALIGNED_IN_PAIRS"],
        "pct_chimeras": pair["PCT_CHIMERAS"],
        "pct_adapter": pair.get("PCT_ADAPTER", 0.0),
        "strand_balance": pair["STRAND_BALANCE"],
        "sentieon_total_bases": quality.get("TOTAL_BASES", 0),
        "sentieon_q20_rate": quality.get("PCT_Q20", 0.0),
        "sentieon_q30_rate": quality.get("PCT_Q30", 0.0),
        "insert_median": insert.get("MEDIAN_INSERT_SIZE"),
        "insert_mad": insert.get("MEDIAN_ABSOLUTE_DEVIATION"),
        "insert_min": insert.get("MIN_INSERT_SIZE"),
        "insert_max": insert.get("MAX_INSERT_SIZE"),
        "insert_mean": insert.get("MEAN_INSERT_SIZE"),
        "insert_std": insert.get("STANDARD_DEVIATION"),
        "insert_read_pairs": insert.get("READ_PAIRS"),
        "at_dropout": gc.get("AT_DROPOUT"),
        "gc_dropout": gc.get("GC_DROPOUT"),
        "gc_nc_0_19": gc.get("GC_NC_0_19"),
        "gc_nc_20_39": gc.get("GC_NC_20_39"),
        "gc_nc_40_59": gc.get("GC_NC_40_59"),
        "gc_nc_60_79": gc.get("GC_NC_60_79"),
        "gc_nc_80_100": gc.get("GC_NC_80_100"),
        "mean_coverage": coverage.get("MEAN_COVERAGE"),
        "sd_coverage": coverage.get("SD_COVERAGE"),
        "median_coverage": coverage.get("MEDIAN_COVERAGE"),
        "mad_coverage": coverage.get("MAD_COVERAGE"),
        "pct_1x": coverage.get("PCT_1X"),
        "pct_5x": coverage.get("PCT_5X"),
        "pct_10x": coverage.get("PCT_10X"),
        "pct_15x": coverage.get("PCT_15X"),
        "pct_20x": coverage.get("PCT_20X"),
        "pct_25x": coverage.get("PCT_25X"),
        "pct_30x": coverage.get("PCT_30X"),
        "pct_40x": coverage.get("PCT_40X"),
        "pct_50x": coverage.get("PCT_50X"),
        "pct_60x": coverage.get("PCT_60X"),
        "pct_70x": coverage.get("PCT_70X"),
        "pct_80x": coverage.get("PCT_80X"),
        "pct_90x": coverage.get("PCT_90X"),
        "pct_100x": coverage.get("PCT_100X"),
        "pct_exc_mapq": coverage.get("PCT_EXC_MAPQ"),
        "pct_exc_dupe": coverage.get("PCT_EXC_DUPE"),
        "pct_exc_unpaired": coverage.get("PCT_EXC_UNPAIRED"),
        "pct_exc_baseq": coverage.get("PCT_EXC_BASEQ"),
        "pct_exc_overlap": coverage.get("PCT_EXC_OVERLAP"),
        "pct_exc_capped": coverage.get("PCT_EXC_CAPPED"),
        "pct_exc_total": coverage.get("PCT_EXC_TOTAL"),
        # Computed below from the depth/quality histograms; bskryb-qc has no theoretical
        # sensitivity model of its own, so these keys are never present in the JSON.
        "het_snp_sensitivity": None,
        "het_snp_q": None,
    }
    summary.update({name: None for name in TARGET_NULL_FIELDS})
    if exome:
        hs = parse_hsmetrics(args.hsmetrics)
        if not hs:
            raise SystemExit(
                f"ERROR: {args.sample}: exome mode needs CollectHsMetrics output; "
                f"--hsmetrics {args.hsmetrics!r} is missing or empty."
            )
        summary.update(hs)
    dedup = parse_dupblaster(args.dupblaster_stats)
    if not dedup:
        dedup = parse_flagstat(args.flagstat)
    summary.update(dedup)
    summary.update(parse_read_counts(args.read_counts, align_total_reads))
    summary["read_pair_optical_duplicates"] = int(
        metrics.get("optical", {}).get("read_pair_optical_duplicates", 0)
    )

    # Theoretical het-SNP sensitivity. Both histograms were added to bskryb-qc alongside
    # this code, so a JSON from an older image simply leaves the two fields null rather
    # than failing the run -- a resume that mixes image generations must not die here.
    sensitivity = het_sensitivity.het_snp_sensitivity(
        coverage.get("DEPTH_HISTOGRAM"), quality.get("QUALITY_HISTOGRAM")
    )
    if sensitivity is None:
        print(
            f"{args.sample}: het_snp_sensitivity not computed "
            "(DEPTH_HISTOGRAM/QUALITY_HISTOGRAM absent; pre-0.1.2 bskryb-qc image)"
        )
    summary["het_snp_sensitivity"] = sensitivity
    summary["het_snp_q"] = het_sensitivity.het_snp_q(sensitivity)

    spec = qc_scoring.spec_wes() if exome else qc_scoring.spec_wgs()
    read_gate = spec["4"]["total_reads"]["ge"]
    below_floor = qc_scoring.depth_floor(
        summary.get("total_reads"),
        min_value=DEPTH_FLOOR_FRAC_OF_TARGET * read_gate,
    )
    verdict = qc_scoring.assign_tier(summary, spec, below_floor=below_floor)
    summary["qc_score"] = verdict["qc_score"]
    summary["qc_status"] = verdict["qc_status"]

    if int(qc_scoring.qc_numeric(summary.get("total_reads")) or 0) == 0:
        summary["qc_status"] = "FAIL"
        summary["qc_score"] = None
        verdict["qc_label"] = None
        verdict["blocking_metric"] = "no_reads"

    print(
        f"{args.sample}: qc_score={summary['qc_score']} "
        f"qc_status={summary['qc_status']} "
        f"blocking={verdict['blocking_metric']!r}"
    )

    canonical_row = {field.name: summary.get(field.name) for field in WGSQC_SCHEMA}
    table = pa.Table.from_pylist([canonical_row], schema=WGSQC_SCHEMA)

    out_dir = (
        f"{prefix}qc_summary/workspace={args.workspace}/workflow_id={args.workflow_id}"
        f"/biosample={args.sample}"
    )
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(table, os.path.join(out_dir, "output.parquet"))

    table.to_pandas().to_csv(
        f"{args.sample}_{prefix}qc_metrics.tsv", sep="\t", index=False, na_rep="NA"
    )
    pd.DataFrame([
        {
            "biosampleName": args.sample,
            "qc_status": summary["qc_status"],
            "pipeline": summary["pipeline"],
            "pipeline_version": summary["pipeline_version"],
        }
    ]).to_csv(f"{args.sample}_status.csv", index=False)

    print(f"wrote {out_dir}/output.parquet")


if __name__ == "__main__":
    main()
