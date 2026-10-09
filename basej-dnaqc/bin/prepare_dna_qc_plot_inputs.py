#!/usr/bin/env python3
"""
Adapt basej-dnaqc outputs to the inputs of the shared DNA QC R plotting scripts
(custom_r_qcplots: dna_qc_plot.R, function_cnv_quadrants_qc.R).

The R scripts were written against basej-dnaqc's combined_selected_metrics table, so this
only renames dnaqc_summary columns to the names those scripts read and writes a metadata
CSV keyed by biosampleName. No metric is recomputed here.

  --summary   dnaqc_all_metrics.tsv (METRICS_TO_PARQUET)
  --input-csv the pipeline input CSV (biosampleName|biosample [,group|groups])
  --out-metrics  metrics table for the R scripts (tab-separated)
  --out-metadata metadata table for dna_qc_plot.R (comma-separated)

Column mapping (dnaqc_summary -> R):
  biosample          -> sample_name
  preseq_count       -> preseq_count
  align_pct_chimeras -> PCT_CHIMERAS
  pct_chrm           -> chrM
  ginkgo_cnv_mapd    -> MAPD_CNV_Log2
  ginkgo_cnv_skew    -> SKEW_CNV
  total_reads        -> total_reads (depth floor) and TotalReads (legacy name)
  final_reads        -> FinalReads
  insert_size        -> insert_size
  fastp_q30_rate     -> fastp_q30_rate

A missing SKEW_CNV is filled with the legacy 999 sentinel (as basej-dnaqc's QC_PLOTS does),
so the R tier for such a cell fails the skew gate. qc_scoring.py treats the same cell as
Inconclusive; compare_qc_scores.py reports any resulting disagreement. The Python verdict
is the one written to Parquet / per_biosample_status.
"""
import argparse
import sys

import pandas as pd

RENAME = {
    "biosample": "sample_name",
    "preseq_count": "preseq_count",
    "align_pct_chimeras": "PCT_CHIMERAS",
    "pct_chrm": "chrM",
    "ginkgo_cnv_mapd": "MAPD_CNV_Log2",
    "ginkgo_cnv_skew": "SKEW_CNV",
    "total_reads": "total_reads",
    "final_reads": "FinalReads",
    "insert_size": "insert_size",
    "fastp_q30_rate": "fastp_q30_rate",
}
SKEW_SENTINEL = 999
NUMERIC = ["preseq_count", "PCT_CHIMERAS", "chrM", "MAPD_CNV_Log2", "SKEW_CNV",
           "total_reads", "FinalReads", "insert_size", "fastp_q30_rate"]


def build_metrics(summary_path):
    df = pd.read_csv(summary_path, sep="\t")
    missing = [c for c in ("biosample", "preseq_count", "align_pct_chimeras", "pct_chrm",
                           "ginkgo_cnv_mapd", "ginkgo_cnv_skew", "total_reads")
               if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: {summary_path} is missing column(s): {', '.join(missing)}")
    out = pd.DataFrame({dst: df[src] if src in df.columns else pd.NA
                        for src, dst in RENAME.items()})
    for col in NUMERIC:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    # cnvSummarizer.R leaves SKEW_CNV NA when a cell has too few segments to measure
    # unevenness (e.g. chr22-only data). basej-dnaqc's QC_PLOTS fills it with the 999
    # sentinel before calling the R scripts (qc_require_metric aborts on an all-NA
    # column); keep that behaviour so the plots match the legacy pipeline.
    out["SKEW_CNV"] = out["SKEW_CNV"].fillna(SKEW_SENTINEL)
    out["TotalReads"] = out["total_reads"]
    return out


def build_metadata(input_csv, samples):
    """biosampleName [+ group] for every sample that reached the summary."""
    meta = pd.read_csv(input_csv, dtype=str).fillna("")
    if "biosampleName" not in meta.columns and "biosample" in meta.columns:
        meta = meta.rename(columns={"biosample": "biosampleName"})
    if "biosampleName" not in meta.columns:
        raise SystemExit(f"ERROR: {input_csv} has no biosampleName/biosample column")
    group_col = "group" if "group" in meta.columns else ("groups" if "groups" in meta.columns else None)
    meta["biosampleName"] = meta["biosampleName"].str.strip()
    meta = meta[meta["biosampleName"].isin(samples)].drop_duplicates("biosampleName")
    out = pd.DataFrame({"biosampleName": meta["biosampleName"]})
    # dna_qc_plot.R drops rows with a blank group, so default blanks to Group1.
    if group_col:
        out["group"] = meta[group_col].where(meta[group_col].str.strip() != "", "Group1")
    else:
        out["group"] = "Group1"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--out-metrics", required=True)
    ap.add_argument("--out-metadata", required=True)
    args = ap.parse_args(argv)

    metrics = build_metrics(args.summary)
    metadata = build_metadata(args.input_csv, set(metrics["sample_name"].astype(str)))
    if metadata.empty:
        raise SystemExit("ERROR: no input_csv biosample matches a dnaqc_summary sample")

    metrics.to_csv(args.out_metrics, sep="\t", index=False, na_rep="NA")
    metadata.to_csv(args.out_metadata, index=False)
    print(f"Wrote {len(metrics)} metric rows and {len(metadata)} metadata rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
