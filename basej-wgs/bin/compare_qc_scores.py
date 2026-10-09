#!/usr/bin/env python3
"""
Concordance check between the authoritative Python QC tier (qc_scoring.py, written to the
summary TSV as qc_score) and the R consensus score (dna_qc_plot.R / wgs_qc_plot.R /
wes_qc_plot.R, *-QC_ConsensusScores.txt CompositeScore).

The Python verdict is what lands in Parquet / per_biosample_status; R only drives the plots.
Both implement the same specs, so they should agree. This writes a per-sample table and
reports mismatches. By default a mismatch is a warning (the plots still publish); pass
--strict to fail the task instead (used by nf-tests).

Zero-read samples are skipped: their Python score is intentionally null (FAIL, no score).
"""
import argparse
import sys

import pandas as pd


def compare(summary_path, scores_path, sample_col="biosample"):
    summary = pd.read_csv(summary_path, sep="\t")
    scores = pd.read_csv(scores_path, sep="\t")
    if "SampleId" not in scores.columns or "CompositeScore" not in scores.columns:
        raise SystemExit(f"ERROR: {scores_path} lacks SampleId/CompositeScore columns")
    py = summary[[sample_col, "qc_score", "total_reads"]].rename(columns={sample_col: "SampleId"})
    py["SampleId"] = py["SampleId"].astype(str)
    r = scores[["SampleId", "CompositeScore"]].copy()
    r["SampleId"] = r["SampleId"].astype(str).str.replace("_sorted$", "", regex=True)
    merged = py.merge(r, on="SampleId", how="left")
    merged["qc_score"] = pd.to_numeric(merged["qc_score"], errors="coerce")
    merged["CompositeScore"] = pd.to_numeric(merged["CompositeScore"], errors="coerce")
    merged["total_reads"] = pd.to_numeric(merged["total_reads"], errors="coerce").fillna(0)

    def verdict(row):
        if row["total_reads"] == 0:
            return "skipped_zero_reads"
        if pd.isna(row["CompositeScore"]):
            return "missing_in_r"
        if pd.isna(row["qc_score"]):
            return "missing_in_python"
        return "match" if int(row["qc_score"]) == int(row["CompositeScore"]) else "mismatch"

    merged["concordance"] = merged.apply(verdict, axis=1)
    return merged


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--sample-col", default="biosample")
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args(argv)

    merged = compare(args.summary, args.scores, args.sample_col)
    merged.to_csv(args.out, sep="\t", index=False, na_rep="NA")
    bad = merged[merged["concordance"].isin(["mismatch", "missing_in_python"])]
    for _, row in bad.iterrows():
        print(f"WARNING: QC score {row['concordance']} for {row['SampleId']}: "
              f"python={row['qc_score']} R={row['CompositeScore']}", file=sys.stderr)
    print(f"QC score concordance: {(merged['concordance'] == 'match').sum()}/{len(merged)} match")
    if args.strict and not bad.empty:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
