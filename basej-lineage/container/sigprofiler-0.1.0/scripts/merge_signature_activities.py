#!/usr/bin/env python3
"""
merge_signature_activities.py

Merges per-sample Assignment_Solution_Activities.txt files produced by
SigProfilerAssignment into a single cohort-level SamplesXSignatures matrix.

Each input directory is expected to contain:
    Assignment_Solution/Activities/Assignment_Solution_Activities.txt

The output is a tab-separated file with samples as rows and COSMIC
signatures as columns (union of all signatures found; missing values
filled with 0).
"""

import argparse
import glob
import os
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-sample SigProfilerAssignment activity files into a cohort matrix."
    )
    parser.add_argument(
        "--input_dirs",
        nargs="+",
        required=True,
        help=(
            "One or more SigProfilerAssignment output directories, each containing "
            "Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for the merged SamplesXSignatures matrix (TSV)",
    )
    args = parser.parse_args()

    dfs = []
    for d in args.input_dirs:
        activities_file = os.path.join(
            d,
            "Assignment_Solution",
            "Activities",
            "Assignment_Solution_Activities.txt",
        )
        if not os.path.isfile(activities_file):
            # Fallback: recursive search in case the directory structure differs
            matches = glob.glob(
                os.path.join(d, "**", "Assignment_Solution_Activities.txt"),
                recursive=True,
            )
            if not matches:
                print(
                    f"WARNING: No Assignment_Solution_Activities.txt found in '{d}', skipping.",
                    file=sys.stderr,
                )
                continue
            activities_file = matches[0]

        df = pd.read_csv(activities_file, sep="\t", index_col=0)
        dfs.append(df)

    if not dfs:
        raise RuntimeError(
            "No Assignment_Solution_Activities.txt files were found in any input directory."
        )

    # Concatenate rows; take the union of all signature columns.
    # Samples absent for a given signature get 0.
    merged = pd.concat(dfs, axis=0, sort=False).fillna(0).astype(int)
    merged.index.name = "Samples"
    merged.to_csv(args.output, sep="\t")

    print(
        f"Merged {len(dfs)} sample(s) → {merged.shape[0]} rows x {merged.shape[1]} signatures.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
