import pandas as pd
import sys
from argparse import ArgumentParser


def parse_args(args_list):
    parser = ArgumentParser()
    parser.add_argument("-c", "--clouds", type=str, required=True, help="CNVs per bin")
    parser.add_argument("-v1", "--cnv1", type=str, required=True, help="CNVs per bin")
    parser.add_argument("-v2", "--cnv2", type=str, required=True, help="CNVs per bin")
    parser.add_argument("-b", "--bins", type=str, required=True, help="positions of bins")
    return parser.parse_args(args_list)


def parse_clouds(filename):
    """
    parse clouds to dataframe
    input file: CNVs per bin with multiple samples
    input columns:
        x - number of bin
        y - CNV of bin
        SampleId - sample identifier
    """
    df = pd.read_csv(filename, sep="\t")
    return df


def parse_cnv(filename):
    """
    parse CNV file to dataframe
    """
    df = pd.read_csv(filename, sep="\t")
    return df


def parse_bins(filename):
    """
    parse bin positions file
    """
    df = pd.read_csv(filename, sep="\t", header=None)
    return df


def main():
    args = parse_args(sys.argv[1:])

    clouds = parse_clouds(args.clouds)
    cnv1 = parse_cnv(args.cnv1)
    cnv2 = parse_cnv(args.cnv2)
    bins = parse_bins(args.bins)

    # Get unique sample IDs
    samples = clouds["SampleId"].unique()

    for sample in samples:
        sample_clouds = clouds[clouds["SampleId"] == sample]
        sample_cnv1 = cnv1[[c for c in cnv1.columns if sample in c or c in ["CHR", "START", "END"]]] if len(cnv1.columns) > 0 else pd.DataFrame()
        sample_cnv2 = cnv2[[c for c in cnv2.columns if sample in c or c in ["CHR", "START", "END"]]] if len(cnv2.columns) > 0 else pd.DataFrame()

        # Write per-sample outputs
        sample_clouds.to_csv(f"{sample}_clouds.tsv", sep="\t", index=False)
        if not sample_cnv1.empty:
            sample_cnv1.to_csv(f"{sample}_CNV1.tsv", sep="\t", index=False)
        if not sample_cnv2.empty:
            sample_cnv2.to_csv(f"{sample}_CNV2.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
