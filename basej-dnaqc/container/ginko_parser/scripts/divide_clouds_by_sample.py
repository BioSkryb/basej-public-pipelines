import pandas as pd
import sys
from argparse import ArgumentParser

def parse_args(args_list):
    parser = ArgumentParser()
    parser.add_argument("-c","--clouds", type=str, required=True, help="CNVs per bin")
    parser.add_argument("-v1","--cnv1", type=str, required=True, help="CNVs per bin")
    parser.add_argument("-v2","--cnv2", type=str, required=True, help="CNVs per bin")
    parser.add_argument("-b","--bins", type=str, required=True, help="positions of bins")
    return parser.parse_args(args_list)


def parse_clouds(filename):
    """
    parse clouds to dataframe
    input file: CNVs per bin with multiple samples
    input columns:
        x - number of bin
        y - CNV of bin
        SampleId - ID of biosample
        chr - number of chromosome
        
    output: multiple dfs grouped by SampleId
    """
    if(filename):
        clouds_df = pd.read_csv(filename, sep="\t")
        dfs = [x for _, x in clouds_df.groupby(clouds_df['SampleId'])]
        # for df in dfs:
            # print(df.head(5))
        return dfs
    else:
        return None
    
def parse_bins(filename):
    """
    parse bins to df
    input files: positions of bins
    input columns:
        CHR - number of chromosome => will need to rename it
        END - position on chromosome where bin ends
    """
    if filename:
        bins_df = pd.read_csv(filename, sep="\t")
        # rename columns for merging
        bins_df.rename(columns={"CHR" : "chr", "END" : "end"}, inplace=True)
        # add bin number in order to merge properly
        bins_df['x'] = range(1, len(bins_df) + 1)
        # print(bins_df.head(5))
        return bins_df
    else:
        return None
    
def parse_cnvs(cnv1_file, cnv2_file):
    if(cnv1_file and cnv2_file):
        cnv1_df = pd.read_csv(cnv1_file, sep="\t", names=["chr", "start", "end", "SampleId", "cnv"])
        cnv1_dfs = [x for _, x in cnv1_df.groupby(cnv1_df['SampleId'])]
        cnv2_df = pd.read_csv(cnv2_file, sep="\t", names=["chr", "start", "end", "SampleId", "cnv"])
        cnv2_dfs = [x for _, x in cnv2_df.groupby(cnv2_df['SampleId'])]
        # for df in cnv2_dfs:
        #     print(df.head(5))
        return cnv1_dfs, cnv2_dfs
    else:
        return None

# def save_tsvs(dfs):
#     """
#     save multiple dfs into tsvs
#     input: dfs
#     output: none

#     """
#     for df in dfs:
#         sample_id = str(df[['SampleId']]['SampleId'].unique()[0])
#         # print(sample_id)
#         output_file_name = sample_id + "_clouds.tsv"
#         df.to_csv(output_file_name, sep="\t", index=False)
        
        
def join_clouds_and_bins(clouds, bins_df):
    for df in clouds:
        if (df is None) or (bins_df is None):
            raise Exception("Could not find any dfs to join")
        sample_id = str(df[['SampleId']]['SampleId'].unique()[0])
        output_file_name = sample_id + "_dots.tsv"
        df.drop(columns=["SampleId"], inplace=True)
        df['chr'] = df['chr'].astype(str)
        df["chr"] = "chr" + df["chr"]
        joined_df = pd.merge(df, bins_df, on=["x", "chr"], how="inner")
        joined_df.rename(columns={"x" : "bin_number", "y" : "cnv"}, inplace=True)
        joined_df = joined_df[["chr", "end", "cnv", "bin_number"]]
        joined_df.to_csv(output_file_name, sep="\t", index=False)
    
        
def save_cnvs(cnv1_dfs, cnv2_dfs):
    for df in cnv1_dfs:
        sample_id = str(df[['SampleId']]['SampleId'].unique()[0])
        # print(sample_id)
        output_file_name = sample_id + "_CNV1.tsv"
        df.to_csv(output_file_name, sep="\t", index=False)
        
    for df in cnv2_dfs:
        sample_id = str(df[['SampleId']]['SampleId'].unique()[0])
        # print(sample_id)
        output_file_name = sample_id + "_CNV2.tsv"
        df.to_csv(output_file_name, sep="\t", index=False)


def main():
    args = parse_args(sys.argv[1:])
    clouds = parse_clouds(args.clouds)
    bins_df = parse_bins(args.bins)
    cnv1_dfs, cnv2_dfs = parse_cnvs(args.cnv1, args.cnv2)
    join_clouds_and_bins(clouds, bins_df)
    # save_tsvs(clouds)
    save_cnvs(cnv1_dfs, cnv2_dfs)

if __name__ == '__main__':
    sys.exit(main())