#!/usr/bin/env python3
"""Stage per-sample mutsig SNP tables under mat_files/ as .txt and run SigProfilerMatrixGenerator."""
import glob
import os
import shutil
import sys

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen


def main():
    if len(sys.argv) < 6:
        print(
            "Usage: run_matrix_generator.py <project> <reference_genome> <exome_true_false> <tsv_glob> <output_dir_name>",
            file=sys.stderr,
        )
        sys.exit(1)

    project = sys.argv[1]
    reference_genome = sys.argv[2]
    exome = sys.argv[3].lower() in ("1", "true", "yes")
    tsv_glob = sys.argv[4]
    out_dir_name = sys.argv[5]

    os.makedirs("mat_files", exist_ok=True)

    paths = sorted(glob.glob(tsv_glob))
    if not paths:
        print(f"No TSV files matched: {tsv_glob}", file=sys.stderr)
        sys.exit(1)

    for p in paths:
        base = os.path.basename(p)
        if base.endswith(".tsv"):
            dst = os.path.join("mat_files", base[:-4] + ".txt")
        else:
            dst = os.path.join("mat_files", base + ".txt")
        shutil.copy2(p, dst)

    mat_dir = os.path.abspath("mat_files") + "/"
    out_dir = os.path.abspath(out_dir_name) + "/"

    matGen.SigProfilerMatrixGeneratorFunc(
        project,
        reference_genome,
        mat_dir,
        exome=exome,
        bed_file=None,
        chrom_based=False,
        plot=False,
        tsb_stat=False,
        seqInfo=False,
        output_directory=out_dir,
    )


if __name__ == "__main__":
    main()
