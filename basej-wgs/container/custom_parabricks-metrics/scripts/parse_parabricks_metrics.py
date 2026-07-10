import argparse
from multiprocessing import Pool, cpu_count
import os
import glob
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def parse_alignment_metrics(folder_path):
    """
    Parse the alignment.txt file in a given folder.
    """
    metrics = {}
    file_path = os.path.join(folder_path, "alignment.txt")

    if os.path.exists(file_path):
        try:
            with open(file_path, "r") as f:
                lines = f.readlines()

            # Find the header (first non-comment line)
            header_line = next(
                line.strip() for line in lines if line.strip() and not line.startswith("#")
            )

            # Find the data line (line starting with 'PAIR' or 'UNPAIRED')
            data_line = next(
                line.strip()
                for line in lines
                if line.startswith("PAIR") or line.startswith("UNPAIRED")
            )

            # Convert the header and data into a dictionary
            headers = header_line.split("\t")
            data = data_line.split("\t")
            metrics = dict(zip(headers, data))

        except Exception as e:
            print(f"Warning: Could not parse alignment metrics from {file_path}: {e}")

    return metrics


def parse_dedup_metrics(folder_path):
    """
    Parse the dedup_metrics.txt file in a given folder.
    """
    metrics = {}
    file_path = os.path.join(folder_path, "dedup_metrics.txt")

    if os.path.exists(file_path):
        try:
            with open(file_path, "r") as f:
                lines = f.readlines()

            header_line = next(
                line.strip()
                for line in lines
                if line.startswith("LIBRARY")
            )
            data_idx = lines.index(header_line + "\n") + 1
            data_line = lines[data_idx].strip()

            headers = header_line.split("\t")
            data = data_line.split("\t")
            metrics = dict(zip(headers, data))

        except Exception as e:
            print(f"Warning: Could not parse dedup metrics from {file_path}: {e}")

    return metrics


def parse_insert_size_metrics(folder_path):
    """
    Parse the insert_size_metrics.txt file in a given folder.
    """
    metrics = {}
    file_path = os.path.join(folder_path, "insert_size_metrics.txt")

    if os.path.exists(file_path):
        try:
            with open(file_path, "r") as f:
                lines = f.readlines()

            header_line = next(
                line.strip()
                for line in lines
                if line.startswith("MEDIAN_INSERT_SIZE")
            )
            data_idx = lines.index(header_line + "\n") + 1
            data_line = lines[data_idx].strip()

            headers = header_line.split("\t")
            data = data_line.split("\t")
            metrics = dict(zip(headers, data))

        except Exception as e:
            print(f"Warning: Could not parse insert size metrics from {file_path}: {e}")

    return metrics


def parse_sample_folder(folder_path):
    """
    Parse all metrics files in a sample folder and return combined metrics.
    """
    sample_name = os.path.basename(folder_path)
    metrics = {"sample_name": sample_name}

    alignment = parse_alignment_metrics(folder_path)
    dedup = parse_dedup_metrics(folder_path)
    insert_size = parse_insert_size_metrics(folder_path)

    metrics.update(alignment)
    metrics.update(dedup)
    metrics.update(insert_size)

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Parse Parabricks/Sentieon metrics")
    parser.add_argument(
        "-i", "--input_dir", required=True, help="Input directory with sample folders"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output parquet file path"
    )
    parser.add_argument(
        "--threads", type=int, default=cpu_count(), help="Number of threads"
    )
    args = parser.parse_args()

    # Find all sample folders
    sample_folders = sorted(glob.glob(os.path.join(args.input_dir, "*")))
    sample_folders = [f for f in sample_folders if os.path.isdir(f)]

    if not sample_folders:
        print(f"No sample folders found in {args.input_dir}")
        return

    # Parse metrics in parallel
    with Pool(args.threads) as pool:
        results = pool.map(parse_sample_folder, sample_folders)

    # Convert to DataFrame and write Parquet
    df = pd.DataFrame(results)
    table = pa.Table.from_pandas(df)
    pq.write_table(table, args.output)
    print(f"Wrote {len(df)} samples to {args.output}")


if __name__ == "__main__":
    main()
