import os
import glob
import argparse
from collections import Counter
import pandas as pd


def parse_kofamscan_file(file_path):
    file_name = os.path.basename(file_path)
    genome = file_name.replace("_kofamscan.txt", "").replace(".tsv", "").replace(".txt", "")
    ko_counts = Counter()

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if line.startswith("*"):
                parts = line[1:].split(maxsplit=5)
                if len(parts) >= 2:
                    ko = parts[1]
                    ko_counts[ko] += 1

    return {"genome": genome, **ko_counts}


def collect_input_files(input_dir=None, input_files=None):
    if input_dir is not None:
        files = sorted(glob.glob(os.path.join(input_dir, "*.tsv")))
    else:
        files = sorted(input_files)

    if not files:
        raise ValueError("No input TSV files were found.")

    return files


def build_ko_matrix(input_dir=None, input_files=None):
    files = collect_input_files(input_dir=input_dir, input_files=input_files)
    rows = [parse_kofamscan_file(file_path) for file_path in files]

    df = pd.DataFrame(rows).fillna(0)

    if not df.empty:
        ko_cols = [c for c in df.columns if c != "genome"]
        df[ko_cols] = df[ko_cols].astype(int)

        nonzero_cols = ["genome"] + [c for c in ko_cols if (df[c] != 0).any()]
        df = df[nonzero_cols]

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Build a KO count matrix from KofamScan TSV files."
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-i", "--input-dir",
        help="Directory containing KofamScan TSV files"
    )
    input_group.add_argument(
        "-f", "--input-files",
        nargs="+",
        help="One or more KofamScan TSV files"
    )

    parser.add_argument(
        "-o", "--output-file",
        required=True,
        help="Output TSV file path"
    )

    args = parser.parse_args()

    df = build_ko_matrix(
        input_dir=args.input_dir,
        input_files=args.input_files
    )
    df.to_csv(args.output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()