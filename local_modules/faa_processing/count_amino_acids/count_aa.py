#!/usr/bin/env python3

from __future__ import annotations
import argparse
import csv
from collections import Counter
from pathlib import Path

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")  # standard 20
EXTRA = []
# EXTRA = ["X", "B", "Z", "J", "U", "O", "*", "-"]  # optional/rare/stop/gap
ALL_AA = AA_ORDER + EXTRA

def iter_fasta_seqs(path: Path):
    """Yield sequence strings from a FASTA file (no headers)."""
    seq_chunks = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_chunks:
                    yield "".join(seq_chunks)
                    seq_chunks = []
            else:
                seq_chunks.append(line)
        if seq_chunks:
            yield "".join(seq_chunks)

def count_amino_acids_in_faa(path: Path) -> Counter[str]:
    counts = Counter()
    for seq in iter_fasta_seqs(path):
        counts.update(seq.strip().upper())
    return counts

def genome_name_from_path(path: Path) -> str:
    name = path.name
    for suf in [".faa.gz", ".gz", ".faa", ".fasta", ".fa"]:
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name

def main():
    ap = argparse.ArgumentParser(
        description="Count amino acids in .faa files and write a CSV (rows=genomes, cols=amino acids)."
    )
    ap.add_argument(
        "inputs",
        nargs="+",
        help="Input .faa files and/or directories containing .faa files",
    )
    ap.add_argument(
        "-o",
        "--out",
        default="amino_acid_counts.csv",
        help="Output CSV path (default: amino_acid_counts.csv)",
    )
    ap.add_argument(
        "--include-nonstandard",
        action="store_true",
        help="Include non-standard symbols found (not just 20 AAs + common extras)",
    )
    ap.add_argument(
        "--mode",
        choices=["counts", "percent"],
        default="counts",
        help="Output raw counts or percent of total amino acids per genome (default: counts).",
    )
    ap.add_argument(
        "--decimals",
        type=int,
        default=3,
        help="Decimal places for percent mode (default: 3).",
    )
    args = ap.parse_args()

    # collect .faa files
    faa_files: list[Path] = []
    for item in args.inputs:
        p = Path(item)
        if p.is_dir():
            faa_files.extend(sorted(p.glob("*.faa")))
        elif p.is_file():
            faa_files.append(p)
        else:
            raise FileNotFoundError(f"Not found: {p}")

    if not faa_files:
        raise SystemExit("No .faa files found.")

    # count per genome
    per_genome: dict[str, Counter[str]] = {}
    all_symbols = set()

    for fp in faa_files:
        genome = genome_name_from_path(fp)
        c = count_amino_acids_in_faa(fp)
        per_genome[genome] = c
        all_symbols |= set(c.keys())

    # choose column order (amino acids)
    if args.include_nonstandard:
        others = sorted(s for s in all_symbols if s not in ALL_AA and s not in AA_ORDER)
        aa_columns = AA_ORDER + [s for s in ALL_AA if s not in AA_ORDER] + others
    else:
        aa_columns = ALL_AA

    genomes_sorted = sorted(per_genome.keys())

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)

        # header
        w.writerow(["genome"] + aa_columns)

        # rows = genomes
        for genome in genomes_sorted:
            counts = per_genome[genome]

            if args.mode == "counts":
                row_vals = [counts.get(aa, 0) for aa in aa_columns]
            else:
                total = sum(counts.get(aa, 0) for aa in aa_columns)
                if total == 0:
                    row_vals = [0.0 for _ in aa_columns]
                else:
                    row_vals = [
                        round((counts.get(aa, 0) / total) * 100.0, args.decimals)
                        for aa in aa_columns
                    ]

            w.writerow([genome] + row_vals)

    print(f"Wrote: {out_path}  ({len(genomes_sorted)} genomes, mode={args.mode})")

if __name__ == "__main__":
    main()