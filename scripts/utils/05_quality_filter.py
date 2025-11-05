#!/usr/bin/env python3
"""
Filter phage sequences based on CheckV quality metrics and extract high-quality contigs.

Usage:
    python 05_quality_filter.py \
        --checkv_summary /path/to/sample/quality_summary.tsv \
        --filtered_summary /path/to/sample/filtered_quality_summary.tsv \
        --viruses_fna /path/to/sample/viruses.fna \
        --proviruses_fna /path/to/sample/proviruses.fna \
        --out_dir /path/to/sample/
"""

import pandas as pd
from Bio import SeqIO
import argparse
import os
import re


# ------------------------------------------------------------
# Filtering logic
# ------------------------------------------------------------
def filter_phage_sequences(input_file, output_file):
    """
    Filter phage sequences based on length, completeness, and viral gene criteria.

    Criteria:
    1. Sequences must be:
       - ≥ 5 kbp long (proviral_length if provirus is 'Yes', else contig_length)
       OR
       - Predicted as ≥ 50% complete
       - If < 5 kbp, must be at least 1 kbp long
    2. Must have at least one viral gene
    """
    df = pd.read_csv(input_file, sep="\t")

    filtered_df = df[
        (
            (
                ((df["provirus"] == "Yes") & (df["proviral_length"].fillna(0) >= 5000))
                | ((df["provirus"] != "Yes") & (df["contig_length"].fillna(0) >= 5000))
            )
            | (
                (df["completeness"].fillna(0) >= 50)
                & (
                    ((df["provirus"] == "Yes") & (df["proviral_length"].fillna(0) >= 1000))
                    | ((df["provirus"] != "Yes") & (df["contig_length"].fillna(0) >= 1000))
                )
            )
        )
        & (df["viral_genes"].fillna(0) > 0)
    ]

    filtered_df.to_csv(output_file, sep="\t", index=False)

    print(f"[INFO] Total input sequences: {len(df)}")
    print(f"[INFO] Filtered sequences: {len(filtered_df)}")

    return filtered_df


# ------------------------------------------------------------
# Sequence extraction
# ------------------------------------------------------------
def extract_and_rename_contigs(filtered_df, viruses_file, proviruses_file, output_viruses_file, output_proviruses_file):
    """
    Extract and rename contigs from viruses and proviruses files based on filtered dataframe.
    """
    filtered_ids = set(filtered_df["contig_id"])

    # --- Proviruses ---
    with open(proviruses_file, "r") as in_handle, open(output_proviruses_file, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            match = re.match(r"(.*?)_\d+\s(\d+)-(\d+)/\d+", record.description)
            if match:
                base_id, start, end = match.groups()
                new_base_id = base_id.rsplit("_", 2)[0]
                if base_id in filtered_ids:
                    SeqIO.write(record, out_handle, "fasta")

    # --- Viruses ---
    with open(viruses_file, "r") as in_handle, open(output_viruses_file, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            if record.id in filtered_ids:
                SeqIO.write(record, out_handle, "fasta")


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Filter CheckV phage contigs based on completeness and length thresholds."
    )
    parser.add_argument("--checkv_summary", required=True, help="Path to CheckV quality_summary.tsv")
    parser.add_argument("--filtered_summary", required=True, help="Path to save filtered summary TSV")
    parser.add_argument("--viruses_fna", required=True, help="Path to viruses.fna file")
    parser.add_argument("--proviruses_fna", required=True, help="Path to proviruses.fna file")
    parser.add_argument("--out_dir", required=True, help="Output directory for filtered FASTAs")

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print(f"[INFO] Filtering sequences from {args.checkv_summary}")
    filtered_df = filter_phage_sequences(args.checkv_summary, args.filtered_summary)

    print(f"[INFO] Extracting filtered sequences...")
    extract_and_rename_contigs(
        filtered_df,
        args.viruses_fna,
        args.proviruses_fna,
        os.path.join(args.out_dir, "filtered_viruses.fna"),
        os.path.join(args.out_dir, "filtered_proviruses.fna"),
    )

    print(f"[✓] Completed filtering and extraction for {args.checkv_summary}")


if __name__ == "__main__":
    main()
