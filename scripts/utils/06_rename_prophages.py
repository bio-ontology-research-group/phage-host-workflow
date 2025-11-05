#!/usr/bin/env python3
"""
Rename provirus sequences based on CheckV coordinate logic and concatenate
filtered virus + provirus FASTAs into a single output.

Usage:
    python 06_rename.py \
        --input_dir /path/to/05_checkv/sample_dir \
        --output_dir /path/to/05_checkv/sample_dir
"""

import os
from pathlib import Path
from Bio import SeqIO
import argparse


# ------------------------------------------------------------
# Functions
# ------------------------------------------------------------
def rename_proviruses(input_fasta, output_fasta):
    """
    Rename provirus sequences using coordinate logic from CheckV headers.
    """
    if not input_fasta.exists():
        print(f"[WARN] Missing {input_fasta}")
        return False

    print(f"[INFO] Renaming {input_fasta.name}")
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(str(input_fasta), "fasta"):
            try:
                contig_info, coords_info = record.description.split()
                parts = contig_info.split("_")

                if "prophage" in parts:
                    # Case 1: already has "prophage" in name
                    contig_id = "_".join(parts[:-3])
                    start = int(parts[-3])

                    rel_range, total_len = coords_info.split("/")
                    rel_start, rel_end = map(int, rel_range.split("-"))

                    actual_start = start + (rel_start - 1)
                    actual_end = start + (rel_end - 1)
                    new_id = f"{contig_id}_{actual_start}_{actual_end}"
                else:
                    # Case 2: plain contig without "prophage"
                    contig_id = "_".join(parts[:-1])
                    abs_range, total_len = coords_info.split("/")
                    actual_start, actual_end = map(int, abs_range.split("-"))
                    new_id = f"{contig_id}_prophage_{actual_start}_{actual_end}"

                record.id = new_id
                record.description = ""
                SeqIO.write(record, out_handle, "fasta")

            except Exception as e:
                print(f"  ⚠️ Skipped {record.id}: {e}")
                continue

    print(f"  ✅ Renamed sequences written to {output_fasta}")
    return True


def concatenate_fasta(file1, file2, output):
    """
    Concatenate two FASTA files if they exist.
    """
    if not file1.exists() and not file2.exists():
        print(f"[WARN] No input files found for {output.name.replace('_checkv.fasta', '')}")
        return False

    with open(output, "w") as outfile:
        for f in [file1, file2]:
            if f.exists():
                with open(f) as infile:
                    outfile.write(infile.read())

    print(f"  🧬 Concatenated to {output}")
    return True


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Rename provirus sequences and concatenate CheckV FASTAs."
    )
    parser.add_argument("--input_dir", required=True, help="Directory containing filtered_*.fna files")
    parser.add_argument("--output_dir", required=True, help="Directory to write renamed and merged FASTAs")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Input and output paths
    provirus_in = input_dir / "filtered_proviruses.fna"
    provirus_out = input_dir / "filtered_proviruses_renamed.fna"
    virus_fa = input_dir / "filtered_viruses.fna"
    checkv_out = output_dir / f"{input_dir.name}_checkv.fasta"

    print(f"\n=== Processing {input_dir.name} ===")

    renamed = rename_proviruses(provirus_in, provirus_out)
    concatenate_fasta(provirus_out, virus_fa, checkv_out)

    print(f"\n[✓] Completed renaming and concatenation for {input_dir.name}")


if __name__ == "__main__":
    main()
