#!/usr/bin/env python3
"""
Extract phage contigs and prophage subsequences.

Usage:
    python 04_extract_prophages.py \
        --assembly_dir /path/to/02_qc_assemblies/filtered \
        --cons_dir /path/to/04_consolidated \
        --tech illumina --assembler megahit
"""

import csv
import argparse
from pathlib import Path

# ----------------------------
# Helper functions
# ----------------------------
def load_fasta_into_memory(fasta_path):
    """Load FASTA into dict[name] = sequence."""
    print(f"  Loading FASTA: {fasta_path.name}")
    seqs = {}
    name, chunks = None, []
    try:
        with open(fasta_path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(chunks).upper()
                    name = line[1:].strip().split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip())
            if name is not None:
                seqs[name] = "".join(chunks).upper()
    except FileNotFoundError:
        print(f"  [WARN] FASTA not found: {fasta_path}")
        return {}
    print(f"  Loaded {len(seqs)} sequences")
    return seqs


def write_fasta_entry(ofh, header, seq, width=80):
    ofh.write(f">{header}\n")
    for i in range(0, len(seq), width):
        ofh.write(seq[i:i+width] + "\n")


def parse_prophage_coordinates(contig_id):
    """Extract base contig name and coordinates if '_prophage_' is present."""
    if "_prophage_" not in contig_id:
        return contig_id, None, None

    parts = contig_id.split("_prophage_")
    if len(parts) != 2:
        return contig_id, None, None

    base_contig = parts[0]
    coords = parts[1].split("_")
    if len(coords) >= 2:
        try:
            start = int(coords[0])
            end = int(coords[1])
            return base_contig, start, end
        except ValueError:
            pass

    return contig_id, None, None


# ----------------------------
# Main
# ----------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly_dir", required=True, help="Path to assembly directory (02_qc_assemblies/filtered)")
    parser.add_argument("--cons_dir", required=True, help="Path to 04_consolidated folder")
    parser.add_argument("--tech", required=True, help="Technology: illumina, pacbio, or ont")
    parser.add_argument("--assembler", required=True, help="Assembler: megahit, spades, flye, hifiasm, autocycler")
    args = parser.parse_args()

    combo = f"{args.tech}.{args.assembler}"
    
    # Find assembly file
    assembly_dir = Path(args.assembly_dir)
    fasta_files = list(assembly_dir.glob(f"{combo}*.fa"))
    
    if not fasta_files:
        print(f"[WARN] No FASTA file found for {combo} in {assembly_dir}")
        return

    fasta_path = fasta_files[0]
    seqs = load_fasta_into_memory(fasta_path)

    if not seqs:
        print(f"[WARN] No sequences loaded from {fasta_path}")
        return

    # Find score matrix
    combo_dir = Path(args.cons_dir) / combo
    scores_path = combo_dir / f"{combo}_merged_tool_scores.tsv"
    
    if not scores_path.exists():
        print(f"[WARN] Score matrix not found: {scores_path}")
        return

    # Output
    out_dir = Path(args.cons_dir) / "phage_contigs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{combo}_contigs.fasta"

    n_ok = n_missing = n_bad = n_prophage = n_full_contig = 0

    print(f"[INFO] Extracting phage contigs for {combo}...")
    with open(scores_path, newline="") as tfh, open(out_path, "w") as ofh:
        rdr = csv.DictReader(tfh, delimiter="\t")
        if not {"Contig_ID", "Consensus_Label"}.issubset(rdr.fieldnames):
            print(f"  [WARN] Missing required columns in {scores_path}")
            return

        for row in rdr:
            contig_id = row["Contig_ID"]
            consensus_label = row["Consensus_Label"]

            # Only export phage predictions
            if consensus_label != "phage":
                continue

            base_contig, start, end = parse_prophage_coordinates(contig_id)

            if start is not None and end is not None:
                seq = seqs.get(base_contig)
                if seq is None:
                    print(f"  [WARN] Missing base contig: {base_contig}")
                    n_missing += 1
                    continue

                start0, end0 = start - 1, end
                if start0 < 0 or end0 <= start0 or end0 > len(seq):
                    print(f"  [WARN] Invalid coordinates for '{contig_id}' ({start}-{end})")
                    n_bad += 1
                    continue

                subseq = seq[start0:end0]
                write_fasta_entry(ofh, contig_id, subseq)
                n_prophage += 1
                n_ok += 1

            else:
                seq = seqs.get(contig_id)
                if seq is None:
                    print(f"  [WARN] Missing contig: {contig_id}")
                    n_missing += 1
                    continue

                write_fasta_entry(ofh, contig_id, seq)
                n_full_contig += 1
                n_ok += 1

    print(f"  ✅ Output: {out_path.name}")
    print(f"    - Total sequences:    {n_ok}")
    print(f"    - Prophages:          {n_prophage}")
    print(f"    - Full contigs:       {n_full_contig}")
    print(f"    - Missing:            {n_missing}")
    print(f"    - Bad coordinates:    {n_bad}")
    print(f"[✓] Completed {combo}")


if __name__ == "__main__":
    main()