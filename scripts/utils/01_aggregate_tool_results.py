#!/usr/bin/env python3
"""
Consolidate and standardize phage/plasmid prediction results
for a specific technology and assembler combination.

Usage:
    python consolidate_predictions.py \
        --pred_dir /path/to/03_predictions \
        --out_dir /path/to/04_consolidated \
        --tech ont --assembler flye
"""

import os
import glob
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path

def process_one_hot_file(input_path, output_path):
    """Compute softmax scores and best_choice for a one-hot TSV file."""
    try:
        df = pd.read_csv(input_path, sep="\t", header=None, skiprows=1)
        df.columns = [
            "Sequence Name",
            "Eukaryote",
            "EukaryoteVirus",
            "Plasmid",
            "Prokaryote",
            "ProkaryoteVirus",
        ]
        score_cols = ["Eukaryote", "EukaryoteVirus", "Plasmid", "Prokaryote", "ProkaryoteVirus"]

        # Ensure numeric
        df[score_cols] = df[score_cols].astype(float)

        # 1) Best choice index (1..5)
        raw = df[score_cols].to_numpy()
        best_choice = np.argmax(raw, axis=1) + 1
        df["best_choice"] = best_choice

        # 2) Apply softmax
        shifted = raw - np.max(raw, axis=1, keepdims=True)
        exp_scores = np.exp(shifted)
        softmax = exp_scores / np.sum(exp_scores, axis=1, keepdims=True)
        df[score_cols] = softmax

        # 3) Compute margin (confidence)
        sorted_softmax = np.sort(softmax, axis=1)
        df["softmax_margin"] = sorted_softmax[:, -1] - sorted_softmax[:, -2]

        # Save
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df.to_csv(output_path, sep="\t", index=False)
        print(f"  ✓ Processed one-hot file: {os.path.basename(input_path)}")
        return True

    except Exception as e:
        print(f"  ⚠️ Error processing {input_path}: {e}")
        return False

def consolidate_genomad(pred_dir, out_dir, tech, assembler):
    """
    Consolidate Genomad virus and plasmid summaries.
    For viruses, standardize provirus names to contig_prophage_start_end.
    """
    base_dir = Path(pred_dir) / f"genomad/{tech}.{assembler}"
    if not base_dir.exists():
        print(f"[WARN] No Genomad results found for {tech}.{assembler}")
        return

    print(f"[INFO] Processing Genomad results for {tech}.{assembler}...")
    virus_files = glob.glob(str(base_dir / "*_summary/*virus_summary.tsv"))
    plasmid_files = glob.glob(str(base_dir / "*_summary/*plasmid_summary.tsv"))

    outpath = ensure_dir(Path(out_dir) / f"{tech}.{assembler}")

    # --- VIRUSES ---
    if virus_files:
        all_dfs = []
        for file in virus_files:
            sample = Path(file).parts[-3]
            try:
                df = pd.read_csv(file, sep="\t")
                df["sample"] = sample

                # Standardize prophage names
                for idx, row in df.iterrows():
                    seq_name = str(row.get("seq_name", ""))
                    topology = str(row.get("topology", ""))

                    if topology.lower() == "provirus":
                        start = end = None

                        # Case 1: Coordinates column exists
                        if pd.notna(row.get("coordinates")) and "-" in str(row["coordinates"]):
                            coords = str(row["coordinates"]).split("-")
                            start, end = coords[0], coords[1]

                        # Case 2: Coordinates encoded in seq_name (e.g., contig|provirus_1234_5678)
                        elif "|provirus_" in seq_name:
                            parts = seq_name.split("|provirus_")[1].split("_")
                            if len(parts) >= 2:
                                start, end = parts[0], parts[1]

                        # Rename only if valid coordinates found
                        if start and end:
                            base_contig = seq_name.split("|")[0]
                            standardized_name = f"{base_contig}_prophage_{start}_{end}"
                            df.at[idx, "seq_name"] = standardized_name

                all_dfs.append(df)
                print(f"  ✓ Processed {sample}: {len(df)} entries")

            except Exception as e:
                print(f"  ⚠️ Error processing {file}: {e}")

        if all_dfs:
            virus_out = Path(outpath) / f"genomad_{tech}_{assembler}_virus_summary.tsv"
            pd.concat(all_dfs, ignore_index=True).to_csv(virus_out, sep="\t", index=False)
            print(f"  ✓ Saved {virus_out}")

    # --- PLASMIDS ---
    if plasmid_files:
        all_dfs = []
        for file in plasmid_files:
            sample = Path(file).parts[-3]
            try:
                df = pd.read_csv(file, sep="\t")
                df["sample"] = sample
                all_dfs.append(df)
                print(f"  ✓ Processed {sample}: {len(df)} plasmid entries")
            except Exception as e:
                print(f"  ⚠️ Error processing {file}: {e}")

        if all_dfs:
            plasmid_out = Path(outpath) / f"genomad_{tech}_{assembler}_plasmid_summary.tsv"
            pd.concat(all_dfs, ignore_index=True).to_csv(plasmid_out, sep="\t", index=False)
            print(f"  ✓ Saved {plasmid_out}")

def consolidate_vibrant(pred_dir, out_dir, tech, assembler):
    base_dir = Path(pred_dir) / f"vibrant/{tech}.{assembler}"
    if not base_dir.exists():
        print(f"[WARN] No VIBRANT results for {tech}.{assembler}")
        return

    print(f"[INFO] Processing VIBRANT results for {tech}.{assembler}...")
    summary_files = glob.glob(str(base_dir / "VIBRANT*/VIBRANT_results*/VIBRANT_summary_normalized*.tsv"))
    coord_files = glob.glob(str(base_dir / "VIBRANT*/VIBRANT_results*/VIBRANT_integrated_prophage_coordinates*.tsv"))

    if not summary_files:
        print("  ⚠ No VIBRANT summary found")
        return

    coord_map = {}
    if coord_files:
        coord_df = pd.read_csv(coord_files[0], sep="\t")
        for _, row in coord_df.iterrows():
            try:
                frag, scaf = row["fragment"], row["scaffold"]
                start, end = row["nucleotide start"], row["nucleotide stop"]
                coord_map[frag] = f"{scaf}_prophage_{start}_{end}"
            except KeyError:
                continue

    dfs = []
    for file in summary_files:
        df = pd.read_csv(file, sep="\t")
        df["sample"] = f"{tech}.{assembler}"
        if coord_map:
            df["scaffold"] = df["scaffold"].apply(lambda x: coord_map.get(x, x))
        dfs.append(df)

    vibrant_out = Path(out_dir) / f"{tech}.{assembler}/vibrant_{tech}_{assembler}_summary.tsv"
    ensure_dir(vibrant_out.parent)
    pd.concat(dfs, ignore_index=True).to_csv(vibrant_out, sep="\t", index=False)
    print(f"  ✓ Saved {vibrant_out}")

def consolidate_virsorter2(pred_dir, out_dir, tech, assembler):
    base_dir = Path(pred_dir) / f"virsorter2/{tech}.{assembler}"
    score_file = base_dir / "final-viral-score.tsv"
    boundary_file = base_dir / "final-viral-boundary.tsv"
    if not score_file.exists():
        print(f"[WARN] No VirSorter2 output for {tech}.{assembler}")
        return

    print(f"[INFO] Processing VirSorter2 results for {tech}.{assembler}...")
    score_df = pd.read_csv(score_file, sep="\t")
    score_df["sample"] = f"{tech}.{assembler}"

    if boundary_file.exists():
        boundary_df = pd.read_csv(boundary_file, sep="\t")
        prophage_map = {}
        for _, row in boundary_df.iterrows():
            if "trim_bp_start" in row and "trim_bp_end" in row:
                start, end = row["trim_bp_start"], row["trim_bp_end"]
                base = row["seqname"]
                prophage_map[row["seqname_new"]] = f"{base}_prophage_{start}_{end}"
        score_df["seqname"] = score_df["seqname"].apply(lambda x: prophage_map.get(x, x))

    vir_out = Path(out_dir) / f"{tech}.{assembler}/virsorter2_{tech}_{assembler}_summary.tsv"
    ensure_dir(vir_out.parent)
    score_df.to_csv(vir_out, sep="\t", index=False)
    print(f"  ✓ Saved {vir_out}")

def consolidate_deepmc(pred_dir, out_dir, tech, assembler):
    base_dir = Path(pred_dir) / f"deepmicroclass/{tech}.{assembler}"
    
    # Look for one-hot files first
    one_hot_files = list(base_dir.glob("*one-hot*.tsv"))

    print(f"[INFO] Processing DeepMicroClass for {tech}.{assembler}...")
    
    # Process one-hot files to generate softmax
    softmax_files = []
    for one_hot_file in one_hot_files:
        # Generate softmax output path
        softmax_file = base_dir / f"{tech}.{assembler}_softmax.tsv"
        
        # Process if softmax doesn't exist
        if not softmax_file.exists():
            if process_one_hot_file(str(one_hot_file), str(softmax_file)):
                softmax_files.append(softmax_file)
        else:
            print(f"  ⊗ Softmax file exists, skipping: {softmax_file.name}")
            softmax_files.append(softmax_file)
    
    # Now consolidate the softmax files
    if not softmax_files:
        print(f"  ⚠️ No softmax files generated")
        return
    
    dfs = []
    for f in softmax_files:
        try:
            df = pd.read_csv(f, sep="\t")
            df["sample"] = f"{tech}.{assembler}"
            dfs.append(df)
        except Exception as e:
            print(f"  ⚠️ Error reading {f}: {e}")
    
    if dfs:
        out = Path(out_dir) / f"{tech}.{assembler}/deepmicroclass_{tech}_{assembler}_softmax.tsv"
        ensure_dir(out.parent)
        pd.concat(dfs, ignore_index=True).to_csv(out, sep="\t", index=False)
        print(f"  ✓ Saved consolidated: {out}")

def consolidate_phamer(pred_dir, out_dir, tech, assembler):
    base_dir = Path(pred_dir) / f"phamer/{tech}.{assembler}/final_prediction"
    file = base_dir / "phamer_prediction.tsv"
    if not file.exists():
        print(f"[WARN] No Phamer results for {tech}.{assembler}")
        return

    print(f"[INFO] Processing Phamer for {tech}.{assembler}...")
    df = pd.read_csv(file, sep="\t")
    df["sample"] = f"{tech}.{assembler}"
    out = Path(out_dir) / f"{tech}.{assembler}/phamer_{tech}_{assembler}_summary.tsv"
    ensure_dir(out.parent)
    df.to_csv(out, sep="\t", index=False)
    print(f"  ✓ Saved {out}")

def consolidate_plasme(pred_dir, out_dir, tech, assembler):
    base_dir = Path(pred_dir) / f"plasme/{tech}.{assembler}"
    files = list(base_dir.glob("*_plasmids.fasta_report.csv"))
    if not files:
        print(f"[WARN] No PLASMe outputs for {tech}.{assembler}")
        return

    print(f"[INFO] Processing PLASMe for {tech}.{assembler}...")
    dfs = []
    for f in files:
        df = pd.read_csv(f)
        df["sample"] = f"{tech}.{assembler}"
        dfs.append(df.iloc[:, :6])
    out = Path(out_dir) / f"{tech}.{assembler}/plasme_{tech}_{assembler}_summary.csv"
    ensure_dir(out.parent)
    pd.concat(dfs, ignore_index=True).to_csv(out, index=False)
    print(f"  ✓ Saved {out}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pred_dir", required=True, help="Path to 03_predictions folder")
    parser.add_argument("--out_dir", required=True, help="Path to 04_consolidated folder")
    parser.add_argument("--tech", required=True)
    parser.add_argument("--assembler", required=True)
    args = parser.parse_args()

    consolidate_genomad(args.pred_dir, args.out_dir, args.tech, args.assembler)
    consolidate_vibrant(args.pred_dir, args.out_dir, args.tech, args.assembler)
    consolidate_virsorter2(args.pred_dir, args.out_dir, args.tech, args.assembler)
    consolidate_deepmc(args.pred_dir, args.out_dir, args.tech, args.assembler)
    consolidate_phamer(args.pred_dir, args.out_dir, args.tech, args.assembler)
    consolidate_plasme(args.pred_dir, args.out_dir, args.tech, args.assembler)

    print(f"\n[✓] Finished consolidating all tools for {args.tech}.{args.assembler}")
    print(f"Results saved in: {args.out_dir}/{args.tech}.{args.assembler}/")

if __name__ == "__main__":
    main()
