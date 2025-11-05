#!/usr/bin/env python3
"""
Create score matrix from tool predictions.

Usage:
    python 03_score_matrix.py \
        --cons_dir /path/to/04_consolidated \
        --tech illumina --assembler megahit
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse

# ---------- Loader helpers ----------
def load_virsorter_group(path, group, colname):
    df = pd.read_csv(path, sep="\t")
    if "max_score_group" not in df.columns or "seqname" not in df.columns:
        return pd.DataFrame(columns=["Contig_ID", colname])
    df = df[df["max_score_group"] == group].copy()
    df["Contig_ID"] = df["seqname"]
    return df[["Contig_ID", "max_score"]].rename(columns={"max_score": colname})

def load_phamer(path):
    df = pd.read_csv(path, sep="\t")
    if not {"Accession", "PhaMerScore", "Proportion"}.issubset(df.columns):
        return pd.DataFrame(columns=["Contig_ID", "Phabox"])
    df["Phabox"] = df["PhaMerScore"] * df["Proportion"]
    df = df[df["Phabox"] > 0.7].copy()
    return df[["Accession", "Phabox"]].rename(columns={"Accession": "Contig_ID"})

def load_genomad_virus(path):
    df = pd.read_csv(path, sep="\t")
    if not {"seq_name", "virus_score"}.issubset(df.columns):
        return pd.DataFrame(columns=["Contig_ID", "GenomadPH"])
    df["Contig_ID"] = df["seq_name"]
    return df[["Contig_ID", "virus_score"]].rename(columns={"virus_score": "GenomadPH"})

def load_genomad_plasmid(path):
    df = pd.read_csv(path, sep="\t")
    if not {"seq_name", "plasmid_score"}.issubset(df.columns):
        return pd.DataFrame(columns=["Contig_ID", "GenomadPL"])
    df["Contig_ID"] = df["seq_name"]
    return df[["Contig_ID", "plasmid_score"]].rename(columns={"plasmid_score": "GenomadPL"})

def load_plasme(path):
    df = pd.read_csv(path)
    if not {"contig", "score"}.issubset(df.columns):
        return pd.DataFrame(columns=["Contig_ID", "Plasme"])
    df = df[df["score"] >= 0.6].copy()
    return df[["contig", "score"]].rename(columns={"contig": "Contig_ID", "score": "Plasme"})

def load_deepmc(path):
    df = pd.read_csv(path, sep="\t")
    needed = {"best_choice", "ProkaryoteVirus", "Plasmid", "Sequence Name"}
    if not needed.issubset(df.columns):
        return pd.DataFrame(columns=["Contig_ID", "DeepmcPH", "DeepmcPL"])

    df["DeepmcPH"] = 0.0
    df["DeepmcPL"] = 0.0

    ph_mask = (df["best_choice"] == 5) & (df["ProkaryoteVirus"] > 0.4)
    pl_mask = (df["best_choice"] == 3) & (df["Plasmid"] > 0.6)
    df.loc[ph_mask, "DeepmcPH"] = df.loc[ph_mask, "ProkaryoteVirus"]
    df.loc[pl_mask, "DeepmcPL"] = df.loc[pl_mask, "Plasmid"]

    df["Contig_ID"] = df["Sequence Name"]
    return df[["Contig_ID", "DeepmcPH", "DeepmcPL"]]

def load_vibrant(path):
    df = pd.read_csv(path, sep="\t")
    if not {"scaffold", "all VOG"}.issubset(df.columns):
        return pd.DataFrame(columns=["Contig_ID", "Vibrant"])
    df["Contig_ID"] = df["scaffold"]
    return df[["Contig_ID", "all VOG"]].rename(columns={"all VOG": "Vibrant"})

# ---------- Classification ----------
def classify_consensus(df):
    phage_cols   = ["Phabox", "GenomadPH", "DeepmcPH", "VirsorterDS", "Vibrant"]
    plasmid_cols = ["GenomadPL", "DeepmcPL", "Plasme"]
    other_cols   = ["VirsorterSS", "VirsorterRNA", "VirsorterNCLDV", "VirsorterLAV"]

    for c in phage_cols + plasmid_cols + other_cols:
        if c not in df.columns:
            df[c] = 0.0

    df["phage_count"]   = (df[phage_cols]   > 0).sum(axis=1)
    df["plasmid_count"] = (df[plasmid_cols] > 0).sum(axis=1)
    df["other_count"]   = (df[other_cols]   > 0).sum(axis=1)

    def classify(row):
        ph, pl, ot = row["phage_count"], row["plasmid_count"], row["other_count"]
        if ph >= 1 and pl == 0 and ot == 0:   return "phage"
        if ph == 0 and pl >= 1 and ot == 0:   return "plasmid"
        if ph >= 1 and pl >= 1 and ot == 0:   return "PP"
        if ph >= 1 and pl >= 1 and ot >= 1:   return "PPV"
        if ph == 0 and pl == 0 and ot >= 1:   return "virus"
        return "uncertain"

    df["Consensus_Label"] = df.apply(classify, axis=1)
    return df

def deduplicate(df):
    score_cols = [c for c in df.columns if c not in ["Contig_ID", "phage_count", "plasmid_count", "other_count", "Consensus_Label"]]
    agg = {c: 'max' for c in score_cols}
    df2 = df.groupby("Contig_ID", as_index=False).agg(agg)
    df2 = classify_consensus(df2)
    return df2

# --- Resolve parent vs prophage conflicts ---
def resolve_parent_prophage_conflicts(df):
    # Map child prophages to parents
    df["Parent_ID"] = df["Contig_ID"].apply(lambda x: x.split("_prophage_")[0] if "_prophage_" in x else x)

    to_drop = set()
    for parent, group in df.groupby("Parent_ID"):
        # Find main contig
        parent_row = group[group["Contig_ID"] == parent]
        if parent_row.empty:
            continue

        parent_label = parent_row["Consensus_Label"].iloc[0]
        parent_phage_count = parent_row["phage_count"].iloc[0]

        # Find prophage fragments
        prophages = group[group["Contig_ID"].str.contains("_prophage_")]
        for _, row in prophages.iterrows():
            if (
                parent_label in ["phage", "PP", "PPV"]
                and row["phage_count"] <= parent_phage_count
            ):
                to_drop.add(row["Contig_ID"])

    # Drop redundant prophages
    df = df[~df["Contig_ID"].isin(to_drop)].copy()
    df.drop(columns=["Parent_ID"], inplace=True)
    return df

# ---------- Main ----------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cons_dir", required=True, help="Path to 04_consolidated folder")
    parser.add_argument("--tech", required=True, help="Technology: illumina, pacbio, or ont")
    parser.add_argument("--assembler", required=True, help="Assembler: megahit, spades, flye, hifiasm, autocycler")
    args = parser.parse_args()

    combo = f"{args.tech}.{args.assembler}"
    combo_dir = Path(args.cons_dir) / combo
    
    if not combo_dir.exists():
        print(f"[WARN] Directory not found: {combo_dir}")
        return

    print(f"[INFO] Processing score matrix for {combo}...")

    files = {
        "virsorter":   next(iter(combo_dir.glob("virsorter2_*_summary_homogenized.tsv")), None),
        "phamer":      next(iter(combo_dir.glob("phamer_*_summary.tsv")), None),
        "genomad_v":   next(iter(combo_dir.glob("genomad_*_virus_summary_homogenized.tsv")), None),
        "genomad_p":   next(iter(combo_dir.glob("genomad_*_plasmid_summary.tsv")), None),
        "deepmc":      next(iter(combo_dir.glob("deepmicroclass_*_softmax.tsv")), None),
        "plasme":      next(iter(combo_dir.glob("plasme_*_summary.csv")), None),
        "vibrant":     next(iter(combo_dir.glob("vibrant_*_summary_homogenized.tsv")), None),
    }

    dfs = []

    if files["phamer"]: dfs.append(load_phamer(files["phamer"]))
    if files["genomad_v"]: dfs.append(load_genomad_virus(files["genomad_v"]))
    if files["virsorter"]:
        dfs.append(load_virsorter_group(files["virsorter"], "dsDNAphage", "VirsorterDS"))
        dfs.append(load_virsorter_group(files["virsorter"], "ssDNA", "VirsorterSS"))
        dfs.append(load_virsorter_group(files["virsorter"], "RNA", "VirsorterRNA"))
        dfs.append(load_virsorter_group(files["virsorter"], "NCLDV", "VirsorterNCLDV"))
        dfs.append(load_virsorter_group(files["virsorter"], "lavidaviridae", "VirsorterLAV"))
    if files["vibrant"]: dfs.append(load_vibrant(files["vibrant"]))
    if files["deepmc"]: dfs.append(load_deepmc(files["deepmc"]))
    if files["genomad_p"]: dfs.append(load_genomad_plasmid(files["genomad_p"]))
    if files["plasme"]: dfs.append(load_plasme(files["plasme"]))

    if not dfs:
        print(f"  ⚠️ No input files found")
        return

    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on="Contig_ID", how="outer")

    merged = merged.fillna(0)
    merged = merged.round(3)
    merged = classify_consensus(merged)
    merged = resolve_parent_prophage_conflicts(merged)

    cols = [c for c in merged.columns if c != "Contig_ID"]
    merged = merged[["Contig_ID"] + cols]

    out_tsv = combo_dir / f"{combo}_merged_tool_scores.tsv"

    merged.to_csv(out_tsv, sep="\t", index=False)

    print(f"  ✅ Saved: {out_tsv.name}")
    print("  Consensus label counts:")
    for label, count in merged["Consensus_Label"].value_counts().items():
        print(f"    {label}: {count}")
    
    print(f"[✓] Completed {combo}")

if __name__ == "__main__":
    main()