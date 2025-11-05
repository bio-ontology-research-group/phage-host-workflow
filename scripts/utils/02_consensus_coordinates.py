#!/usr/bin/env python3
"""
Script to homogenize prophage prediction results using consensus coordinates.

Usage:
    python 02_consensus_coordinates.py \
        --cons_dir /path/to/04_consolidated \
        --tech illumina --assembler megahit
"""

from pathlib import Path
import pandas as pd
import numpy as np
import re
import argparse

# ------------------------- Parsing -------------------------

PROPHAGE_RE = re.compile(r"^(?P<contig>.+?)_prophage_(?P<start>\d+)_(?P<end>\d+)$")

def parse_prophage_id(s: str):
    m = PROPHAGE_RE.match(str(s))
    if not m:
        return None, None, None
    contig = m.group("contig")
    a, b = int(m.group("start")), int(m.group("end"))
    return contig, min(a, b), max(a, b)

def load_from_virsorter(fp: Path):
    rows = []
    if not (fp and fp.exists()):
        return pd.DataFrame(rows)
    df = pd.read_csv(fp, sep="\t", dtype=str).fillna("")
    if "seqname" not in df.columns:
        return pd.DataFrame(rows)
    for _, r in df[["seqname"]].drop_duplicates().iterrows():
        contig, s, e = parse_prophage_id(r["seqname"])
        if contig:
            rows.append({"tool": "VirSorter2", "contig": contig, "start": s, "end": e, "id_raw": r["seqname"]})
    return pd.DataFrame(rows)

def load_from_vibrant(fp: Path):
    rows = []
    if not (fp and fp.exists()):
        return pd.DataFrame(rows)
    df = pd.read_csv(fp, sep="\t", dtype=str).fillna("")
    if "scaffold" not in df.columns:
        return pd.DataFrame(rows)
    for _, r in df[["scaffold"]].drop_duplicates().iterrows():
        contig, s, e = parse_prophage_id(r["scaffold"])
        if contig:
            rows.append({"tool": "VIBRANT", "contig": contig, "start": s, "end": e, "id_raw": r["scaffold"]})
    return pd.DataFrame(rows)

def load_from_genomad(fp: Path):
    rows = []
    if not (fp and fp.exists()):
        return pd.DataFrame(rows)
    df = pd.read_csv(fp, sep="\t", dtype=str).fillna("")
    name_col = "seq_name" if "seq_name" in df.columns else ("seqname" if "seqname" in df.columns else None)
    if not name_col:
        return pd.DataFrame(rows)
    for _, r in df[[name_col]].drop_duplicates().iterrows():
        contig, s, e = parse_prophage_id(r[name_col])
        if contig:
            rows.append({"tool": "GeNomad", "contig": contig, "start": s, "end": e, "id_raw": r[name_col]})
    return pd.DataFrame(rows)

# ------------------------- Clustering -------------------------

def cluster_by_contig(calls, overlap_thresh=0.5, max_dist=1000):
    member_rows = []
    for contig, dfc in calls.groupby("contig", sort=False):
        dfc = dfc.sort_values(["start","end"]).reset_index(drop=True)
        cluster_id = 0
        cur_start = cur_end = None
        cur_members = []

        def flush():
            nonlocal cur_members, cur_start, cur_end, cluster_id
            if cur_members:
                sub = dfc.loc[cur_members].copy()
                sub["cluster_idx"] = cluster_id
                sub["contig"] = contig
                member_rows.append(sub)

        for i, r in dfc.iterrows():
            s, e = int(r["start"]), int(r["end"])
            if cur_start is None:
                cluster_id += 1
                cur_start, cur_end = s, e
                cur_members = [i]
                continue
            inter = max(0, min(e, cur_end) - max(s, cur_start))
            len_a = e - s
            len_b = cur_end - cur_start
            min_recip = 0 if (len_a==0 or len_b==0) else min(inter/len_a, inter/len_b)
            gap = 0 if inter>0 else max(s - cur_end, cur_start - e, 0)
            if (min_recip >= overlap_thresh) or (gap <= max_dist):
                cur_members.append(i)
                cur_start = min(cur_start, s)
                cur_end   = max(cur_end, e)
            else:
                flush()
                cluster_id += 1
                cur_start, cur_end = s, e
                cur_members = [i]
        flush()

    if not member_rows:
        return pd.DataFrame(columns=list(calls.columns) + ["cluster_idx"])
    return pd.concat(member_rows, ignore_index=True)

# ------------------------- Consensus -------------------------

def best_group(values, tools, tol):
    n = len(values)
    adj = {i:set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if abs(values[i] - values[j]) <= tol:
                adj[i].add(j); adj[j].add(i)

    seen = set()
    groups = []
    for i in range(n):
        if i in seen: continue
        stack = [i]
        comp = set()
        while stack:
            u = stack.pop()
            if u in seen: continue
            seen.add(u); comp.add(u)
            for v in adj[u]:
                if v not in seen:
                    stack.append(v)
        if len(comp) >= 2:
            groups.append(comp)

    if not groups:
        return None, []

    best = None
    best_key = None
    for comp in groups:
        idx = sorted(list(comp))
        vals = [values[k] for k in idx]
        width = max(vals) - min(vals)
        key = (-len(idx), width)
        if best_key is None or key < best_key:
            best_key = key
            best = idx

    med = int(np.median([values[k] for k in best]))
    used = [tools[k] for k in best]
    return med, used

def consensus_edgewise(members_df, tol=200):
    rows = []
    rep = (members_df.groupby(["contig","cluster_idx","tool"], as_index=False)
       .agg(start=("start", lambda x: int(np.median(x))),
            end=("end",   lambda x: int(np.median(x)))))

    for (contig, cid), g in rep.groupby(["contig","cluster_idx"], sort=False):
        tools  = list(g["tool"])
        starts = g["start"].to_numpy()
        ends   = g["end"].to_numpy()
        n_tools = len(tools)

        if n_tools == 1:
            s_cons = int(starts[0]); e_cons = int(ends[0])
            method = "single_tool"
            used_s = tools; used_e = tools
        else:
            s_group, used_s = best_group(starts, tools, tol)
            if s_group is None:
                s_cons = int(starts.min())
                used_s = tools
                start_part = "start_min"
            else:
                s_cons = int(s_group)
                start_part = f"start_consensus_{len(used_s)}"

            e_group, used_e = best_group(ends, tools, tol)
            if e_group is None:
                e_cons = int(ends.max())
                used_e = tools
                end_part = "end_max"
            else:
                e_cons = int(e_group)
                end_part = f"end_consensus_{len(used_e)}"

            method = f"{start_part}+{end_part}"

        consensus_id = f"{contig}_prophage_{s_cons}_{e_cons}"

        rows.append({
            "consensus_id": consensus_id,
            "contig": contig,
            "n_tools": int(n_tools),
            "tools": ",".join(sorted(tools)),
            "consensus_start": int(s_cons),
            "consensus_end": int(e_cons),
            "method": method,
            "used_tools_start": ",".join(sorted(used_s)),
            "used_tools_end": ",".join(sorted(used_e)),
        })
    out = pd.DataFrame(rows).sort_values(["contig","consensus_start"]).reset_index(drop=True)
    return out

# ------------------------- Homogenization -------------------------

def create_mapping(consensus_df, members_df):
    mapping = {}
    consensus_with_clusters = consensus_df.copy()
    consensus_with_clusters["cluster_idx"] = range(1, len(consensus_df) + 1)
    
    for contig in members_df["contig"].unique():
        contig_members = members_df[members_df["contig"] == contig]
        contig_consensus = consensus_with_clusters[consensus_with_clusters["contig"] == contig]
        
        for cluster_idx in contig_members["cluster_idx"].unique():
            cluster_members = contig_members[contig_members["cluster_idx"] == cluster_idx]
            
            if len(contig_consensus) >= cluster_idx:
                consensus_id = contig_consensus.iloc[cluster_idx - 1]["consensus_id"]
                
                for _, member in cluster_members.iterrows():
                    mapping[member["id_raw"]] = consensus_id
    
    return mapping

def homogenize_file(input_path, output_path, mapping, id_column):
    if not input_path.exists():
        print(f"  ⚠️ {input_path.name} does not exist, skipping")
        return
    
    df = pd.read_csv(input_path, sep="\t")
    
    if id_column not in df.columns:
        print(f"  ⚠️ Column '{id_column}' not found, skipping")
        return
    
    df[f"{id_column}_consensus"] = df[id_column].map(mapping).fillna(df[id_column])
    df[id_column] = df[f"{id_column}_consensus"]
    df.drop(f"{id_column}_consensus", axis=1, inplace=True)
    
    df.to_csv(output_path, sep="\t", index=False)
    print(f"  ✅ Homogenized {input_path.name}")

# ------------------------- Main -------------------------

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

    print(f"[INFO] Processing consensus for {combo}...")

    # Find files
    virsorter_files = list(combo_dir.glob("virsorter2_*_summary.tsv"))
    vibrant_files = list(combo_dir.glob("vibrant_*_summary.tsv"))
    genomad_files = list(combo_dir.glob("genomad_*_virus_summary.tsv"))
    
    VIRSORTER_PATH = virsorter_files[0] if virsorter_files else None
    VIBRANT_PATH = vibrant_files[0] if vibrant_files else None
    GENOMAD_PATH = genomad_files[0] if genomad_files else None

    if not any([VIRSORTER_PATH, VIBRANT_PATH, GENOMAD_PATH]):
        print(f"  ⚠️ No prophage prediction files found")
        return

    # Load calls
    parts = []
    if VIRSORTER_PATH: parts.append(load_from_virsorter(VIRSORTER_PATH))
    if VIBRANT_PATH: parts.append(load_from_vibrant(VIBRANT_PATH))
    if GENOMAD_PATH: parts.append(load_from_genomad(GENOMAD_PATH))
    
    parts = [x for x in parts if not x.empty]
    if not parts:
        print(f"  ⚠️ No prophage calls parsed")
        return
    
    calls = pd.concat(parts, ignore_index=True)
    calls["start"] = calls["start"].astype(int)
    calls["end"] = calls["end"].astype(int)
    print(f"  Loaded {len(calls)} prophage calls")

    # Cluster
    members = cluster_by_contig(calls)
    if members.empty:
        print(f"  ⚠️ No clusters found")
        return
    print(f"  Created {members['cluster_idx'].nunique()} clusters")

    # Consensus
    consensus = consensus_edgewise(members)
    consensus_output = consensus[[
        "consensus_id", "contig", "n_tools", "tools",
        "consensus_start", "consensus_end", "method",
        "used_tools_start", "used_tools_end"
    ]]
    
    out_dir = Path(args.cons_dir) / "consensus_coordinates"
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / f"{combo}_prophage_consensus.tsv"
    consensus_output.to_csv(out_file, sep="\t", index=False)
    print(f"  ✅ Consensus: {out_file.name} ({len(consensus)} prophages)")

    # Mapping
    cluster_mapping = create_mapping(consensus, members)
    print(f"  Created {len(cluster_mapping)} mappings")

    # Homogenize
    print(f"  Homogenizing files...")
    if VIRSORTER_PATH:
        homogenize_file(VIRSORTER_PATH, 
                       combo_dir / f"{VIRSORTER_PATH.stem}_homogenized.tsv",
                       cluster_mapping, "seqname")
    
    if VIBRANT_PATH:
        homogenize_file(VIBRANT_PATH,
                       combo_dir / f"{VIBRANT_PATH.stem}_homogenized.tsv",
                       cluster_mapping, "scaffold")
    
    if GENOMAD_PATH:
        genomad_id_col = "seq_name" if "seq_name" in pd.read_csv(GENOMAD_PATH, sep="\t", nrows=1).columns else "seqname"
        homogenize_file(GENOMAD_PATH,
                       combo_dir / f"{GENOMAD_PATH.stem}_homogenized.tsv",
                       cluster_mapping, genomad_id_col)
    
    print(f"[✓] Completed {combo}")

if __name__ == "__main__":
    main()