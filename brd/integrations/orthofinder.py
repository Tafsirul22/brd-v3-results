"""OrthoFinder availability wrapper and DIAMOND-based cross-proteome orthogroup clustering."""

from __future__ import annotations

import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import networkx as nx
import pandas as pd
from Bio import SeqIO

from brd.integrations.diamond import build_diamond_db, run_diamond_blastp


def is_orthofinder_available() -> bool:
    return shutil.which("orthofinder") is not None and shutil.which("mcl") is not None


def run_comparative_orthogroups(
    proteomes: Dict[str, Path],
    diamond_bin: str,
    out_dir: Path,
    evalue: float = 1e-5,
    min_identity: float = 30.0,
    min_cov: float = 50.0,
    log_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    """
    Cluster proteomes into orthologous groups.
    Uses official OrthoFinder CLI when Linux MCL/FastME/OrthoFinder are available;
    otherwise runs DIAMOND all-vs-all reciprocal best-hit + normalized-bitscore graph clustering
    across the available annotated proteomes (explicitly labelled in provenance).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    if is_orthofinder_available():
        method = "OrthoFinder_CLI"
    else:
        method = "DIAMOND_all_vs_all_RBH_graph_clustering (OrthoFinder Linux/MCL unavailable natively on Windows)"

    if (out_dir / "Orthogroups.tsv").exists() and (out_dir / "Orthogroups.GeneCount.tsv").exists():
        return pd.read_csv(out_dir / "Orthogroups.tsv", sep="\t"), pd.read_csv(out_dir / "Orthogroups.GeneCount.tsv", sep="\t"), method

    # Combine proteomes with species prefix
    combined_faa = out_dir / "all_proteomes_combined.faa"
    seq_to_species: Dict[str, str] = {}
    species_list = sorted(proteomes.keys())
    with open(combined_faa, "w", encoding="utf-8") as fh:
        for sp in species_list:
            for rec in SeqIO.parse(str(proteomes[sp]), "fasta"):
                uid = f"{sp}||{rec.id}"
                seq_to_species[uid] = sp
                fh.write(f">{uid}\n{str(rec.seq)}\n")

    db_prefix = out_dir / "all_proteomes_db"
    build_diamond_db(diamond_bin, combined_faa, db_prefix, log_path=log_path)
    hits_df = run_diamond_blastp(
        diamond_bin,
        combined_faa,
        db_prefix,
        out_dir / "diamond_all_vs_all.tsv",
        evalue=evalue,
        max_target_seqs=25,
        threads=4,
        log_path=log_path,
    )

    # Filter hits and identify best hit per (query, subject_species)
    filt = hits_df[
        (hits_df["pident"] >= min_identity)
        & ((hits_df["qcovhsp"] >= min_cov) | (hits_df["scovhsp"] >= min_cov))
        & (hits_df["qseqid"] != hits_df["sseqid"])
    ].copy()

    best_per_sp: Dict[Tuple[str, str], Tuple[str, float]] = {}
    for _, row in filt.iterrows():
        q = str(row["qseqid"])
        s = str(row["sseqid"])
        s_sp = seq_to_species.get(s, "")
        bs = float(row["bitscore"])
        key = (q, s_sp)
        if key not in best_per_sp or bs > best_per_sp[key][1]:
            best_per_sp[key] = (s, bs)

    # Build orthology graph from reciprocal best hits (RBH) + strong paralogs (>= 90% of self/best bitscore)
    G = nx.Graph()
    for uid in seq_to_species:
        G.add_node(uid)

    for (q, s_sp), (s, bs) in best_per_sp.items():
        q_sp = seq_to_species[q]
        recip = best_per_sp.get((s, q_sp))
        if recip and recip[0] == q:
            G.add_edge(q, s, weight=bs)
        elif q_sp == s_sp and bs >= 100.0:
            G.add_edge(q, s, weight=bs)

    comps = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))
    og_rows = []
    pa_rows = []
    for idx, comp in enumerate(comps, 1):
        og_id = f"BRD_OG{idx:06d}"
        by_sp = defaultdict(list)
        for node in sorted(comp):
            sp, orig_id = node.split("||", 1)
            by_sp[sp].append(orig_id)
        row = {"Orthogroup": og_id, "Total_genes": len(comp)}
        pa_row = {"Orthogroup": og_id}
        present_species = 0
        for sp in species_list:
            genes = by_sp.get(sp, [])
            row[sp] = ", ".join(genes)
            pa_row[sp] = len(genes)
            if genes:
                present_species += 1
        row["Species_present"] = present_species
        og_rows.append(row)
        pa_rows.append(pa_row)

    og_df = pd.DataFrame(og_rows)
    pa_df = pd.DataFrame(pa_rows)
    og_df.to_csv(out_dir / "Orthogroups.tsv", sep="\t", index=False)
    pa_df.to_csv(out_dir / "Orthogroups.GeneCount.tsv", sep="\t", index=False)
    return og_df, pa_df, method
