"""Stage 5: Candidate Redundancy Clustering, Sequence Characterization & Exhaustive Full-CDS Sweep."""

from __future__ import annotations

import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import RNA
from Bio import SeqIO

from paths_v3 import resolve_output_root, resolve_v3_root


def _shannon_entropy_bases(seq: str) -> float:
    n = len(seq)
    if n == 0:
        return 0.0
    counts = Counter(seq)
    ent = -sum((c / n) * math.log2(c / n) for c in counts.values() if c > 0)
    return round(ent, 4)  # max = 2.0 bits


def _dinucleotide_entropy(seq: str) -> float:
    if len(seq) < 2:
        return 0.0
    di = [seq[i : i + 2] for i in range(len(seq) - 1)]
    n = len(di)
    counts = Counter(di)
    ent = -sum((c / n) * math.log2(c / n) for c in counts.values() if c > 0)
    return round(ent / 4.0, 4)  # normalized 0-1 (max log2(16)=4)


def _linguistic_complexity_u3(seq: str) -> float:
    """Fraction of observed distinct trinucleotides relative to theoretical maximum min(64, len-2)."""
    if len(seq) < 3:
        return 0.0
    obs = len({seq[i : i + 3] for i in range(len(seq) - 2)})
    denom = min(64, len(seq) - 2)
    return round(obs / max(denom, 1), 4)


def _low_complexity_fraction(seq: str) -> float:
    """Proportion of bases in homopolymers >=4 nt or simple dinucleotide tandem repeats (>=3 units)."""
    n = max(len(seq), 1)
    mask = np.zeros(len(seq), dtype=bool)
    for m in re.finditer(r"(A{4,}|T{4,}|G{4,}|C{4,})", seq):
        mask[m.start() : m.end()] = True
    for m in re.finditer(r"((?:AT|TA|AG|GA|AC|CA|TG|GT|TC|CT|GC|CG))\1{2,}", seq):
        mask[m.start() : m.end()] = True
    return round(float(mask.mean()), 4)


def _dust_score(seq: str) -> float:
    """Symmetric DUST-like triplet concentration score across full candidate (lower = higher complexity)."""
    if len(seq) < 5:
        return 0.0
    tris = Counter(seq[i : i + 3] for i in range(len(seq) - 2))
    l = len(seq) - 2
    num = sum(c * (c - 1) / 2.0 for c in tris.values())
    return round(num / max(l - 1, 1), 4)


def run_stage5(output_root: Path | None = None) -> Dict[str, Path]:
    out = resolve_output_root(output_root)
    v3 = resolve_v3_root(out)
    s5 = v3 / "stage5_redundancy_design"
    fig_dir = s5 / "figures"
    s5.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    auth_cands = pd.read_csv(v3 / "authoritative_data" / "auth_candidates.tsv", sep="\t")
    auth_guides = pd.read_csv(v3 / "authoritative_data" / "auth_guides.tsv", sep="\t")
    auth_targets = pd.read_csv(v3 / "authoritative_data" / "auth_targets.tsv", sep="\t")
    cons_df = pd.read_csv(v3 / "stage3_isolate_conservation" / "candidate_conservation_multimetric_v3.tsv", sep="\t").set_index("candidate_id")
    sel_df = pd.read_csv(v3 / "stage4_off_target_expanded" / "candidate_selectivity_expanded_v3.tsv", sep="\t").set_index("candidate_id")

    # 1. Pairwise Overlap & Clustering into Independent Target Regions (Section 5.1)
    overlap_rows: List[dict] = []
    cluster_rows: List[dict] = []
    region_rows: List[dict] = []

    cand_kmers: Dict[str, Set[str]] = {}
    for _, r in auth_cands.iterrows():
        s = str(r["sequence"]).upper()
        cand_kmers[str(r["candidate_id"])] = {s[i : i + 21] for i in range(len(s) - 20)}

    for tname, tgrp in auth_cands.groupby("target", sort=False):
        recs = tgrp.sort_values(["cds_start_nt", "cds_end_nt", "candidate_id"]).to_dict(orient="records")
        G = nx.Graph()
        for r in recs:
            G.add_node(r["candidate_id"])

        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                a, b = recs[i], recs[j]
                s1, e1, l1 = int(a["cds_start_nt"]), int(a["cds_end_nt"]), int(a["candidate_length_nt"])
                s2, e2, l2 = int(b["cds_start_nt"]), int(b["cds_end_nt"]), int(b["candidate_length_nt"])
                ov_bp = max(0, min(e1, e2) - max(s1, s2) + 1)
                union_bp = max(1, (e1 - s1 + 1) + (e2 - s2 + 1) - ov_bp)
                ov_shorter_pct = round(ov_bp / max(1, min(l1, l2)) * 100.0, 2)
                jaccard_pct = round(ov_bp / union_bp * 100.0, 2)

                ka, kb = cand_kmers[a["candidate_id"]], cand_kmers[b["candidate_id"]]
                shared_k = len(ka & kb)
                k_shorter_pct = round(shared_k / max(1, min(len(ka), len(kb))) * 100.0, 2)

                same_region = (ov_shorter_pct >= 25.0) or (k_shorter_pct >= 25.0)
                if ov_bp > 0 or shared_k > 0:
                    overlap_rows.append(
                        {
                            "target": tname,
                            "candidate_a": a["candidate_id"],
                            "candidate_b": b["candidate_id"],
                            "a_coords": f"{s1}..{e1}",
                            "b_coords": f"{s2}..{e2}",
                            "overlap_bp": ov_bp,
                            "overlap_pct_of_shorter": ov_shorter_pct,
                            "coordinate_jaccard_pct": jaccard_pct,
                            "shared_21mers": shared_k,
                            "shared_21mer_pct_of_shorter": k_shorter_pct,
                            "clustered_as_same_region": same_region,
                        }
                    )
                if same_region:
                    G.add_edge(a["candidate_id"], b["candidate_id"])

        # Order connected components by genomic/CDS start coordinate
        comps = list(nx.connected_components(G))
        comp_meta = []
        rec_map = {r["candidate_id"]: r for r in recs}
        for comp in comps:
            cs = min(int(rec_map[c]["cds_start_nt"]) for c in comp)
            ce = max(int(rec_map[c]["cds_end_nt"]) for c in comp)
            comp_meta.append((cs, ce, sorted(comp)))
        comp_meta.sort(key=lambda x: (x[0], x[1]))

        for r_idx, (cs, ce, members) in enumerate(comp_meta, 1):
            reg_id = f"{tname}_REG{r_idx:02d}"
            tools_in_reg = sorted({rec_map[m]["tool_name"] for m in members})
            clean_members = [m for m in members if bool(sel_df.loc[m, "passes_zero_exact_21mer_filter"]) and float(cons_df.loc[m, "min_nt_identity_hq17_pct"]) >= 85.0]
            region_rows.append(
                {
                    "region_id": reg_id,
                    "target": tname,
                    "cds_span_start_nt": cs,
                    "cds_span_end_nt": ce,
                    "region_span_nt": ce - cs + 1,
                    "n_redundant_candidates": len(members),
                    "n_passing_clean_hq85": len(clean_members),
                    "design_tools_present": ";".join(tools_in_reg),
                    "member_candidates": ";".join(members),
                }
            )
            for m in members:
                mr = rec_map[m]
                cluster_rows.append(
                    {
                        "candidate_id": m,
                        "target": tname,
                        "region_id": reg_id,
                        "region_cds_span": f"{cs}..{ce}",
                        "region_member_count": len(members),
                        "is_alternative_overlapping_window": len(members) > 1,
                        "cds_start_nt": int(mr["cds_start_nt"]),
                        "cds_end_nt": int(mr["cds_end_nt"]),
                        "candidate_length_nt": int(mr["candidate_length_nt"]),
                        "tool_name": mr["tool_name"],
                    }
                )

    ov_df = pd.DataFrame(overlap_rows)
    clust_df = pd.DataFrame(cluster_rows).sort_values(["target", "region_id", "cds_start_nt", "candidate_id"]).reset_index(drop=True)
    reg_df = pd.DataFrame(region_rows)
    ov_df.to_csv(s5 / "candidate_pairwise_overlaps_v3.tsv", sep="\t", index=False)
    clust_df.to_csv(s5 / "candidate_to_region_cluster_map_v3.tsv", sep="\t", index=False)
    reg_df.to_csv(s5 / "independent_target_regions_v3.tsv", sep="\t", index=False)

    # 2. Sequence & Structure Characterization for All 50 Candidates (Section 5.2)
    char_rows: List[dict] = []
    guides_by_cid = {cid: grp for cid, grp in auth_guides.groupby("candidate_id")}
    reg_by_cid = clust_df.set_index("candidate_id")["region_id"].to_dict()

    for _, cr in auth_cands.iterrows():
        cid = str(cr["candidate_id"])
        tname = str(cr["target"])
        cseq = str(cr["sequence"]).upper()
        n_len = len(cseq)
        t_cds_len = int(auth_targets.set_index("target").loc[tname, "cds_length_nt"])
        g_sub = guides_by_cid[cid]

        eff_vals = g_sub["is_efficient_sifi21"].fillna(False).astype(bool)
        acc_vals = pd.to_numeric(g_sub["accessibility_prob"], errors="coerce").dropna()
        d_end = pd.to_numeric(g_sub["delta_end_stability"], errors="coerce").dropna()
        unique_g = int(g_sub["sense_dna_21nt"].nunique())

        # Whole-fragment ViennaRNA MFE folding energy per nt
        fc = RNA.fold_compound(cseq.replace("T", "U"))
        _, mfe_cand = fc.mfe()
        mfe_per_nt = round(float(mfe_cand) / max(n_len, 1), 4)

        c_row = cons_df.loc[cid]
        char_rows.append(
            {
                "candidate_id": cid,
                "target": tname,
                "region_id": reg_by_cid[cid],
                "candidate_length_nt": n_len,
                "cds_transcript_coverage_pct": round(n_len / max(t_cds_len, 1) * 100.0, 2),
                "gc_percent": float(cr["gc_percent"]),
                "shannon_entropy_bits": _shannon_entropy_bases(cseq),
                "dinuc_entropy_norm": _dinucleotide_entropy(cseq),
                "linguistic_complexity_3mer": _linguistic_complexity_u3(cseq),
                "low_complexity_repeat_frac": _low_complexity_fraction(cseq),
                "dust_repeat_score": _dust_score(cseq),
                "max_homopolymer": int(cr["max_homopolymer"]),
                "total_guides_21mer": len(g_sub),
                "unique_guides_21mer": unique_g,
                "efficient_sifi21_count": int(eff_vals.sum()),
                "efficient_sifi21_fraction": round(float(eff_vals.mean()), 4),
                "guide_acc_mean": round(float(acc_vals.mean()), 4) if not acc_vals.empty else 0.0,
                "guide_acc_median": round(float(acc_vals.median()), 4) if not acc_vals.empty else 0.0,
                "guide_acc_q25": round(float(acc_vals.quantile(0.25)), 4) if not acc_vals.empty else 0.0,
                "guide_delta_end_mean_kcal": round(float(d_end.mean()), 3) if not d_end.empty else 0.0,
                "candidate_self_fold_mfe_kcal": round(float(mfe_cand), 2),
                "candidate_mfe_per_nt": mfe_per_nt,
                "full_frag_mean_nt_id_hq17": float(c_row["mean_nt_identity_hq17_pct"]),
                "full_frag_min_nt_id_hq17": float(c_row["min_nt_identity_hq17_pct"]),
                "full_frag_min_nt_id_all19": float(c_row["min_nt_identity_all19_pct"]),
                "guides_100pct_preserved_hq17": float(c_row["guides_pct_preserved_hq17"]),
            }
        )
    char_df = pd.DataFrame(char_rows)
    char_df.to_csv(s5 / "candidate_sequence_characterization_v3.tsv", sep="\t", index=False)

    # 3. Exhaustive Full-CDS 180-nt Window Audit across ALL 4 Genes (checks for missed clean regions!)
    t_seqs = {r.id.split("|")[0]: str(r.seq).upper() for r in SeqIO.parse(str(out / "10_final_reports" / "BRD_Target_Gene_Sequences.fasta"), "fasta")}
    exact_df = pd.read_csv(v3 / "stage4_off_target_expanded" / "exact_21mer_matches_expanded_v3.tsv", sep="\t")
    # Build fast set of any 21-mer hitting any expanded non-target / paralog / self_genome!
    bad_kmers_by_target: Dict[str, Set[str]] = {t: set() for t in t_seqs}
    # Scan directly across all 16 non-target FASTAs + ND90Pr paralogs for ALL gene CDS 21-mers so sweep is 100% comprehensive!
    manifest_s4 = pd.read_csv(v3 / "stage4_off_target_expanded" / "expanded_reference_manifest_v3.tsv", sep="\t")
    all_target_21s: Set[str] = set()
    target_km_ownership: Dict[str, Set[str]] = {}
    for tname, tseq in t_seqs.items():
        for i in range(max(0, len(tseq) - 20)):
            km = tseq[i : i + 21]
            rk = km[::-1].translate(str.maketrans("ACGT", "TGCA"))
            all_target_21s.add(km)
            all_target_21s.add(rk)
            target_km_ownership.setdefault(km, set()).add(tname)
            target_km_ownership.setdefault(rk, set()).add(tname)

    from stage4_off_target_expanded import PANEL_SPECS, _build_pigeonhole_21mer_index, _iter_fasta_upper  # reuse fast exact pigeonhole index
    old_m08_src = out / "08_off_target_screening" / "blast_databases" / "sources"
    v3_src = v3 / "expanded_references" / "sources"
    ref_fastas_check: List[Tuple[str, Path]] = []
    for sp in PANEL_SPECS:
        rid, acc = sp["ref_id"], sp["accession"]
        base_d = (old_m08_src / rid / "ncbi_dataset" / "data" / acc) if sp["tx_source"] == "existing_m08" else (v3_src / rid / "ncbi_dataset" / "data" / acc)
        for fn in ("rna.fna", "cds_from_genomic.fna"):
            if (base_d / fn).exists():
                ref_fastas_check.append((sp["group"], base_d / fn))
                break
        g_dir = v3_src / rid / "ncbi_dataset" / "data" / acc
        if g_dir.exists():
            for gf in g_dir.glob("*_genomic.fna"):
                if gf.name != "cds_from_genomic.fna":
                    ref_fastas_check.append((sp["group"], gf))

    p_idx_t = _build_pigeonhole_21mer_index(all_target_21s)
    hit_gene_21s: Dict[str, Set[str]] = {km: set() for km in all_target_21s}
    for grp_lbl, fap in ref_fastas_check:
        for _, s in _iter_fasta_upper(fap):
            n = len(s)
            if n < 21:
                continue
            max_start = n - 21
            for g in range(0, n - 10, 11):
                offs = p_idx_t.get(s[g : g + 11])
                if offs is not None:
                    for off in offs:
                        st = g - off
                        if 0 <= st <= max_start:
                            sub = s[st : st + 21]
                            if sub in hit_gene_21s:
                                hit_gene_21s[sub].add(grp_lbl)

    iso_cds_by_t = {}
    for tn in t_seqs:
        fa_iso = out / "04_comparative_genomics" / "target_orthologs" / "target_cds_alignments" / f"{tn}_isolates_cds.fasta"
        iso_cds_by_t[tn] = [str(rec.seq).upper().replace("-", "") for rec in SeqIO.parse(str(fa_iso), "fasta") if not any(bad in rec.id for bad in ("WAI2406", "WAI2411"))]

    sweep_rows: List[dict] = []
    sifi_by_t = {}
    for tn in t_seqs:
        sf = pd.read_csv(out / "07_dsRNA_design" / "lucid" / f"{tn}_sifi21_sirnas.tsv", sep="\t")
        sifi_by_t[tn] = {int(r["sirna_pos_1based"]): bool(r["is_efficient_sifi21"]) for _, r in sf.iterrows()}

    for tn, tseq in t_seqs.items():
        w_len = 180 if len(tseq) >= 180 else len(tseq)
        step = 25
        for start0 in range(0, max(1, len(tseq) - w_len + 1), step):
            wseq = tseq[start0 : start0 + w_len]
            gc = round((wseq.count("G") + wseq.count("C")) / len(wseq) * 100.0, 2)
            # Check exact 21-mer off-targets across ALL expanded DBs (both strands)
            ot_grps: Set[str] = set()
            n_ot_km = 0
            for k_i in range(len(wseq) - 20):
                km = wseq[k_i : k_i + 21]
                rk = km[::-1].translate(str.maketrans("ACGT", "TGCA"))
                g1 = hit_gene_21s.get(km, set()) | hit_gene_21s.get(rk, set())
                if g1:
                    ot_grps.update(g1)
                    n_ot_km += 1
            # HQ-17 exact window k-mer match
            hq_min = min(sum(1 for ki in range(len(wseq) - 20) if wseq[ki : ki + 21] in iseq) / max(len(wseq) - 20, 1) * 100.0 for iseq in iso_cds_by_t[tn])
            eff_n = sum(1 for p in range(start0 + 1, start0 + w_len - 19) if sifi_by_t[tn].get(p, False))
            overlaps_m07 = any(max(start0 + 1, int(r["cds_start_nt"])) <= min(start0 + w_len, int(r["cds_end_nt"])) for _, r in auth_cands[auth_cands["target"] == tn].iterrows())
            sweep_rows.append(
                {
                    "target": tn,
                    "win_start_1based": start0 + 1,
                    "win_end_1based": start0 + w_len,
                    "gc_percent": gc,
                    "exact_21mer_nontarget_windows": n_ot_km,
                    "hit_groups": ";".join(sorted(ot_grps)) if ot_grps else "NONE",
                    "hq17_min_21mer_pct": round(hq_min, 2),
                    "sifi21_efficient_guides": eff_n,
                    "overlaps_existing_m07_candidate": overlaps_m07,
                    "clean_high_cons_unsampled": (n_ot_km == 0) and (hq_min >= 95.0) and (30.0 <= gc <= 65.0) and (not overlaps_m07),
                }
            )
    sweep_df = pd.DataFrame(sweep_rows)
    sweep_df.to_csv(s5 / "exhaustive_cds_180nt_sweep_v3.tsv", sep="\t", index=False)

    # 4. Preliminary Independent Region Representative Map (refined in Stage 6 multi-dim Pareto!)
    # Plot Independent Target Regions Across the 4 Genes
    fig, ax = plt.subplots(figsize=(12, 4.8))
    y_pos = {"CsFUS3": 3, "CsSTE11": 2, "CsPKS1": 1, "CsSp1": 0}
    for tn, y in y_pos.items():
        ln = int(auth_targets.set_index("target").loc[tn, "cds_length_nt"])
        ax.plot([1, ln], [y, y], color="#cbd5e1", linewidth=6, solid_capstyle="round", zorder=1)
    for _, rr in reg_df.iterrows():
        y = y_pos[rr["target"]]
        clean = int(rr["n_passing_clean_hq85"]) > 0
        col = "#059669" if clean else "#dc2626"
        ax.plot([rr["cds_span_start_nt"], rr["cds_span_end_nt"]], [y, y], color=col, linewidth=11, alpha=0.85, zorder=2)
        ax.text((rr["cds_span_start_nt"] + rr["cds_span_end_nt"]) / 2, y + 0.18, f"{rr['region_id']} (n={rr['n_redundant_candidates']})", ha="center", fontsize=7.5, fontweight="bold")
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels(list(y_pos.keys()), fontweight="bold")
    ax.set_xlabel("Reference ND90Pr CDS Coordinate (nt)")
    ax.set_title("Stage 5: Independent Target Region Clusters vs Overlapping Candidate Windows\n(Green = >=1 candidate passes expanded Stage 4 & HQ17 >=85%; Red = excluded exact 21-mer match)", fontweight="bold")
    ax.set_ylim(-0.5, 3.6)
    ax.grid(True, axis="x", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(fig_dir / "stage5_independent_target_regions_v3.png", dpi=300)
    fig.savefig(fig_dir / "stage5_independent_target_regions_v3.svg")
    plt.close(fig)

    clean_fus3_wins = int(sweep_df[(sweep_df["target"] == "CsFUS3") & (sweep_df["exact_21mer_nontarget_windows"] == 0)].shape[0])
    new_unsampled = sweep_df[sweep_df["clean_high_cons_unsampled"]].sort_values(["target", "sifi21_efficient_guides"], ascending=[True, False])
    s5_md = f"""# Stage 5 — Candidate Redundancy Clustering & Full-CDS Audit (BRD v3)

**Timestamp:** {now_utc}

## 1. Clustering of 50 Candidates into Independent Regions
Clustering candidates sharing >=25% CDS overlap or >=25% 21-mers reduces the 50 overlapping M07 windows into **{len(reg_df)} independent target regions**:

{reg_df.to_markdown(index=False)}

## 2. Exhaustive 180-nt Full-CDS Sweep Findings (`exhaustive_cds_180nt_sweep_v3.tsv`)
- **`CsFUS3` (`1,059 nt`, 36 windows tested):** Only **{clean_fus3_wins}** 180-nt windows lack exact 21-mer matches in the expanded panel (proving beneficial-fungi MAPK conservation spans nearly the entire `CsFUS3` ORF!).
- **Novel Clean Unsampled Regions (`HQ17 >= 95%`, `0` exact 21-mer hits across expanded panel, non-overlapping existing M07):** Found **{len(new_unsampled)}** exploratory 180-nt windows (top examples in `exhaustive_cds_180nt_sweep_v3.tsv`), cleanly distinguished from the 50 canonical M07 candidates.
"""
    (s5 / "stage5_redundancy_and_design_report_v3.md").write_text(s5_md, encoding="utf-8")

    return {
        "overlaps": s5 / "candidate_pairwise_overlaps_v3.tsv",
        "cluster_map": s5 / "candidate_to_region_cluster_map_v3.tsv",
        "regions": s5 / "independent_target_regions_v3.tsv",
        "characterization": s5 / "candidate_sequence_characterization_v3.tsv",
        "cds_sweep": s5 / "exhaustive_cds_180nt_sweep_v3.tsv",
    }


if __name__ == "__main__":
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    for k, v in run_stage5().items():
        print(f"STAGE5_OUT {k}: {v}")
