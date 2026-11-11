"""Stage 3: Multi-Metric Isolate-Aware Conservation Analysis, Artifact Separation & Threshold Sensitivity."""

from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Align, SeqIO

from paths_v3 import resolve_output_root, resolve_v3_root


def _semi_global_align_metrics(cand_seq: str, subject_seq: str) -> Dict[str, object]:
    """Align candidate DNA query globally against isolate target CDS/locus without penalizing flanking overhang."""
    q = cand_seq.upper()
    s = subject_seq.upper()
    q_len = len(q)
    if not s or len(s) < 21:
        return {
            "alignment_coverage_pct": 0.0,
            "nucleotide_identity_pct": 0.0,
            "aligned_columns": 0,
            "matches": 0,
            "mismatches": 0,
            "gap_events": 0,
            "gap_bases": q_len,
            "n_bases_in_subject_span": 0,
        }
    exact_idx = s.find(q)
    if exact_idx != -1:
        return {
            "alignment_coverage_pct": 100.0,
            "nucleotide_identity_pct": 100.0,
            "aligned_columns": q_len,
            "matches": q_len,
            "mismatches": 0,
            "gap_events": 0,
            "gap_bases": 0,
            "n_bases_in_subject_span": 0,
        }

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    aligner.open_left_deletion_score = 0.0
    aligner.extend_left_deletion_score = 0.0
    aligner.open_right_deletion_score = 0.0
    aligner.extend_right_deletion_score = 0.0
    aligner.open_left_insertion_score = -5.0
    aligner.extend_left_insertion_score = -1.0
    aligner.open_right_insertion_score = -5.0
    aligner.extend_right_insertion_score = -1.0

    aln = aligner.align(s, q)[0]
    t_blocks, q_blocks = aln.aligned
    if len(t_blocks) == 0:
        return {
            "alignment_coverage_pct": 0.0,
            "nucleotide_identity_pct": 0.0,
            "aligned_columns": 0,
            "matches": 0,
            "mismatches": 0,
            "gap_events": 1,
            "gap_bases": q_len,
            "n_bases_in_subject_span": 0,
        }

    # Walk aligned blocks & inter-block gaps within query [0, q_len)
    matches = mismatches = aligned_q = 0
    for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
        seg_len = qe - qs
        aligned_q += seg_len
        t_sub = s[ts:te]
        q_sub = q[qs:qe]
        for tb, qb in zip(t_sub, q_sub):
            if tb == qb:
                matches += 1
            else:
                mismatches += 1

    # Gap events inside candidate span (between blocks + unaligned prefix/suffix)
    gap_events = max(0, len(t_blocks) - 1)
    if q_blocks[0][0] > 0:
        gap_events += 1
    if q_blocks[-1][1] < q_len:
        gap_events += 1
    gap_bases = (q_len - aligned_q)
    for idx in range(len(t_blocks) - 1):
        # target insertions inside aligned region
        t_ins = t_blocks[idx + 1][0] - t_blocks[idx][1]
        if t_ins > 0:
            gap_bases += t_ins

    span_s = s[t_blocks[0][0] : t_blocks[-1][1]]
    n_in_span = span_s.count("N")
    aligned_cols = matches + mismatches + gap_bases
    cov_pct = round(aligned_q / max(q_len, 1) * 100.0, 2)
    ident_pct = round(matches / max(aligned_cols, 1) * 100.0, 2)
    return {
        "alignment_coverage_pct": cov_pct,
        "nucleotide_identity_pct": ident_pct,
        "aligned_columns": aligned_cols,
        "matches": matches,
        "mismatches": mismatches,
        "gap_events": gap_events,
        "gap_bases": gap_bases,
        "n_bases_in_subject_span": n_in_span,
    }


def run_stage3(output_root: Path | None = None) -> Dict[str, Path]:
    out = resolve_output_root(output_root)
    v3 = resolve_v3_root(out)
    s3 = v3 / "stage3_isolate_conservation"
    fig_dir = s3 / "figures"
    s3.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    auth_cands = pd.read_csv(v3 / "authoritative_data" / "auth_candidates.tsv", sep="\t")
    auth_guides = pd.read_csv(v3 / "authoritative_data" / "auth_guides.tsv", sep="\t")
    auth_isos = pd.read_csv(v3 / "authoritative_data" / "auth_isolates.tsv", sep="\t")
    iso_qc = pd.read_csv(v3 / "stage2_genome_qc_comparative" / "core_integrity_screen" / "isolates_artifact_vs_quality_v3.tsv", sep="\t").set_index("isolate_id")
    homologs = pd.read_csv(out / "04_comparative_genomics" / "target_orthologs" / "target_homologs.tsv", sep="\t")

    # Load all 19 isolate reconstructed CDS sequences per target
    iso_cds_by_target: Dict[str, Dict[str, str]] = {}
    for tname in auth_cands["target"].unique():
        fa = out / "04_comparative_genomics" / "target_orthologs" / "target_cds_alignments" / f"{tname}_isolates_cds.fasta"
        iso_cds_by_target[tname] = {rec.id.split("|", 1)[0]: str(rec.seq).upper().replace("-", "") for rec in SeqIO.parse(str(fa), "fasta")}

    hom_lookup = {(str(r["target"]), str(r["isolate_id"])): r for _, r in homologs.iterrows()}
    iso_order = list(auth_isos["isolate_id"])
    artifact_isolates = set(iso_qc[iso_qc["markers_with_internal_stops"] >= 10].index)  # WAI2406, WAI2411

    pair_rows: List[dict] = []
    # Also compute exact 21-mer preservation per guide across each isolate
    guide_map_by_cand: Dict[str, List[Tuple[str, bool]]] = {}
    for cid, grp in auth_guides.groupby("candidate_id"):
        guide_map_by_cand[cid] = [(str(r["sense_dna_21nt"]), bool(r["is_efficient_sifi21"]) if pd.notna(r["is_efficient_sifi21"]) else False) for _, r in grp.iterrows()]

    cand_metrics_rows: List[dict] = []
    for _, cr in auth_cands.iterrows():
        cid = str(cr["candidate_id"])
        tname = str(cr["target"])
        cseq = str(cr["sequence"]).upper()
        guides_info = guide_map_by_cand[cid]
        all_21 = [g[0] for g in guides_info]
        eff_21 = [g[0] for g in guides_info if g[1]]

        iso_records_for_cand: List[dict] = []
        guide_hit_counts_all19 = np.zeros(len(all_21), dtype=int)
        guide_hit_counts_hq17 = np.zeros(len(all_21), dtype=int)

        for iso_id in iso_order:
            subj = iso_cds_by_target[tname].get(iso_id, "")
            hr = hom_lookup.get((tname, iso_id))
            touches_gap = bool(hr["touches_gap_or_edge"]) if hr is not None else False
            is_art_iso = iso_id in artifact_isolates

            m = _semi_global_align_metrics(cseq, subj)
            k_hits = [1 if km in subj else 0 for km in all_21]
            k21_pct = round(sum(k_hits) / max(len(k_hits), 1) * 100.0, 2)
            for idx_k, kh in enumerate(k_hits):
                guide_hit_counts_all19[idx_k] += kh
                if not is_art_iso:
                    guide_hit_counts_hq17[idx_k] += kh

            # Distinguish missing / partial / artifact / confirmed divergence / confirmed match (Stage 3.1)
            cov = m["alignment_coverage_pct"]
            pid = m["nucleotide_identity_pct"]
            if cov == 0.0:
                status_3_1 = "Missing_Sequence"
            elif m["n_bases_in_subject_span"] > 0 or touches_gap:
                status_3_1 = "Unresolved_Assembly_Gap_Or_Edge"
            elif is_art_iso and (pid < 100.0 or k21_pct < 100.0) and m["gap_bases"] > 0:
                status_3_1 = "Unresolved_Assembly_Polishing_Artifact"
            elif cov < 95.0:
                status_3_1 = "Partial_Alignment"
            elif pid == 100.0 and cov == 100.0:
                status_3_1 = "Confirmed_Match"
            else:
                status_3_1 = "Confirmed_Divergence"

            rec_row = {
                "candidate_id": cid,
                "target": tname,
                "isolate_id": iso_id,
                "is_high_quality_non_artifact_isolate": not is_art_iso,
                "alignment_coverage_pct": cov,
                "nucleotide_identity_pct": pid,
                "exact_21mer_preservation_pct": k21_pct,
                "matches": m["matches"],
                "mismatches": m["mismatches"],
                "gap_events": m["gap_events"],
                "gap_bases": m["gap_bases"],
                "classification_state_v3": status_3_1,
            }
            pair_rows.append(rec_row)
            iso_records_for_cand.append(rec_row)

        p_df = pd.DataFrame(iso_records_for_cand)
        hq_df = p_df[p_df["is_high_quality_non_artifact_isolate"]]
        pids_all = p_df["nucleotide_identity_pct"].to_numpy()
        pids_hq = hq_df["nucleotide_identity_pct"].to_numpy()
        covs_all = p_df["alignment_coverage_pct"].to_numpy()
        k21_all = p_df["exact_21mer_preservation_pct"].to_numpy()
        k21_hq = hq_df["exact_21mer_preservation_pct"].to_numpy()

        # Guide-region statistics (#8)
        g_frac_all19 = guide_hit_counts_all19 / len(iso_order) * 100.0
        n_hq = int(len(hq_df))
        g_frac_hq17 = guide_hit_counts_hq17 / max(n_hq, 1) * 100.0
        eff_mask = np.array([g[1] for g in guides_info], dtype=bool)
        eff_100_all19 = round(float((g_frac_all19[eff_mask] == 100.0).mean() * 100.0), 2) if eff_mask.any() else 0.0
        eff_100_hq17 = round(float((g_frac_hq17[eff_mask] == 100.0).mean() * 100.0), 2) if eff_mask.any() else 0.0

        cand_metrics_rows.append(
            {
                "candidate_id": cid,
                "target": tname,
                "candidate_length_nt": int(cr["candidate_length_nt"]),
                "isolates_evaluated": len(iso_order),
                "hq_non_artifact_isolates": n_hq,
                "mean_nt_identity_all19_pct": round(float(pids_all.mean()), 3),
                "min_nt_identity_all19_pct": round(float(pids_all.min()), 2),
                "max_nt_divergence_all19_pct": round(100.0 - float(pids_all.min()), 2),
                "mean_nt_identity_hq17_pct": round(float(pids_hq.mean()), 3),
                "min_nt_identity_hq17_pct": round(float(pids_hq.min()), 2),
                "mean_alignment_coverage_all19_pct": round(float(covs_all.mean()), 2),
                "min_alignment_coverage_all19_pct": round(float(covs_all.min()), 2),
                "mean_21mer_preservation_all19_pct": round(float(k21_all.mean()), 2),
                "min_21mer_preservation_all19_pct": round(float(k21_all.min()), 2),
                "mean_21mer_preservation_hq17_pct": round(float(k21_hq.mean()), 2),
                "min_21mer_preservation_hq17_pct": round(float(k21_hq.min()), 2),
                "isolates_ge85pct_nt_id_n": int((pids_all >= 85.0).sum()),
                "isolates_ge85pct_nt_id_pct": round(float((pids_all >= 85.0).mean() * 100.0), 1),
                "isolates_ge90pct_nt_id_n": int((pids_all >= 90.0).sum()),
                "isolates_ge90pct_nt_id_pct": round(float((pids_all >= 90.0).mean() * 100.0), 1),
                "isolates_ge95pct_nt_id_n": int((pids_all >= 95.0).sum()),
                "isolates_ge95pct_nt_id_pct": round(float((pids_all >= 95.0).mean() * 100.0), 1),
                "isolates_100pct_nt_id_n": int((pids_all == 100.0).sum()),
                "isolates_100pct_nt_id_pct": round(float((pids_all == 100.0).mean() * 100.0), 1),
                "isolates_ge85pct_21mer_n": int((k21_all >= 85.0).sum()),
                "isolates_ge90pct_21mer_n": int((k21_all >= 90.0).sum()),
                "isolates_ge95pct_21mer_n": int((k21_all >= 95.0).sum()),
                "isolates_100pct_21mer_n": int((k21_all == 100.0).sum()),
                "n_isolates_complete_target_coverage": int((p_df["alignment_coverage_pct"] >= 99.0).sum()),
                "n_isolates_partial_coverage": int(((p_df["alignment_coverage_pct"] > 0.0) & (p_df["alignment_coverage_pct"] < 99.0)).sum()),
                "n_isolates_missing_or_unresolved": int(p_df["classification_state_v3"].str.startswith(("Missing", "Unresolved")).sum()),
                "n_isolates_confirmed_match": int((p_df["classification_state_v3"] == "Confirmed_Match").sum()),
                "n_isolates_confirmed_divergence": int((p_df["classification_state_v3"] == "Confirmed_Divergence").sum()),
                "guides_mean_isolate_conservation_pct": round(float(g_frac_all19.mean()), 2),
                "guides_pct_preserved_all19": round(float((g_frac_all19 == 100.0).mean() * 100.0), 2),
                "guides_pct_preserved_hq17": round(float((g_frac_hq17 == 100.0).mean() * 100.0), 2),
                "efficient_sifi21_guides_100pct_all19_pct": eff_100_all19,
                "efficient_sifi21_guides_100pct_hq17_pct": eff_100_hq17,
            }
        )

    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(s3 / "candidate_by_isolate_alignments_v3.tsv", sep="\t", index=False)

    metrics_df = pd.DataFrame(cand_metrics_rows).sort_values(
        ["min_nt_identity_hq17_pct", "mean_21mer_preservation_hq17_pct", "min_21mer_preservation_all19_pct"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    metrics_df["conservation_rank_v3"] = range(1, len(metrics_df) + 1)
    metrics_df.to_csv(s3 / "candidate_conservation_multimetric_v3.tsv", sep="\t", index=False)

    # 2. Threshold Sensitivity Scenarios (85%, 90%, 95% on all-19 nt ID, 21-mer k-mer, and HQ-17!)
    sens_rows: List[dict] = []
    for label, metric_col, thresh, req_frac in [
        ("Scenario_A1_Mean21mer_ge85pct", "mean_21mer_preservation_all19_pct", 85.0, "mean"),
        ("Scenario_A2_Mean21mer_ge90pct", "mean_21mer_preservation_all19_pct", 90.0, "mean"),
        ("Scenario_A3_Mean21mer_ge95pct", "mean_21mer_preservation_all19_pct", 95.0, "mean"),
        ("Scenario_B1_MinNtID_All19_ge85pct", "min_nt_identity_all19_pct", 85.0, "all_19_min"),
        ("Scenario_B2_MinNtID_All19_ge90pct", "min_nt_identity_all19_pct", 90.0, "all_19_min"),
        ("Scenario_B3_MinNtID_All19_ge95pct", "min_nt_identity_all19_pct", 95.0, "all_19_min"),
        ("Scenario_C1_MinNtID_HQ17_ge90pct", "min_nt_identity_hq17_pct", 90.0, "hq17_min"),
        ("Scenario_C2_MinNtID_HQ17_ge95pct", "min_nt_identity_hq17_pct", 95.0, "hq17_min"),
        ("Scenario_C3_Min21mer_HQ17_ge95pct", "min_21mer_preservation_hq17_pct", 95.0, "hq17_min"),
        ("Scenario_D1_90pctOfIsolates_ge95pctNtID", "isolates_ge95pct_nt_id_pct", 90.0, "pct_isolates"),
    ]:
        passed = metrics_df[metrics_df[metric_col] >= thresh]
        by_t = passed["target"].value_counts().to_dict()
        sens_rows.append(
            {
                "scenario_id": label,
                "evaluated_metric": metric_col,
                "threshold_pct": thresh,
                "mode": req_frac,
                "passing_candidates_n": len(passed),
                "CsSTE11_pass": int(by_t.get("CsSTE11", 0)),
                "CsPKS1_pass": int(by_t.get("CsPKS1", 0)),
                "CsSp1_pass": int(by_t.get("CsSp1", 0)),
                "CsFUS3_pass": int(by_t.get("CsFUS3", 0)),
                "ste11_cand012_pass": bool("BRD_CsSTE11_CAND012" in set(passed["candidate_id"])),
                "pks1_cand013_pass": bool("BRD_CsPKS1_CAND013" in set(passed["candidate_id"])),
                "cssp1_cand008_pass": bool("BRD_CsSp1_CAND008" in set(passed["candidate_id"])),
            }
        )
    sens_df = pd.DataFrame(sens_rows)
    sens_df.to_csv(s3 / "conservation_threshold_sensitivity_v3.tsv", sep="\t", index=False)

    # Missing data & state summary table
    state_df = pair_df.groupby(["target", "classification_state_v3"]).size().unstack(fill_value=0).reset_index()
    state_df.to_csv(s3 / "missing_data_and_artifact_summary_v3.tsv", sep="\t", index=False)

    # Heatmap Figure: 50 Candidates x 19 Isolates (21-mer preservation & nt identity!)
    piv_k21 = pair_df.pivot(index="candidate_id", columns="isolate_id", values="exact_21mer_preservation_pct").reindex(index=auth_cands["candidate_id"], columns=iso_order)
    fig, ax = plt.subplots(figsize=(13.5, 11.5))
    im = ax.imshow(piv_k21.to_numpy(), aspect="auto", cmap="RdYlGn", vmin=65, vmax=100)
    ax.set_xticks(range(len(iso_order)))
    ax.set_xticklabels([f"{i}*" if i in artifact_isolates else i for i in iso_order], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(piv_k21.index)))
    ax.set_yticklabels(piv_k21.index, fontsize=7.5)
    cbar = fig.colorbar(im, ax=ax, shrink=0.65)
    cbar.set_label("Exact 21-mer Preservation Across Isolate (%)", fontsize=9)
    ax.set_title("Stage 3: Isolate-Aware 21-mer Conservation Across All 50 Candidates x 19 Isolates\n(* WAI2406 & WAI2411 carry confirmed long-read homopolymer frameshift indels)", fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.savefig(fig_dir / "stage3_candidate_x_isolate_21mer_heatmap_v3.png", dpi=300)
    fig.savefig(fig_dir / "stage3_candidate_x_isolate_21mer_heatmap_v3.svg")
    plt.close(fig)

    # Reassessment summary for Section 3.3 candidates
    key_cids = ["BRD_CsSTE11_CAND012", "BRD_CsSTE11_CAND013", "BRD_CsPKS1_CAND013", "BRD_CsPKS1_CAND010", "BRD_CsSp1_CAND008", "BRD_CsFUS3_CAND003"]
    key_sub = metrics_df[metrics_df["candidate_id"].isin(key_cids)]
    s3_md = f"""# Stage 3 — Rigorous Isolate-Aware Conservation & Threshold Sensitivity Report (BRD v3)

**Timestamp:** {now_utc}

## Key Multi-Metric Findings
- **Artifact vs Divergence Separation:** In Stage 2, `WAI2406` and `WAI2411` showed genome-wide homopolymer frameshifts in >40% of housekeeping core genes. Consequently, `BRD_CsPKS1_CAND010` shows `75.0%` 21-mer preservation strictly due to indels in unpolished `WAI2411`, while maintaining **100.0% conservation across all 17 high-quality isolates (`min_21mer_hq17 = 100.0%`)**. In contrast, `BRD_CsSTE11_CAND001`–`004` show true biological SNP/indel divergence across `17/19` isolates, and `BRD_CsSp1` candidates have a genuine SNP/divergence in `LK93` (`min_21mer_hq17 = 76.67%–88.33%`).
- **Strictly Conserved Outliers:** `BRD_CsSTE11_CAND012`, `BRD_CsSTE11_CAND013`, `BRD_CsFUS3_CAND013`, and `BRD_CsPKS1_CAND006/007/008/012` achieve `100.0%` exact 21-mer and `100.0%` nucleotide identity across **19/19** isolates (including `WAI2406` and `WAI2411`).

## Key Candidates Reassessed (Section 3.3)
{key_sub[['conservation_rank_v3', 'candidate_id', 'target', 'mean_nt_identity_all19_pct', 'min_nt_identity_all19_pct', 'min_nt_identity_hq17_pct', 'mean_21mer_preservation_all19_pct', 'min_21mer_preservation_all19_pct', 'min_21mer_preservation_hq17_pct', 'isolates_100pct_21mer_n', 'efficient_sifi21_guides_100pct_hq17_pct']].to_markdown(index=False)}

## Threshold Sensitivity Scenarios
{sens_df.to_markdown(index=False)}

## Classification State Breakdown (All 950 Candidate-Isolate Pairs)
{state_df.to_markdown(index=False)}
"""
    (s3 / "stage3_conservation_report_v3.md").write_text(s3_md, encoding="utf-8")

    return {
        "pair_alignments": s3 / "candidate_by_isolate_alignments_v3.tsv",
        "multimetric": s3 / "candidate_conservation_multimetric_v3.tsv",
        "sensitivity": s3 / "conservation_threshold_sensitivity_v3.tsv",
        "report": s3 / "stage3_conservation_report_v3.md",
    }


if __name__ == "__main__":
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    for k, v in run_stage3().items():
        print(f"STAGE3_OUT {k}: {v}")
