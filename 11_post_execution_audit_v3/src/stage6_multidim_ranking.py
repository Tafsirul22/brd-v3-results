"""Stage 6: Transparent Multi-Dimensional Candidate Ranking, Pareto Analysis & Evidence Cards."""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from paths_v3 import resolve_output_root, resolve_v3_root


def _compute_pareto_tiers(df: pd.DataFrame, cols: List[str]) -> List[int]:
    """Non-dominated sorting (higher is better on all `cols`), returning 1-based Pareto front tier."""
    vals = df[cols].to_numpy(dtype=float)
    n = len(vals)
    remaining = set(range(n))
    tiers = [0] * n
    tier = 1
    while remaining:
        front: List[int] = []
        rem_list = sorted(remaining)
        for i in rem_list:
            dominated = False
            vi = vals[i]
            for j in rem_list:
                if i == j:
                    continue
                vj = vals[j]
                if np.all(vj >= vi - 1e-12) and np.any(vj > vi + 1e-12):
                    dominated = True
                    break
            if not dominated:
                front.append(i)
        for idx in front:
            tiers[idx] = tier
            remaining.remove(idx)
        tier += 1
    return tiers


def run_stage6(output_root: Path | None = None) -> Dict[str, Path]:
    out = resolve_output_root(output_root)
    v3 = resolve_v3_root(out)
    s6 = v3 / "stage6_multidim_ranking"
    fig_dir = s6 / "figures"
    s6.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # Load authoritative Stage 1, 2, 3, 4, 5 tables
    auth_t = pd.read_csv(v3 / "authoritative_data" / "auth_targets.tsv", sep="\t").set_index("target")
    auth_c = pd.read_csv(v3 / "authoritative_data" / "auth_candidates.tsv", sep="\t").set_index("candidate_id")
    cons_df = pd.read_csv(v3 / "stage3_isolate_conservation" / "candidate_conservation_multimetric_v3.tsv", sep="\t").set_index("candidate_id")
    pair_df = pd.read_csv(v3 / "stage3_isolate_conservation" / "candidate_by_isolate_alignments_v3.tsv", sep="\t")
    sel_df = pd.read_csv(v3 / "stage4_off_target_expanded" / "candidate_selectivity_expanded_v3.tsv", sep="\t").set_index("candidate_id")
    exact_df = pd.read_csv(v3 / "stage4_off_target_expanded" / "exact_21mer_matches_expanded_v3.tsv", sep="\t")
    near_df = pd.read_csv(v3 / "stage4_off_target_expanded" / "near_guide_matches_expanded_v3.tsv", sep="\t")
    long_df = pd.read_csv(v3 / "stage4_off_target_expanded" / "long_nucleotide_alignments_expanded_v3.tsv", sep="\t")
    clust_df = pd.read_csv(v3 / "stage5_redundancy_design" / "candidate_to_region_cluster_map_v3.tsv", sep="\t").set_index("candidate_id")
    char_df = pd.read_csv(v3 / "stage5_redundancy_design" / "candidate_sequence_characterization_v3.tsv", sep="\t").set_index("candidate_id")
    cov_gaps_df = pd.read_csv(v3 / "stage4_off_target_expanded" / "reference_panel_coverage_and_gaps_v3.tsv", sep="\t")

    # Verified M05 biological evidence files + validated v2 GSE130517 expression stats
    phi_df = pd.read_csv(out / "05_functional_evidence" / "phi_base" / "phi_base_evidence.tsv", sep="\t")
    deg_df = pd.read_csv(out / "05_functional_evidence" / "phi_base" / "targets_vs_deg_diamond.tsv", sep="\t", names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp"])
    ipr_df = pd.read_csv(out / "05_functional_evidence" / "interpro" / "interpro_annotations.tsv", sep="\t")
    bio_v2 = pd.read_csv(out / "11_post_execution_audit" / "BRD_Target_Biological_Prioritisation_v2.tsv", sep="\t").set_index("target")

    # 1. DIMENSION 1: Biological Target Relevance (strictly gene-level, context-aware PHI-base)
    bio_rows: List[dict] = []
    # Published direct perturbation in B. sorokiniana with documented evidence record:
    direct_exp_rubric = {
        "CsFUS3": (1.00, "Direct B. sorokiniana KO loses appressoria & barley/wheat pathogenicity (Leng & Zhong 2015)"),
        "CsSTE11": (1.00, "Direct B. sorokiniana KO abolishes pathogenicity (PHI:4778; Leng & Zhong 2015)"),
        "CsSp1": (1.00, "Direct B. sorokiniana KO reduces virulence on wheat leaves (Zhang et al. 2022, PMC8743017)"),
        "CsPKS1": (0.25, "Direct B. sorokiniana RNAi/KO causes melanin defect with UNAFFECTED pathogenicity (PHI:2260; Leng et al. 2011)"),
    }
    for tname, tr in auth_t.iterrows():
        ref_prot = str(tr["refseq_protein_id"])
        phi_sub = phi_df[phi_df["target"] == tname]
        direct_phi = phi_sub[phi_sub["evidence_class"].str.contains("Direct B. sorokiniana", na=False)]
        phenos_all = " | ".join(phi_sub["mutant_phenotype"].dropna().astype(str).tolist())
        # Context-aware PHI score: direct B. sorokiniana unaffected_pathogenicity overrides heterologous predictions!
        if direct_phi.empty and phi_sub.empty:
            phi_score = 0.0
            phi_context = "Absent from PHI-base v4.19 (effect relies on Zhang et al. 2022 primary paper, not PHI curation)"
        elif not direct_phi.empty and any("loss of pathogenicity" in str(p) for p in direct_phi["mutant_phenotype"]):
            phi_score = 1.0
            n_loss = int(phi_sub["mutant_phenotype"].str.contains("loss of pathogenicity|reduced virulence", na=False).sum())
            caveat_m = "; note M. robertsii ortholog unaffected" if "Metarhizium" in " ".join(phi_sub["phi_species"].astype(str)) else ""
            phi_context = f"Direct/close PHI loss_of_pathogenicity ({', '.join(direct_phi['phi_id'].astype(str))}); {n_loss}/{len(phi_sub)} plant-fungal homologs show virulence reduction{caveat_m}"
        elif not direct_phi.empty and all("unaffected pathogenicity" in str(p) and "loss" not in str(p) for p in direct_phi["mutant_phenotype"]):
            phi_score = 0.25
            phi_context = f"Direct B. sorokiniana entry {', '.join(direct_phi['phi_id'].astype(str))} records 'unaffected pathogenicity' (secondary DHN melanin phenotype; mixed toxin roles only in distant homologs)"
        else:
            phi_score = 0.5
            phi_context = f"Mixed heterologous PHI phenotypes ({phenos_all[:120]})"

        deg_hits = deg_df[deg_df["qseqid"] == ref_prot]
        deg_score = 1.0 if len(deg_hits) > 0 else 0.0

        ipr_sub = ipr_df[ipr_df["target"] == tname]
        n_ipr = int(ipr_sub["interpro_accession"].nunique())
        dom_score = round(min(1.0, n_ipr / 3.0), 4)

        bv2 = bio_v2.loc[tname]
        iso_pres_n = int(str(bv2["isolate_homolog_presence"]).split("/")[0])
        cds_comp_n = int(str(bv2["complete_cds_availability"]).split("/")[0])
        ortho_score = round(0.5 * (iso_pres_n / 19.0) + 0.5 * (cds_comp_n / 19.0), 4)

        ver_stat = str(tr["verification_status"])
        ver_score = 1.0 if ("Full RefSeq/GenBank Verified" in ver_stat) else 0.5

        exp_score, exp_note = direct_exp_rubric[tname]
        bio_rel = round(
            0.35 * exp_score
            + 0.20 * phi_score
            + 0.10 * deg_score
            + 0.10 * dom_score
            + 0.15 * ortho_score
            + 0.10 * ver_score,
            4,
        )
        bio_rows.append(
            {
                "target": tname,
                "refseq_protein_id": ref_prot,
                "locus_tag": str(tr["locus_tag"]),
                "cds_length_nt": int(tr["cds_length_nt"]),
                "subscore_direct_ko_pathogenicity": exp_score,
                "subscore_phi_context_aware": phi_score,
                "subscore_deg_heterologous_essential": deg_score,
                "subscore_interpro_domains": dom_score,
                "subscore_orthology_19iso": ortho_score,
                "subscore_sequence_traceability": ver_score,
                "biological_target_relevance_score": bio_rel,
                "direct_experimental_context": exp_note,
                "phi_base_context_interpretation": phi_context,
                "deg_homolog_hits_n": len(deg_hits),
                "interpro_domains_n": n_ipr,
                "gse130517_valid_3b4b5b_mean_rpkm_separate": float(bv2["gse130517_valid_3B4B5B_inplanta_mean_rpkm"]),
                "gse130517_valid_3b4b5b_raw_reads_separate": int(bv2["gse130517_valid_3B4B5B_raw_reads_sum"]),
                "gse130517_fold_change_vs_axenic_separate": float(bv2["gse130517_inplanta_fold_change"]),
                "expression_context_caveat": "Barley root RNA-seq (GSE130517 samples 3B/4B/5B) excluded from composite score; excludes 1-read 3A_I/2B_III artifacts; barley-root expression may differ from wheat foliar infection.",
            }
        )
    bio_df = pd.DataFrame(bio_rows).sort_values("biological_target_relevance_score", ascending=False).reset_index(drop=True)
    bio_df["biological_relevance_rank_v3"] = range(1, len(bio_df) + 1)
    bio_df.to_csv(s6 / "biological_target_relevance_v3.tsv", sep="\t", index=False)
    bio_map = bio_df.set_index("target").to_dict(orient="index")

    # 2. DIMENSION 2: dsRNA Targetability (strictly sequence/conservation/guide/accessibility/independence)
    def _calc_targetability(min_nt_hq: float, min_21_hq: float, g100_hq: float, eff_frac: float, g_acc: float, dinuc_ent: float, low_comp: float, gc_pct: float, max_hp: int) -> Tuple[float, float]:
        s_min_nt = max(0.0, min(1.0, (min_nt_hq - 75.0) / 25.0))
        s_min_21 = max(0.0, min(1.0, min_21_hq / 100.0))
        s_g100 = max(0.0, min(1.0, g100_hq / 100.0))
        s_eff = max(0.0, min(1.0, eff_frac))
        s_acc = max(0.0, min(1.0, g_acc / 0.5))
        gc_pen = 0.0 if 35.0 <= gc_pct <= 60.0 else 0.25
        hp_pen = 0.0 if max_hp <= 5 else 0.15
        s_comp = max(0.0, min(1.0, dinuc_ent * (1.0 - low_comp) - gc_pen - hp_pen))
        score = round(0.25 * s_min_nt + 0.25 * s_min_21 + 0.15 * s_g100 + 0.15 * s_eff + 0.10 * s_acc + 0.10 * s_comp, 4)
        return score, round(s_comp, 4)

    # Panel coverage fraction from Stage 4
    cov_ok = int((cov_gaps_df["transcript_screened"].astype(bool) | cov_gaps_df["genome_screened"].astype(bool)).sum())
    cov_both = int((cov_gaps_df["transcript_screened"].astype(bool) & cov_gaps_df["genome_screened"].astype(bool)).sum())
    cov_tot = len(cov_gaps_df)
    panel_cov_note = (
        f"{cov_ok}/{cov_tot} planned reference entries covered ({cov_both}/{cov_tot} with both transcriptome+genome; "
        "wheat/barley/maize screened via complete RefSeq transcriptomes only due to multi-GB whole-genome size; unsequenced soil/phyllosphere taxa remain uncovered)"
    )

    master_rows: List[dict] = []
    for cid in auth_c.index:
        ac = auth_c.loc[cid]
        tname = str(ac["target"])
        co = cons_df.loc[cid]
        se = sel_df.loc[cid]
        ch = char_df.loc[cid]
        cl = clust_df.loc[cid]
        bm = bio_map[tname]

        targ_score, comp_score = _calc_targetability(
            float(co["min_nt_identity_hq17_pct"]),
            float(co["min_21mer_preservation_hq17_pct"]),
            float(co["guides_pct_preserved_hq17"]),
            float(ch["efficient_sifi21_fraction"]),
            float(ch["guide_acc_mean"]),
            float(ch["dinuc_entropy_norm"]),
            float(ch["low_complexity_repeat_frac"]),
            float(ac["gc_percent"]),
            int(ac["max_homopolymer"]),
        )
        bio_score = float(bm["biological_target_relevance_score"])
        sel_score = float(se["heuristic_sequence_selectivity_index"])
        comp_heur = round(0.30 * bio_score + 0.35 * targ_score + 0.35 * sel_score, 4)

        master_rows.append(
            {
                "candidate_id": cid,
                "target": tname,
                "region_id": str(cl["region_id"]),
                "region_member_count": int(cl["region_member_count"]),
                "cds_coords": f"{int(ac['cds_start_nt'])}..{int(ac['cds_end_nt'])}",
                "nd90pr_genomic_span": f"{ac['genomic_scaffold']}:{int(ac['genomic_start_1based'])}..{int(ac['genomic_end_1based'])}({ac['genomic_strand']})",
                "candidate_length_nt": int(ac["candidate_length_nt"]),
                "cds_transcript_coverage_pct": float(ch["cds_transcript_coverage_pct"]),
                "gc_percent": float(ac["gc_percent"]),
                "max_homopolymer": int(ac["max_homopolymer"]),
                "tool_provenance_audit_only": str(ac["tool_name"]),
                "biological_target_relevance_score": bio_score,
                "dsrna_targetability_score": targ_score,
                "sequence_selectivity_score": sel_score,
                "subscore_seq_composition_complexity": comp_score,
                "mean_nt_identity_hq17_pct": float(co["mean_nt_identity_hq17_pct"]),
                "min_nt_identity_hq17_pct": float(co["min_nt_identity_hq17_pct"]),
                "min_nt_identity_all19_pct": float(co["min_nt_identity_all19_pct"]),
                "min_21mer_preservation_hq17_pct": float(co["min_21mer_preservation_hq17_pct"]),
                "guides_pct_preserved_hq17": float(co["guides_pct_preserved_hq17"]),
                "n_isolates_missing_or_unresolved": int(co["n_isolates_missing_or_unresolved"]),
                "unique_guides_21mer": int(ch["unique_guides_21mer"]),
                "efficient_sifi21_fraction": float(ch["efficient_sifi21_fraction"]),
                "guide_acc_mean": float(ch["guide_acc_mean"]),
                "candidate_self_fold_mfe_kcal": float(ch["candidate_self_fold_mfe_kcal"]),
                "passes_zero_exact_21mer_filter": bool(se["passes_zero_exact_21mer_filter"]),
                "selectivity_classification_v3": str(se["selectivity_classification_v3"]),
                "exact_21_crops_total": int(se["wheat_exact_21"]) + int(se["barley_exact_21"]) + int(se["rice_exact_21"]) + int(se["maize_exact_21"]),
                "exact_21_beneficial_fungi_tx": int(se["beneficial_fungi_tx_exact_21"]),
                "exact_21_beneficial_fungi_gen": int(se["beneficial_fungi_genome_exact_21"]),
                "exact_21_bacteria_and_self_paralog": int(se["bacteria_exact_21"]) + int(se["self_paralog_or_genome_exact_21"]),
                "near_20of21_guides_total": int(se["crops_unique_20of21_guides"]) + int(se["beneficial_fungi_unique_20of21_guides"]) + int(se["bacteria_unique_20of21_guides"]),
                "long_homology_ge30nt_hits": int(se["crops_long_homology_ge30nt"]) + int(se["beneficial_fungi_long_homology_ge30nt"]) + int(se["bacteria_long_homology_ge30nt"]),
                "heuristic_composite_prioritization_index": comp_heur,
            }
        )

    rank_df = pd.DataFrame(master_rows)
    dim_cols = ["biological_target_relevance_score", "dsrna_targetability_score", "sequence_selectivity_score"]
    rank_df["pareto_front_3d_all50"] = _compute_pareto_tiers(rank_df, dim_cols)

    # Identify best representative per independent region (Stage 5 clustering):
    # Prefer passing zero-exact-21mer filter first, then heuristic_composite_prioritization_index, then min_nt_identity_hq17_pct!
    rep_ids: Dict[str, str] = {}
    for rid, rgrp in rank_df.groupby("region_id"):
        best_r = rgrp.sort_values(
            ["passes_zero_exact_21mer_filter", "heuristic_composite_prioritization_index", "dsrna_targetability_score"],
            ascending=[False, False, False],
        ).iloc[0]
        rep_ids[rid] = str(best_r["candidate_id"])
    rank_df["is_independent_region_representative"] = rank_df.apply(lambda r: r["candidate_id"] == rep_ids[r["region_id"]], axis=1)

    rank_df = rank_df.sort_values(
        ["passes_zero_exact_21mer_filter", "pareto_front_3d_all50", "heuristic_composite_prioritization_index"],
        ascending=[False, True, False],
    ).reset_index(drop=True)
    rank_df["overall_heuristic_rank_all50"] = rank_df["heuristic_composite_prioritization_index"].rank(ascending=False, method="min").astype(int)

    # Independent Region Shortlist table (1 representative per independent region; also compute 3D Pareto among representatives!)
    reg_reps = rank_df[rank_df["is_independent_region_representative"]].copy().reset_index(drop=True)
    reg_reps["pareto_front_3d_regions"] = _compute_pareto_tiers(reg_reps, dim_cols)
    reg_reps = reg_reps.sort_values(
        ["passes_zero_exact_21mer_filter", "pareto_front_3d_regions", "heuristic_composite_prioritization_index"],
        ascending=[False, True, False],
    ).reset_index(drop=True)
    reg_reps["independent_region_rank_v3"] = range(1, len(reg_reps) + 1)

    rep_tier_map = reg_reps.set_index("candidate_id")["pareto_front_3d_regions"].to_dict()
    rep_rank_map = reg_reps.set_index("candidate_id")["independent_region_rank_v3"].to_dict()
    rank_df["pareto_front_3d_regions"] = rank_df["candidate_id"].map(rep_tier_map).fillna(0).astype(int)
    rank_df["independent_region_rank_v3"] = rank_df["candidate_id"].map(rep_rank_map).fillna(0).astype(int)

    rank_df.to_csv(s6 / "candidate_multidim_ranking_v3.tsv", sep="\t", index=False)
    reg_reps.to_csv(s6 / "independent_region_shortlist_v3.tsv", sep="\t", index=False)

    # 3. Weight Sensitivity Analysis across the Independent Regions (5 transparent weight vectors)
    scenarios = [
        ("S1_Balanced_30_35_35", (0.30, 0.35, 0.35)),
        ("S2_SelectivityFirst_20_30_50", (0.20, 0.30, 0.50)),
        ("S3_BioPathogenicityFirst_50_25_25", (0.50, 0.25, 0.25)),
        ("S4_TargetabilityFirst_20_50_30", (0.20, 0.50, 0.30)),
        ("S5_EqualWeights_33_33_33", (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)),
    ]
    sens_rows: List[dict] = []
    scen_ranks: Dict[str, np.ndarray] = {}
    clean_reps = reg_reps[reg_reps["passes_zero_exact_21mer_filter"]].copy().reset_index(drop=True)
    for sname, (wb, wt, ws) in scenarios:
        sc = wb * clean_reps["biological_target_relevance_score"] + wt * clean_reps["dsrna_targetability_score"] + ws * clean_reps["sequence_selectivity_score"]
        rnk = sc.rank(ascending=False, method="min").astype(int).to_numpy()
        scen_ranks[sname] = rnk
        for idx, cr in clean_reps.iterrows():
            sens_rows.append(
                {
                    "scenario": sname,
                    "w_bio": round(wb, 3),
                    "w_targetability": round(wt, 3),
                    "w_selectivity": round(ws, 3),
                    "candidate_id": cr["candidate_id"],
                    "target": cr["target"],
                    "region_id": cr["region_id"],
                    "scenario_score": round(float(sc.iloc[idx]), 4),
                    "scenario_rank_among_clean_regions": int(rnk[idx]),
                }
            )
    sens_df = pd.DataFrame(sens_rows)
    sens_df.to_csv(s6 / "stage6_weight_sensitivity_v3.tsv", sep="\t", index=False)

    corr_rows: List[dict] = []
    snames = list(scen_ranks.keys())
    for i in range(len(snames)):
        for j in range(i + 1, len(snames)):
            rho, pval = spearmanr(scen_ranks[snames[i]], scen_ranks[snames[j]])
            corr_rows.append({"scenario_a": snames[i], "scenario_b": snames[j], "spearman_rho": round(float(rho), 4), "p_value": float(f"{pval:.3e}")})
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(s6 / "stage6_weight_sensitivity_spearman_v3.tsv", sep="\t", index=False)

    # 4. 1,000-Iteration Isolate-Resampling Bootstrap on Clean Independent Region Representatives
    hq_iso_df = pair_df[~pair_df["isolate_id"].isin(["WAI2406", "WAI2411"])].copy()
    clean_cids = clean_reps["candidate_id"].tolist()
    iso_names = sorted(hq_iso_df["isolate_id"].unique().tolist())
    # Pre-index per-candidate arrays of (nt_id, k21_pres) indexed by isolate position 0..len(iso_names)-1
    iso_pos = {iso: k for k, iso in enumerate(iso_names)}
    nt_mat = np.zeros((len(clean_cids), len(iso_names)), dtype=float)
    k21_mat = np.zeros((len(clean_cids), len(iso_names)), dtype=float)
    cid_idx = {cid: r for r, cid in enumerate(clean_cids)}
    for _, pr in hq_iso_df[hq_iso_df["candidate_id"].isin(clean_cids)].iterrows():
        ci = cid_idx[str(pr["candidate_id"])]
        ii = iso_pos[str(pr["isolate_id"])]
        nt_mat[ci, ii] = float(pr["nucleotide_identity_pct"])
        k21_mat[ci, ii] = float(pr["exact_21mer_preservation_pct"])

    rng = np.random.default_rng(42)
    n_boot = 1000
    boot_ranks = np.zeros((len(clean_cids), n_boot), dtype=float)
    base_bio = clean_reps["biological_target_relevance_score"].to_numpy()
    base_sel = clean_reps["sequence_selectivity_score"].to_numpy()
    base_g100 = (clean_reps["guides_pct_preserved_hq17"].to_numpy() / 100.0)
    base_eff = clean_reps["efficient_sifi21_fraction"].to_numpy()
    base_acc = np.minimum(1.0, clean_reps["guide_acc_mean"].to_numpy() / 0.5)
    base_comp = clean_reps["subscore_seq_composition_complexity"].to_numpy()

    for b in range(n_boot):
        draw = rng.integers(0, len(iso_names), size=len(iso_names))
        min_nt_b = nt_mat[:, draw].min(axis=1)
        min_21_b = k21_mat[:, draw].min(axis=1)
        s_nt_b = np.clip((min_nt_b - 75.0) / 25.0, 0.0, 1.0)
        s_21_b = np.clip(min_21_b / 100.0, 0.0, 1.0)
        targ_b = 0.25 * s_nt_b + 0.25 * s_21_b + 0.15 * base_g100 + 0.15 * base_eff + 0.10 * base_acc + 0.10 * base_comp
        comp_b = 0.30 * base_bio + 0.35 * targ_b + 0.35 * base_sel
        order = np.argsort(-comp_b, kind="mergesort")
        rnk_b = np.empty_like(order, dtype=float)
        rnk_b[order] = np.arange(1, len(clean_cids) + 1)
        boot_ranks[:, b] = rnk_b

    boot_rows: List[dict] = []
    for i, cr in clean_reps.iterrows():
        br = boot_ranks[i]
        boot_rows.append(
            {
                "candidate_id": cr["candidate_id"],
                "target": cr["target"],
                "region_id": cr["region_id"],
                "observed_clean_region_rank": int(cr["independent_region_rank_v3"]),
                "bootstrap_mean_rank": round(float(br.mean()), 2),
                "bootstrap_median_rank": round(float(np.median(br)), 2),
                "bootstrap_rank_sd": round(float(br.std()), 2),
                "bootstrap_rank_ci95_low": round(float(np.percentile(br, 2.5)), 2),
                "bootstrap_rank_ci95_high": round(float(np.percentile(br, 97.5)), 2),
                "prob_top1_pct": round(float((br == 1).mean() * 100.0), 1),
                "prob_top3_pct": round(float((br <= 3).mean() * 100.0), 1),
                "prob_top5_pct": round(float((br <= 5).mean() * 100.0), 1),
            }
        )
    boot_df = pd.DataFrame(boot_rows).sort_values(["observed_clean_region_rank", "bootstrap_mean_rank"]).reset_index(drop=True)
    boot_df.to_csv(s6 / "stage6_bootstrap_stability_v3.tsv", sep="\t", index=False)

    # 5. Per-Candidate Evidence Cards (ALL 50 candidates in structured JSON + human-readable Markdown!)
    cards: List[dict] = []
    card_md_lines = [
        "# BRD v3 — Complete Candidate Evidence Cards (All 50 Candidates)",
        "",
        f"**Generated:** {now_utc}",
        "**Scope Guard:** All metrics below are computational predictions across 19 *B. sorokiniana* isolates and 17 Stage 4 reference organisms. They do **not** represent wet-lab RNAi silencing efficiency or ecological safety.",
        "",
    ]
    boot_by_cid = boot_df.set_index("candidate_id").to_dict(orient="index")
    for _, rr in rank_df.iterrows():
        cid = str(rr["candidate_id"])
        tname = str(rr["target"])
        ex_hits = exact_df[exact_df["candidate_id"] == cid]
        nr_hits = near_df[(near_df["candidate_id"] == cid) & (near_df["match_tier"] == "near_20_of_21")]
        lg_hits = long_df[(long_df["candidate_id"] == cid) & (long_df["is_high_homology_ge30_80pct"])]

        ex_orgs = sorted({f"{r['organism']} ({r['screen_scope']})" for _, r in ex_hits.iterrows()})
        nr_orgs = sorted({f"{r['organism']} ({r['screen_scope']})" for _, r in nr_hits.iterrows()})
        lg_orgs = sorted({f"{r['organism']} ({r['alignment_length_nt']}nt, {r['pident']}%)" for _, r in lg_hits.iterrows()})

        weaknesses: List[str] = []
        if not rr["passes_zero_exact_21mer_filter"]:
            weaknesses.append(f"EXCLUDED BY ZERO-21MER FILTER: Exact 21-mer matches detected in {', '.join(ex_orgs)}")
        if int(rr["near_20of21_guides_total"]) > 0:
            weaknesses.append(f"Carries {int(rr['near_20of21_guides_total'])} unique 20/21-mer near-matches ({', '.join(nr_orgs[:4])})")
        if int(rr["long_homology_ge30nt_hits"]) > 0:
            weaknesses.append(f"Carries {int(rr['long_homology_ge30nt_hits'])} >=30 nt (>=80% id) alignment(s) ({', '.join(lg_orgs[:3])})")
        if float(rr["min_nt_identity_hq17_pct"]) < 95.0:
            weaknesses.append(f"Lower minimum nucleotide conservation across HQ-17 isolates ({rr['min_nt_identity_hq17_pct']}%)")
        if int(rr["n_isolates_missing_or_unresolved"]) > 0:
            weaknesses.append(f"{int(rr['n_isolates_missing_or_unresolved'])} isolate(s) have incomplete/artifactual assemblies (WAI2406/WAI2411)")
        if tname == "CsPKS1":
            weaknesses.append("Target biology caveat: direct B. sorokiniana PHI:2260 Pks1 perturbation showed UNAFFECTED pathogenicity")
        if tname == "CsSp1":
            weaknesses.append("Target traceability caveat: lacks independent GenBank nuccore/PHI-base accession; very short 270-nt CDS")
        if float(rr["gc_percent"]) < 40.0 or float(rr["gc_percent"]) > 58.0:
            weaknesses.append(f"GC content ({rr['gc_percent']}%) approaches design margin")
        weaknesses.append(f"Panel coverage caveat: {panel_cov_note}")

        strengths: List[str] = []
        if float(rr["biological_target_relevance_score"]) >= 0.85:
            strengths.append(f"Strong target biology relevance ({rr['biological_target_relevance_score']}) backed by direct B. sorokiniana knockout loss-of-pathogenicity evidence")
        if float(rr["min_nt_identity_hq17_pct"]) >= 98.0:
            strengths.append(f"High cross-isolate uniformity (HQ-17 min nt identity = {rr['min_nt_identity_hq17_pct']}%, min 21-mer = {rr['min_21mer_preservation_hq17_pct']}%)")
        if rr["passes_zero_exact_21mer_filter"]:
            strengths.append(f"Zero exact 21-mer matches detected across all {cov_ok} covered non-target/self reference entries")
        if int(rr["pareto_front_3d_all50"]) == 1:
            strengths.append("On Tier-1 3D Pareto frontier (non-dominated across biology, targetability, selectivity)")
        if rr["is_independent_region_representative"]:
            strengths.append(f"Selected representative for independent target region {rr['region_id']} (cluster of {int(rr['region_member_count'])} window(s))")

        b_info = boot_by_cid.get(cid, {})
        card = {
            "candidate_id": cid,
            "target": tname,
            "region_id": rr["region_id"],
            "is_independent_region_representative": bool(rr["is_independent_region_representative"]),
            "coordinates": {"cds_1based": rr["cds_coords"], "nd90pr_genomic": rr["nd90pr_genomic_span"], "length_nt": int(rr["candidate_length_nt"])},
            "dimensions": {
                "biological_target_relevance": float(rr["biological_target_relevance_score"]),
                "dsrna_targetability": float(rr["dsrna_targetability_score"]),
                "sequence_selectivity": float(rr["sequence_selectivity_score"]),
                "heuristic_composite_prioritization_index": float(rr["heuristic_composite_prioritization_index"]),
            },
            "ranking": {
                "overall_heuristic_rank_all50": int(rr["overall_heuristic_rank_all50"]),
                "pareto_front_3d_all50": int(rr["pareto_front_3d_all50"]),
                "independent_region_rank_v3": int(rr["independent_region_rank_v3"]),
                "pareto_front_3d_regions": int(rr["pareto_front_3d_regions"]),
                "passes_zero_exact_21mer_filter": bool(rr["passes_zero_exact_21mer_filter"]),
                "selectivity_classification_v3": rr["selectivity_classification_v3"],
                "bootstrap_stability_if_clean_rep": b_info,
            },
            "strengths": strengths,
            "weaknesses_and_uncertainties": weaknesses,
            "missing_data_caveats": [bio_map[tname]["expression_context_caveat"], panel_cov_note],
        }
        cards.append(card)
        card_md_lines.extend(
            [
                f"## {cid} ({tname} | `{rr['region_id']}` | Rep={rr['is_independent_region_representative']})",
                f"- **Coordinates:** CDS `{rr['cds_coords']}` (`{rr['candidate_length_nt']} nt`, {rr['cds_transcript_coverage_pct']}% CDS) | Genomic `{rr['nd90pr_genomic_span']}` | Provenance: `{rr['tool_provenance_audit_only']}`",
                f"- **3D Scores:** Biological Relevance = `{rr['biological_target_relevance_score']}`, Targetability = `{rr['dsrna_targetability_score']}`, Selectivity = `{rr['sequence_selectivity_score']}` | Optional Heuristic Composite = `{rr['heuristic_composite_prioritization_index']}` (Rank `{rr['overall_heuristic_rank_all50']}/50`, Pareto Tier `{rr['pareto_front_3d_all50']}`)",
                f"- **Conservation:** HQ-17 min/mean nt = `{rr['min_nt_identity_hq17_pct']}% / {rr['mean_nt_identity_hq17_pct']}%`; min 21-mer = `{rr['min_21mer_preservation_hq17_pct']}%`; All-19 min nt = `{rr['min_nt_identity_all19_pct']}%`",
                f"- **Selectivity (`{rr['selectivity_classification_v3']}`):** Crops 21-mer=`{rr['exact_21_crops_total']}`, BenFungi(tx/gen)=`{rr['exact_21_beneficial_fungi_tx']}/{rr['exact_21_beneficial_fungi_gen']}`, Bac+Paralog=`{rr['exact_21_bacteria_and_self_paralog']}`, 20/21 guides=`{rr['near_20of21_guides_total']}`, Long>=30nt=`{rr['long_homology_ge30nt_hits']}`",
                f"- **Strengths:** {'; '.join(strengths) if strengths else 'None highlighted'}",
                f"- **Weaknesses & Uncertainty:** {'; '.join(weaknesses)}",
                "",
            ]
        )

    (s6 / "candidate_evidence_cards_v3.json").write_text(json.dumps(cards, indent=2), encoding="utf-8")
    (s6 / "candidate_evidence_cards_v3.md").write_text("\n".join(card_md_lines), encoding="utf-8")

    # 6. Figures (3D Pareto / Selectivity vs Targetability bubble chart + Weight/Bootstrap stability)
    fig, ax = plt.subplots(figsize=(10.5, 6.0))
    colors = {"CsSTE11": "#059669", "CsSp1": "#7c3aed", "CsPKS1": "#d97706", "CsFUS3": "#dc2626"}
    for tname, grp in rank_df.groupby("target"):
        p_pass = grp[grp["passes_zero_exact_21mer_filter"]]
        p_fail = grp[~grp["passes_zero_exact_21mer_filter"]]
        if not p_pass.empty:
            ax.scatter(
                p_pass["dsrna_targetability_score"], p_pass["sequence_selectivity_score"],
                s=p_pass["biological_target_relevance_score"] * 220, color=colors[tname],
                alpha=0.82, edgecolors="#0f172a", linewidth=1.2, label=f"{tname} (0 exact 21-mer hits)"
            )
        if not p_fail.empty:
            ax.scatter(
                p_fail["dsrna_targetability_score"], p_fail["sequence_selectivity_score"],
                s=p_fail["biological_target_relevance_score"] * 170, color=colors[tname],
                marker="X", alpha=0.45, edgecolors="#0f172a", linewidth=0.7, label=f"{tname} (Excluded >=1 exact 21-mer)"
            )
    for _, cr in clean_reps.head(6).iterrows():
        ax.annotate(
            f"{cr['candidate_id'].replace('BRD_', '')}\n({cr['region_id']})",
            (cr["dsrna_targetability_score"], cr["sequence_selectivity_score"]),
            xytext=(6, 4), textcoords="offset points", fontsize=7.5, fontweight="bold",
        )
    ax.set_xlabel("Dimension 2: dsRNA Targetability Score (HQ-17 conservation, guides, accessibility, composition)", fontweight="bold")
    ax.set_ylabel("Dimension 3: Heuristic Sequence Selectivity Index (Expanded 17-organism panel)", fontweight="bold")
    ax.set_title("Stage 6: Multi-Dimensional Tradeoff Across 50 Candidates (Bubble size = Biological Relevance)", fontweight="bold")
    ax.set_xlim(0.25, 1.02)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(loc="lower left", fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(fig_dir / "stage6_3d_pareto_tradeoff_v3.png", dpi=300)
    fig.savefig(fig_dir / "stage6_3d_pareto_tradeoff_v3.svg")
    plt.close(fig)

    # Bootstrap stability plot for clean independent region representatives
    fig, ax = plt.subplots(figsize=(10, 5.2))
    b_top = boot_df.iloc[::-1].reset_index(drop=True)
    y_idx = np.arange(len(b_top))
    ax.errorbar(
        b_top["bootstrap_mean_rank"], y_idx,
        xerr=[b_top["bootstrap_mean_rank"] - b_top["bootstrap_rank_ci95_low"], b_top["bootstrap_rank_ci95_high"] - b_top["bootstrap_mean_rank"]],
        fmt="o", color="#0284c7", ecolor="#64748b", capsize=3, markersize=6,
    )
    ax.set_yticks(y_idx)
    ax.set_yticklabels([f"{r['candidate_id']} ({r['region_id']})" for _, r in b_top.iterrows()], fontsize=8.5)
    ax.invert_xaxis()
    ax.set_xlabel("Bootstrap Rank Among Clean Independent Region Representatives (1 = Best, 95% CI across 1,000 HQ-isolate resamples)", fontweight="bold")
    ax.set_title("Stage 6: Isolate-Resampling Bootstrap Stability of Independent Target Regions", fontweight="bold")
    ax.grid(True, axis="x", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(fig_dir / "stage6_weight_bootstrap_stability_v3.png", dpi=300)
    fig.savefig(fig_dir / "stage6_weight_bootstrap_stability_v3.svg")
    plt.close(fig)

    # 7. Methodology & Shortlist Markdown Report
    meth_md = f"""# Stage 6 — Multi-Dimensional Candidate Ranking Methodology & Shortlist (BRD v3)

**Timestamp:** {now_utc}

## 1. Strict Separation of Evaluation Dimensions
To prevent cross-dimension leakage and premature aggregation:
1. **Dimension 1 — Biological Target Relevance (`biological_target_relevance_score`)**: Gene-level only (`0.35*KO + 0.20*PHI_context + 0.10*DEG + 0.10*InterPro + 0.15*Orthology + 0.10*Traceability`). Notably, `CsPKS1` direct *B. sorokiniana* `PHI:2260` records `unaffected pathogenicity`, dropping `CsPKS1` below `CsFUS3`, `CsSTE11`, and `CsSp1` despite having 21 InterPro PKS domains. Validated `GSE130517` (`3B/4B/5B`) RNA-seq expression is reported separately and excluded from all composite equations.
2. **Dimension 2 — dsRNA Targetability (`dsrna_targetability_score`)**: Sequence/isolate-level only (`0.25*min_nt_HQ17 + 0.25*min_21mer_HQ17 + 0.15*guides_100pct_HQ17 + 0.15*eff_sifi21 + 0.10*guide_accessibility + 0.10*seq_complexity`). Tool provenance (`LUCID`, `dsRNAbuilder`, `dsRNAmax`) is **excluded** from targetability and never counted as biological evidence.
3. **Dimension 3 — Heuristic Sequence Selectivity (`sequence_selectivity_score`)**: Derived from the expanded Stage 4 screen ({cov_ok}/{cov_tot} covered non-target/self entries, both strands). Never interpreted as ecological safety.

## 2. Biological Target Relevance Summary
{bio_df[['target', 'biological_target_relevance_score', 'subscore_direct_ko_pathogenicity', 'subscore_phi_context_aware', 'subscore_sequence_traceability', 'gse130517_valid_3b4b5b_mean_rpkm_separate', 'phi_base_context_interpretation']].to_markdown(index=False)}

## 3. Independent Region Shortlist (Zero Redundant Window Inflation)
Selecting the highest-ranking candidate within each Stage 5 independent target region yields **{len(reg_reps)} independent region representatives** ({len(clean_reps)} pass the zero-exact-21mer filter):

{reg_reps[['independent_region_rank_v3', 'candidate_id', 'target', 'region_id', 'region_member_count', 'biological_target_relevance_score', 'dsrna_targetability_score', 'sequence_selectivity_score', 'heuristic_composite_prioritization_index', 'pareto_front_3d_regions', 'passes_zero_exact_21mer_filter', 'selectivity_classification_v3']].to_markdown(index=False)}

## 4. Weight Sensitivity (Spearman rho across 5 scenarios)
{corr_df.to_markdown(index=False)}

## 5. Mandatory Scientific Interpretation Guard
- `heuristic_composite_prioritization_index` (`0.30*Bio + 0.35*Targ + 0.35*Sel`) is a **heuristic prioritization index**, NOT a probability of RNAi gene silencing, fungal uptake, crop protection, or environmental safety.
- No experimental efficacy against *B. sorokiniana* is claimed for any candidate in this computational pipeline.
"""
    (s6 / "stage6_ranking_methodology_and_shortlist_v3.md").write_text(meth_md, encoding="utf-8")

    return {
        "bio_table": s6 / "biological_target_relevance_v3.tsv",
        "ranking_all50": s6 / "candidate_multidim_ranking_v3.tsv",
        "region_shortlist": s6 / "independent_region_shortlist_v3.tsv",
        "weight_sensitivity": s6 / "stage6_weight_sensitivity_v3.tsv",
        "bootstrap": s6 / "stage6_bootstrap_stability_v3.tsv",
        "cards_json": s6 / "candidate_evidence_cards_v3.json",
        "cards_md": s6 / "candidate_evidence_cards_v3.md",
    }


if __name__ == "__main__":
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    for k, v in run_stage6().items():
        print(f"STAGE6_OUT {k}: {v}")
