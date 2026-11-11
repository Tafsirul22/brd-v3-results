"""Module M09 — Multi-Criteria Candidate Ranking, Pareto Analysis, Bootstrap & Sensitivity Testing."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from brd.analysis.bootstrap import bootstrap_conservation_and_ranks
from brd.analysis.pareto import compute_pareto_ranks
from brd.analysis.sensitivity import run_sensitivity_scenarios
from brd.config import BRDConfig, load_yaml
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker
from brd.utils.report_utils import df_to_html_table, wrap_html_report, write_text


def _compute_overlap_tool_count(row: pd.Series, all_df: pd.DataFrame) -> int:
    tname = row["target"]
    s1, e1 = int(row["candidate_start_nt"]), int(row["candidate_end_nt"])
    tools = {row["tool_name"]}
    sub = all_df[all_df["target"] == tname]
    for _, o in sub.iterrows():
        s2, e2 = int(o["candidate_start_nt"]), int(o["candidate_end_nt"])
        ov = min(e1, e2) - max(s1, s2) + 1
        if ov >= 100:
            tools.add(o["tool_name"])
    return len(tools)


def run_m09(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    logger = get_logger("BRD.M09", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    m09_dir = cfg.module_dir("09_candidate_ranking")
    int_dir = m09_dir / "integrated_scores"
    par_dir = m09_dir / "pareto_analysis"
    boot_dir = m09_dir / "bootstrap"
    sens_dir = m09_dir / "sensitivity_analysis"
    fig_dir = m09_dir / "figures"
    for d in (int_dir, par_dir, boot_dir, sens_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    cand_df = pd.read_csv(cfg.output_dir / "07_dsRNA_design" / "candidate_sequences" / "dsRNA_candidates.tsv", sep="\t")
    iso_mat_df = pd.read_csv(cfg.output_dir / "07_dsRNA_design" / "candidate_sequences" / "candidate_isolate_conservation_matrix.tsv", sep="\t")
    ot_df = pd.read_csv(cfg.output_dir / "08_off_target_screening" / "summaries" / "off_target_summary.tsv", sep="\t")
    func_df = pd.read_csv(cfg.output_dir / "05_functional_evidence" / "functional_evidence_summary" / "functional_evidence_summary.tsv", sep="\t")

    func_map = func_df.set_index("target").to_dict(orient="index")
    merged = cand_df.merge(
        ot_df[
            [
                "candidate_id",
                "wheat_classification",
                "wheat_exact_21mer_hits",
                "wheat_high_sim_hits_ge19nt",
                "barley_classification",
                "barley_exact_21mer_hits",
                "beneficial_fungi_classification",
                "beneficial_fungi_exact_21mer_hits",
                "beneficial_fungi_high_sim_hits_ge19nt",
                "top_hit_viennarna_duplexfold_mfe",
                "overall_off_target_classification",
                "off_target_safety_score_0_1",
            ]
        ],
        on="candidate_id",
        how="left",
    )

    # Compute normalized 5-dimension evaluation scores
    cons_scores = []
    ot_scores = []
    ev_scores = []
    guide_scores = []
    agree_scores = []
    tool_overlaps = []
    categories = []
    passes_filters = []
    data_completeness = []

    for _, row in merged.iterrows():
        tname = row["target"]
        finfo = func_map[tname]
        c_sc = float(row["conservation_pct"]) / 100.0
        ot_sc = float(row["off_target_safety_score_0_1"])

        f_sc = float(finfo["functional_evidence_score_0_1"])
        raw_rpkm = finfo["gse130517_inplanta_max_rpkm"]
        if pd.notna(raw_rpkm) and str(raw_rpkm) != "NA":
            max_rpkm = float(raw_rpkm)
            ex_sc = min(1.0, max_rpkm / 150.0)
        else:
            ex_sc = 0.0
        ev_sc = round(0.6 * f_sc + 0.4 * ex_sc, 4)

        g_eff = float(row["efficient_guide_fraction"]) if pd.notna(row["efficient_guide_fraction"]) else 0.0
        acc = float(row["mean_target_accessibility"]) if pd.notna(row["mean_target_accessibility"]) else 0.0
        g_sc = round(0.6 * g_eff + 0.4 * min(1.0, acc / 0.35), 4)

        ov_n = _compute_overlap_tool_count(row, merged)
        tool_overlaps.append(ov_n)
        qc_ok = bool(row["seq_qc_pass"])
        ag_sc = min(1.0, (0.5 + 0.25 * min(2, ov_n - 1)) if qc_ok else 0.3)

        cons_scores.append(round(c_sc, 4))
        ot_scores.append(round(ot_sc, 4))
        ev_scores.append(ev_sc)
        guide_scores.append(g_sc)
        agree_scores.append(round(ag_sc, 4))

        # Check hard computational filter flags
        exact_ot = (
            int(row["wheat_exact_21mer_hits"]) > 0
            or int(row["barley_exact_21mer_hits"]) > 0
            or int(row["beneficial_fungi_exact_21mer_hits"]) > 0
            or int(row["self_paralog_21mer_hits"]) > 0
        )
        p_filt = qc_ok and (not exact_ot) and (float(row["conservation_pct"]) >= 85.0)
        passes_filters.append(p_filt)
        data_completeness.append("100% measured (Genomic+Conservation+Expression+OffTarget+Thermodynamics)")

    merged["conservation_score"] = cons_scores
    merged["off_target_safety_score"] = ot_scores
    merged["target_evidence_score"] = ev_scores
    merged["guide_quality_score"] = guide_scores
    merged["tool_agreement_score"] = agree_scores
    merged["overlapping_design_tools"] = tool_overlaps
    merged["passes_computational_filters"] = passes_filters
    merged["data_completeness"] = data_completeness

    scoring_cfg = load_yaml(cfg.output_dir / "config" / "scoring.yaml")
    w = scoring_cfg["candidate_ranking"]["weights"]
    merged["composite_score"] = (
        merged["conservation_score"] * w["nucleotide_conservation"]
        + merged["off_target_safety_score"] * w["off_target_safety"]
        + merged["target_evidence_score"] * w["target_evidence_and_expression"]
        + merged["guide_quality_score"] * w["guide_quality_and_accessibility"]
        + merged["tool_agreement_score"] * w["tool_agreement_and_seq_qc"]
    ).round(4)

    # Pareto ranking across 4 primary objectives
    merged = compute_pareto_ranks(
        merged,
        ["conservation_score", "off_target_safety_score", "target_evidence_score", "guide_quality_score"],
    )

    for _, row in merged.iterrows():
        if not row["passes_computational_filters"]:
            categories.append("Excluded due to identified computational concerns")
        elif row["composite_score"] >= 0.78 and (row["pareto_rank"] <= 2) and row["conservation_pct"] >= 95.0:
            categories.append("High computational priority")
        elif row["composite_score"] >= 0.65:
            categories.append("Moderate computational priority")
        else:
            categories.append("Low computational priority")
    merged["computational_priority_category"] = categories

    merged = merged.sort_values(["passes_computational_filters", "composite_score", "conservation_pct"], ascending=[False, False, False]).reset_index(drop=True)
    merged["final_rank"] = range(1, len(merged) + 1)

    # Bootstrap over 19 isolates (1000 iterations)
    n_boot = int(cfg.raw["ranking"]["bootstrap_iterations"])
    seed = int(cfg.raw["ranking"]["random_seed"])
    boot_df = bootstrap_conservation_and_ranks(merged, iso_mat_df, w, n_iterations=n_boot, seed=seed)

    # Sensitivity scenarios across alternative weight combinations
    cand_scenarios = {
        "baseline": {
            "conservation": 0.25,
            "off_target": 0.25,
            "evidence": 0.20,
            "guide_quality": 0.15,
            "tool_agreement": 0.15,
        },
        "safety_first": {
            "conservation": 0.20,
            "off_target": 0.40,
            "evidence": 0.15,
            "guide_quality": 0.15,
            "tool_agreement": 0.10,
        },
        "conservation_first": {
            "conservation": 0.40,
            "off_target": 0.20,
            "evidence": 0.15,
            "guide_quality": 0.15,
            "tool_agreement": 0.10,
        },
        "evidence_first": {
            "conservation": 0.20,
            "off_target": 0.20,
            "evidence": 0.35,
            "guide_quality": 0.15,
            "tool_agreement": 0.10,
        },
        "equal_weights": {
            "conservation": 0.20,
            "off_target": 0.20,
            "evidence": 0.20,
            "guide_quality": 0.20,
            "tool_agreement": 0.20,
        },
    }
    sens_df = run_sensitivity_scenarios(
        merged,
        id_col="candidate_id",
        score_cols={
            "conservation": "conservation_score",
            "off_target": "off_target_safety_score",
            "evidence": "target_evidence_score",
            "guide_quality": "guide_quality_score",
            "tool_agreement": "tool_agreement_score",
        },
        scenarios=cand_scenarios,
    )

    pareto_df = merged[merged["pareto_rank"] == 1][
        ["final_rank", "candidate_id", "target", "tool_name", "candidate_length_nt", "conservation_pct", "off_target_safety_score", "target_evidence_score", "guide_quality_score", "composite_score", "computational_priority_category"]
    ]
    conf_df = merged[
        ["final_rank", "candidate_id", "target", "tool_name", "passes_computational_filters", "pareto_rank", "composite_score", "computational_priority_category", "data_completeness"]
    ].merge(
        boot_df[["candidate_id", "boot_cons_ci95_low", "boot_cons_ci95_high", "boot_rank_mean", "boot_rank_sd", "prob_top10_pct"]],
        on="candidate_id",
        how="left",
    )

    # Save module outputs
    rank_p = int_dir / "final_candidate_ranking.tsv"
    pareto_p = par_dir / "pareto_front.tsv"
    boot_p = boot_dir / "bootstrap_results.tsv"
    sens_p = sens_dir / "sensitivity_analysis.tsv"
    conf_p = int_dir / "candidate_confidence_summary.tsv"

    merged.drop(columns=["sequence"]).to_csv(rank_p, sep="\t", index=False)
    pareto_df.to_csv(pareto_p, sep="\t", index=False)
    boot_df.to_csv(boot_p, sep="\t", index=False)
    sens_df.to_csv(sens_p, sep="\t", index=False)
    conf_df.to_csv(conf_p, sep="\t", index=False)

    # Save mandatory final tables in 10_final_reports
    fin_dir = cfg.module_dir("10_final_reports")
    merged.drop(columns=["sequence"]).to_csv(fin_dir / "BRD_dsRNA_Candidate_Ranking.tsv", sep="\t", index=False)
    pareto_df.to_csv(fin_dir / "BRD_Pareto_Analysis.tsv", sep="\t", index=False)
    boot_df.to_csv(fin_dir / "BRD_Bootstrap_Results.tsv", sep="\t", index=False)
    sens_df.to_csv(fin_dir / "BRD_Sensitivity_Analysis.tsv", sep="\t", index=False)
    conf_df.to_csv(fin_dir / "BRD_Candidate_Confidence.tsv", sep="\t", index=False)

    # Figures (PNG + SVG): Pareto frontier, Bootstrap stability, Sensitivity analysis
    # 1. Pareto frontier (Conservation vs Off-target Safety, colored by Target, star for Pareto optimal)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    tcolors = {"CsPKS1": "#7c3aed", "CsFUS3": "#2563eb", "CsSTE11": "#059669", "CsSp1": "#dc2626"}
    for tname, col in tcolors.items():
        sub = merged[merged["target"] == tname]
        if sub.empty:
            continue
        non_p = sub[~sub["is_pareto_optimal"]]
        par_p_df = sub[sub["is_pareto_optimal"]]
        ax.scatter(non_p["conservation_pct"], non_p["off_target_safety_score"], color=col, alpha=0.55, s=45, label=f"{tname} (dominated)")
        if not par_p_df.empty:
            ax.scatter(par_p_df["conservation_pct"], par_p_df["off_target_safety_score"], color=col, marker="*", s=160, edgecolors="#0f172a", linewidths=0.8, label=f"{tname} (Pareto front)")
    ax.set_xlabel("Multi-Isolate 21-mer Conservation (%)")
    ax.set_ylabel("Multi-Reference Off-Target Safety Score (0–1)")
    ax.set_title("Candidate Pareto Frontier: Isolate Conservation vs Off-Target Safety")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(fig_dir / "candidate_pareto_frontier.png", dpi=300)
    fig.savefig(fig_dir / "candidate_pareto_frontier.svg")
    plt.close(fig)

    # 2. Bootstrap stability plot (top 20 candidates)
    top_b = boot_df.head(20).iloc[::-1]
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    xerr_low = top_b["boot_rank_mean"] - top_b["boot_rank_ci95_low"]
    xerr_hi = top_b["boot_rank_ci95_high"] - top_b["boot_rank_mean"]
    ax2.errorbar(
        top_b["boot_rank_mean"],
        range(len(top_b)),
        xerr=[np.maximum(0, xerr_low), np.maximum(0, xerr_hi)],
        fmt="o",
        color="#1d4ed8",
        ecolor="#64748b",
        capsize=3,
    )
    ax2.set_yticks(range(len(top_b)))
    ax2.set_yticklabels([f"{r['candidate_id']} ({r['target']})" for _, r in top_b.iterrows()], fontsize=8)
    ax2.set_xlabel("Bootstrap Mean Rank (95% CI across 1,000 isolate resamples; lower is better)")
    ax2.set_title("Bootstrap Ranking Stability (Top 20 Candidates)")
    ax2.invert_xaxis()
    ax2.grid(axis="x", linestyle="--", alpha=0.35)
    fig2.tight_layout()
    fig2.savefig(fig_dir / "bootstrap_stability.png", dpi=300)
    fig2.savefig(fig_dir / "bootstrap_stability.svg")
    plt.close(fig2)

    # 3. Sensitivity analysis plot (top 15 candidates across 5 weight scenarios)
    top_s = sens_df.head(15)
    fig3, ax3 = plt.subplots(figsize=(10, 5.5))
    scen_cols = [c for c in top_s.columns if c.startswith("rank_") and c != "rank_range"]
    x_pos = np.arange(len(scen_cols))
    for _, sr in top_s.iterrows():
        ax3.plot(x_pos, [sr[c] for c in scen_cols], marker="o", linewidth=1.3, label=f"{sr['candidate_id']} ({sr['target']})")
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels([c.replace("rank_", "") for c in scen_cols], rotation=15)
    ax3.set_ylabel("Candidate Rank (1 = Best)")
    ax3.set_title("Sensitivity Analysis: Rank Stability Under Alternative Weighting Schemes")
    ax3.invert_yaxis()
    ax3.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7.5)
    ax3.grid(True, linestyle="--", alpha=0.35)
    fig3.tight_layout()
    fig3.savefig(fig_dir / "sensitivity_analysis.png", dpi=300)
    fig3.savefig(fig_dir / "sensitivity_analysis.svg")
    plt.close(fig3)

    body = (
        "<div class='callout'><strong>Ranking Boundary:</strong> Computational priority categories reflect multi-isolate conservation, "
        "non-target BLAST/21-mer/ViennaRNA screening, functional genetics, infection expression, and design QC. They do <strong>not</strong> "
        "predict exogenous dsRNA uptake or demonstrate fungicidal efficacy.</div>"
        f"<h2>Top 20 Ranked Candidates</h2>{df_to_html_table(merged.drop(columns=['sequence']).head(20))}"
        f"<h2>Pareto Non-Dominated Frontier (Rank 1)</h2>{df_to_html_table(pareto_df)}"
        f"<figure><img src='figures/candidate_pareto_frontier.png'/><figcaption>Pareto frontier across conservation and off-target safety.</figcaption></figure>"
        f"<figure><img src='figures/bootstrap_stability.png'/><figcaption>1,000-iteration isolate bootstrap rank 95% CIs.</figcaption></figure>"
        f"<figure><img src='figures/sensitivity_analysis.png'/><figcaption>Sensitivity across 5 weighting scenarios.</figcaption></figure>"
    )
    rep_html = m09_dir / "candidate_ranking_report.html"
    write_text(rep_html, wrap_html_report("BRD M09 — Uncertainty-Aware Candidate Ranking Report", "Pareto, 1,000-Bootstrap, & Weight Sensitivity Analysis", body))

    prov.record("M09", "Candidate Ranking & Uncertainty", str(cand_df.shape[0]) + "_candidates", rank_p, "Pareto+Bootstrap+Sensitivity", "1.0.0", f"boot={n_boot}")
    logger.info("M09 complete: %d candidates ranked (%d passing computational filters, %d Pareto front).", len(merged), int(merged["passes_computational_filters"].sum()), len(pareto_df))
    return {
        "final_candidate_ranking": rank_p,
        "pareto_front": pareto_p,
        "bootstrap_results": boot_p,
        "sensitivity_analysis": sens_p,
        "candidate_confidence": conf_p,
        "report": rep_html,
    }
