"""Isolate-resampling bootstrap stability analysis for conservation and rankings."""

from __future__ import annotations

from typing import Dict, List
import numpy as np
import pandas as pd


def bootstrap_conservation_and_ranks(
    cand_df: pd.DataFrame,
    isolate_conservation_matrix: pd.DataFrame,
    weights: Dict[str, float],
    n_iterations: int = 1000,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Resample the B. sorokiniana isolate panel with replacement (n_iterations)
    to compute 95% confidence intervals on candidate 21-mer conservation
    and Spearman/top-tier rank stability.
    """
    rng = np.random.default_rng(seed)
    iso_cols = [c for c in isolate_conservation_matrix.columns if c != "candidate_id"]
    n_iso = len(iso_cols)
    if n_iso < 2 or cand_df.empty:
        return pd.DataFrame()

    aligned = cand_df.set_index("candidate_id").loc[isolate_conservation_matrix["candidate_id"]]
    # Map candidate_id -> base non-conservation weighted score component (in exact aligned row order!)
    base_non_cons = (
        aligned["off_target_safety_score"].astype(float) * weights["off_target_safety"]
        + aligned["target_evidence_score"].astype(float) * weights["target_evidence_and_expression"]
        + aligned["guide_quality_score"].astype(float) * weights["guide_quality_and_accessibility"]
        + aligned["tool_agreement_score"].astype(float) * weights["tool_agreement_and_seq_qc"]
    ).to_numpy()
    w_cons = weights["nucleotide_conservation"]
    mat = isolate_conservation_matrix[iso_cols].to_numpy(dtype=float)  # (n_cand, n_iso)

    n_cand = mat.shape[0]
    cons_samples = np.zeros((n_cand, n_iterations), dtype=float)
    rank_samples = np.zeros((n_cand, n_iterations), dtype=float)

    for b in range(n_iterations):
        idx = rng.integers(0, n_iso, size=n_iso)
        boot_cons = mat[:, idx].mean(axis=1)
        cons_samples[:, b] = boot_cons
        boot_total = base_non_cons + w_cons * boot_cons
        # Rank descending (1 = best)
        order = np.argsort(-boot_total, kind="mergesort")
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, n_cand + 1)
        rank_samples[:, b] = ranks

    rows = []
    for i, cid in enumerate(isolate_conservation_matrix["candidate_id"]):
        cs = cons_samples[i]
        rs = rank_samples[i]
        rows.append(
            {
                "candidate_id": cid,
                "target": aligned.loc[cid, "target"],
                "n_isolates_resampled": n_iso,
                "bootstrap_iterations": n_iterations,
                "observed_conservation": round(float(mat[i].mean()), 4),
                "boot_cons_mean": round(float(cs.mean()), 4),
                "boot_cons_ci95_low": round(float(np.percentile(cs, 2.5)), 4),
                "boot_cons_ci95_high": round(float(np.percentile(cs, 97.5)), 4),
                "boot_rank_mean": round(float(rs.mean()), 2),
                "boot_rank_sd": round(float(rs.std()), 2),
                "boot_rank_ci95_low": round(float(np.percentile(rs, 2.5)), 1),
                "boot_rank_ci95_high": round(float(np.percentile(rs, 97.5)), 1),
                "prob_top10_pct": round(float((rs <= 10).mean() * 100.0), 2),
                "prob_top20_pct": round(float((rs <= 20).mean() * 100.0), 2),
            }
        )
    return pd.DataFrame(rows).sort_values("boot_rank_mean")
