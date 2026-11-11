"""Target scoring helper with explicit unknown/missing-data transparency."""

from __future__ import annotations

from typing import Dict
import pandas as pd


def compute_weighted_target_scores(df: pd.DataFrame, weights: Dict[str, float]) -> pd.DataFrame:
    out = df.copy()
    score_col_map = {
        "multi_isolate_nucleotide_conservation": "score_conservation",
        "off_target_safety": "score_off_target",
        "direct_functional_pathogenicity_evidence": "score_functional",
        "infection_stage_expression": "score_expression",
        "dsrna_designability": "score_designability",
    }
    total = pd.Series(0.0, index=out.index)
    known_weight = pd.Series(0.0, index=out.index)
    for crit, col in score_col_map.items():
        w = float(weights.get(crit, 0.0))
        vals = pd.to_numeric(out[col], errors="coerce")
        is_known = vals.notna()
        known_weight += w * is_known.astype(float)
        total += w * vals.fillna(0.0)
    out["weighted_score_conservative"] = total.round(4)
    out["known_weight_fraction"] = known_weight.round(4)
    out["normalized_known_score"] = (total / known_weight.replace(0, float("nan"))).round(4)
    return out
