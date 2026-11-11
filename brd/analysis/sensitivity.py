"""Weight-perturbation sensitivity analysis for target and candidate rankings."""

from __future__ import annotations

from typing import Dict
import pandas as pd
from scipy.stats import spearmanr


def run_sensitivity_scenarios(
    df: pd.DataFrame,
    id_col: str,
    score_cols: Dict[str, str],
    scenarios: Dict[str, Dict[str, float]],
    baseline_name: str = "baseline",
) -> pd.DataFrame:
    """Evaluate ranking stability across alternative weighting scenarios."""
    if df.empty:
        return pd.DataFrame()

    scenario_scores: Dict[str, pd.Series] = {}
    scenario_ranks: Dict[str, pd.Series] = {}

    for sname, wmap in scenarios.items():
        total = pd.Series(0.0, index=df.index)
        for crit, col in score_cols.items():
            total += df[col].astype(float) * float(wmap.get(crit, 0.0))
        scenario_scores[sname] = total
        scenario_ranks[sname] = total.rank(ascending=False, method="min")

    base_ranks = scenario_ranks.get(baseline_name, next(iter(scenario_ranks.values())))
    rows = []
    for i, row in df.iterrows():
        rec = {id_col: row[id_col]}
        if "target" in df.columns and id_col != "target":
            rec["target"] = row["target"]
        ranks_list = []
        for sname in scenarios:
            rec[f"score_{sname}"] = round(float(scenario_scores[sname].loc[i]), 4)
            rk = int(scenario_ranks[sname].loc[i])
            rec[f"rank_{sname}"] = rk
            ranks_list.append(rk)
        rec["min_rank_across_scenarios"] = min(ranks_list)
        rec["max_rank_across_scenarios"] = max(ranks_list)
        rec["rank_range"] = max(ranks_list) - min(ranks_list)
        rows.append(rec)

    res_df = pd.DataFrame(rows).sort_values(f"rank_{baseline_name}")
    # Attach pairwise Spearman correlations vs baseline as attributes/metadata if needed
    corrs = {}
    for sname, rseries in scenario_ranks.items():
        rho, pval = spearmanr(base_ranks, rseries)
        corrs[sname] = (round(float(rho), 4), float(pval))
    res_df.attrs["spearman_vs_baseline"] = corrs
    return res_df
