"""Pareto non-dominated sorting for multi-objective candidate ranking."""

from __future__ import annotations

from typing import List
import numpy as np
import pandas as pd


def compute_pareto_ranks(df: pd.DataFrame, objectives: List[str], maximize: List[bool] | None = None) -> pd.DataFrame:
    """Compute Pareto front tiers (layer 1 = non-dominated front) across numeric objectives."""
    if df.empty:
        out = df.copy()
        out["pareto_rank"] = pd.Series(dtype=int)
        out["is_pareto_optimal"] = pd.Series(dtype=bool)
        return out

    if maximize is None:
        maximize = [True] * len(objectives)

    mat = df[objectives].astype(float).to_numpy().copy()
    for j, is_max in enumerate(maximize):
        if not is_max:
            mat[:, j] = -mat[:, j]

    n = mat.shape[0]
    remaining = set(range(n))
    ranks = np.zeros(n, dtype=int)
    current_rank = 1

    while remaining:
        front = []
        rem_list = sorted(remaining)
        for i in rem_list:
            dominated = False
            for j in rem_list:
                if i == j:
                    continue
                if np.all(mat[j] >= mat[i]) and np.any(mat[j] > mat[i]):
                    dominated = True
                    break
            if not dominated:
                front.append(i)
        for i in front:
            ranks[i] = current_rank
            remaining.remove(i)
        current_rank += 1

    out = df.copy()
    out["pareto_rank"] = ranks
    out["is_pareto_optimal"] = out["pareto_rank"] == 1
    return out
