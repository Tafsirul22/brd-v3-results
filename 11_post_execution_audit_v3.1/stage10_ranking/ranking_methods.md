# Stage 10 — Multi-Dimensional Ranking Transparency & Baseline Continuity (BRD v3.1)

**Audited:** 2026-11-11 08:30:11 UTC

## Continuity with v3 Baseline Shortlist
No new validated biological or sequencing evidence altered the target genes or candidate coordinates relative to v3 (`11_post_execution_audit_v3`), and all `8` clean independent-region representatives remain unchanged:
- **Ranks 1–3 (`CsSTE11`):** `BRD_CsSTE11_CAND013` (`0.9418`, Reg01, `3` 20/21 hits), `BRD_CsSTE11_CAND014` (`0.9367`, Reg03, `4` 20/21 hits), `BRD_CsSTE11_CAND010` (`0.8861`, Reg02, `14` 20/21 hits).
- **Ranks 4–7 (`CsPKS1`, PHI:2260 unaffected pathogenicity caveat):** `BRD_CsPKS1_CAND012` (`0.7979`, Reg02, `7` 20/21 hits), `BRD_CsPKS1_CAND010` (`0.7762`, Reg04, `12` 20/21 hits), `BRD_CsPKS1_CAND013` (`0.7841`, Reg03, `6` 20/21 hits), `BRD_CsPKS1_CAND006` (`0.7836`, Reg01, `9` 20/21 hits).
- **Rank 8 (`CsSp1`, short 270-nt ORF & LK93 SNP caveat):** `BRD_CsSp1_CAND008` (`0.7363`, Reg01, `8` 20/21 hits).

## Explicit Ranking Equations & Divergence Explanation (`AUDIT_09`)
- **Heuristic composite prioritization index** (NOT a probability):  
  `0.30 * biological_target_relevance_score + 0.35 * dsrna_targetability_score + 0.35 * sequence_selectivity_score`
- **Why `independent_region_rank_v3_1` differs from `composite_index_rank` at Ranks 5–6:**  
  `independent_region_ranking_v3_1.tsv` orders regions lexicographically by `(passes_zero_exact_21mer_filter DESC, pareto_front_3d_regions ASC, heuristic_composite_prioritization_index DESC)`. Because **`BRD_CsPKS1_CAND010`** (`0.7762`, targetability `0.9009`) sits on **3D Region Pareto Tier 1**, it ranks `#5` ahead of Tier-2 **`BRD_CsPKS1_CAND013`** (`0.7841`, `#6`) under region-selection priority even though `CAND013` has a slightly higher linear weighted sum (`+0.0079`, `#5` by pure composite index).
- **Near-match burden prominence (`AUDIT_08`):** `7/8` clean shortlist leaders (`CAND010`, `CAND012`, `CAND006`, `CAND013`, `CAND008`, plus unranked representatives) retain `Elevated_20of21_Near_Match_Burden` (`3–14` unique `20/21` near matches) directly in all ranking tables. Zero exact 21-mers never equals zero off-target risk.
