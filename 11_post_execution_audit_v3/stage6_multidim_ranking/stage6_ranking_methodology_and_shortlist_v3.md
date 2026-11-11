# Stage 6 — Multi-Dimensional Candidate Ranking Methodology & Shortlist (BRD v3)

**Timestamp:** 2026-11-11 06:55:06 UTC

## 1. Strict Separation of Evaluation Dimensions
To prevent cross-dimension leakage and premature aggregation:
1. **Dimension 1 — Biological Target Relevance (`biological_target_relevance_score`)**: Gene-level only (`0.35*KO + 0.20*PHI_context + 0.10*DEG + 0.10*InterPro + 0.15*Orthology + 0.10*Traceability`). Notably, `CsPKS1` direct *B. sorokiniana* `PHI:2260` records `unaffected pathogenicity`, dropping `CsPKS1` below `CsFUS3`, `CsSTE11`, and `CsSp1` despite having 21 InterPro PKS domains. Validated `GSE130517` (`3B/4B/5B`) RNA-seq expression is reported separately and excluded from all composite equations.
2. **Dimension 2 — dsRNA Targetability (`dsrna_targetability_score`)**: Sequence/isolate-level only (`0.25*min_nt_HQ17 + 0.25*min_21mer_HQ17 + 0.15*guides_100pct_HQ17 + 0.15*eff_sifi21 + 0.10*guide_accessibility + 0.10*seq_complexity`). Tool provenance (`LUCID`, `dsRNAbuilder`, `dsRNAmax`) is **excluded** from targetability and never counted as biological evidence.
3. **Dimension 3 — Heuristic Sequence Selectivity (`sequence_selectivity_score`)**: Derived from the expanded Stage 4 screen (18/18 covered non-target/self entries, both strands). Never interpreted as ecological safety.

## 2. Biological Target Relevance Summary
| target   |   biological_target_relevance_score |   subscore_direct_ko_pathogenicity |   subscore_phi_context_aware |   subscore_sequence_traceability |   gse130517_valid_3b4b5b_mean_rpkm_separate | phi_base_context_interpretation                                                                                                                                     |
|:---------|------------------------------------:|-----------------------------------:|-----------------------------:|---------------------------------:|--------------------------------------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CsFUS3   |                              0.9921 |                               1    |                         1    |                              1   |                                       62.53 | Direct/close PHI loss_of_pathogenicity (PHI:142__PHI:2366__PHI:2398__PHI:125577); 10/10 plant-fungal homologs show virulence reduction                              |
| CsSTE11  |                              0.9921 |                               1    |                         1    |                              1   |                                       22.12 | Direct/close PHI loss_of_pathogenicity (PHI:4778, PHI:2397__PHI:125579); 9/10 plant-fungal homologs show virulence reduction; note M. robertsii ortholog unaffected |
| CsPKS1   |                              0.5796 |                               0.25 |                         0.25 |                              1   |                                       13.19 | Direct B. sorokiniana entry PHI:2260 records 'unaffected pathogenicity' (secondary DHN melanin phenotype; mixed toxin roles only in distant homologs)               |
| CsSp1    |                              0.55   |                               1    |                         0    |                              0.5 |                                       79.3  | Absent from PHI-base v4.19 (effect relies on Zhang et al. 2022 primary paper, not PHI curation)                                                                     |

## 3. Independent Region Shortlist (Zero Redundant Window Inflation)
Selecting the highest-ranking candidate within each Stage 5 independent target region yields **12 independent region representatives** (8 pass the zero-exact-21mer filter):

|   independent_region_rank_v3 | candidate_id        | target   | region_id     |   region_member_count |   biological_target_relevance_score |   dsrna_targetability_score |   sequence_selectivity_score |   heuristic_composite_prioritization_index |   pareto_front_3d_regions | passes_zero_exact_21mer_filter   | selectivity_classification_v3          |
|-----------------------------:|:--------------------|:---------|:--------------|----------------------:|------------------------------------:|----------------------------:|-----------------------------:|-------------------------------------------:|--------------------------:|:---------------------------------|:---------------------------------------|
|                            1 | BRD_CsSTE11_CAND013 | CsSTE11  | CsSTE11_REG01 |                     8 |                              0.9921 |                      0.8855 |                        0.955 |                                     0.9418 |                         1 | True                             | Low_20of21_Near_Match_Signal           |
|                            2 | BRD_CsSTE11_CAND014 | CsSTE11  | CsSTE11_REG03 |                     1 |                              0.9921 |                      0.8859 |                        0.94  |                                     0.9367 |                         1 | True                             | Low_20of21_Near_Match_Signal           |
|                            3 | BRD_CsSTE11_CAND010 | CsSTE11  | CsSTE11_REG02 |                     2 |                              0.9921 |                      0.8913 |                        0.79  |                                     0.8861 |                         1 | True                             | Elevated_20of21_Near_Match_Burden      |
|                            4 | BRD_CsPKS1_CAND012  | CsPKS1   | CsPKS1_REG02  |                     1 |                              0.5796 |                      0.8879 |                        0.895 |                                     0.7979 |                         1 | True                             | Elevated_20of21_Near_Match_Burden      |
|                            5 | BRD_CsPKS1_CAND010  | CsPKS1   | CsPKS1_REG04  |                     6 |                              0.5796 |                      0.9009 |                        0.82  |                                     0.7762 |                         1 | True                             | Elevated_20of21_Near_Match_Burden      |
|                            6 | BRD_CsPKS1_CAND013  | CsPKS1   | CsPKS1_REG03  |                     1 |                              0.5796 |                      0.8335 |                        0.91  |                                     0.7841 |                         2 | True                             | Elevated_20of21_Near_Match_Burden      |
|                            7 | BRD_CsPKS1_CAND006  | CsPKS1   | CsPKS1_REG01  |                     5 |                              0.5796 |                      0.8772 |                        0.865 |                                     0.7836 |                         2 | True                             | Elevated_20of21_Near_Match_Burden      |
|                            8 | BRD_CsSp1_CAND008   | CsSp1    | CsSp1_REG01   |                    10 |                              0.55   |                      0.7524 |                        0.88  |                                     0.7363 |                         3 | True                             | Elevated_20of21_Near_Match_Burden      |
|                            9 | BRD_CsFUS3_CAND001  | CsFUS3   | CsFUS3_REG02  |                     7 |                              0.9921 |                      0.8958 |                        0.1   |                                     0.6462 |                         1 | False                            | Exact_21mer_NonTarget_Or_Paralog_Match |
|                           10 | BRD_CsFUS3_CAND003  | CsFUS3   | CsFUS3_REG03  |                     1 |                              0.9921 |                      0.8836 |                        0.09  |                                     0.6384 |                         2 | False                            | Exact_21mer_NonTarget_Or_Paralog_Match |
|                           11 | BRD_CsFUS3_CAND009  | CsFUS3   | CsFUS3_REG01  |                     5 |                              0.9921 |                      0.858  |                        0.1   |                                     0.6329 |                         2 | False                            | Exact_21mer_NonTarget_Or_Paralog_Match |
|                           12 | BRD_CsSTE11_CAND004 | CsSTE11  | CsSTE11_REG04 |                     3 |                              0.9921 |                      0.7789 |                        0.1   |                                     0.6052 |                         3 | False                            | Exact_21mer_NonTarget_Or_Paralog_Match |

## 4. Weight Sensitivity (Spearman rho across 5 scenarios)
| scenario_a                        | scenario_b                        |   spearman_rho |   p_value |
|:----------------------------------|:----------------------------------|---------------:|----------:|
| S1_Balanced_30_35_35              | S2_SelectivityFirst_20_30_50      |         1      |  0        |
| S1_Balanced_30_35_35              | S3_BioPathogenicityFirst_50_25_25 |         1      |  0        |
| S1_Balanced_30_35_35              | S4_TargetabilityFirst_20_50_30    |         0.9286 |  0.000863 |
| S1_Balanced_30_35_35              | S5_EqualWeights_33_33_33          |         1      |  0        |
| S2_SelectivityFirst_20_30_50      | S3_BioPathogenicityFirst_50_25_25 |         1      |  0        |
| S2_SelectivityFirst_20_30_50      | S4_TargetabilityFirst_20_50_30    |         0.9286 |  0.000863 |
| S2_SelectivityFirst_20_30_50      | S5_EqualWeights_33_33_33          |         1      |  0        |
| S3_BioPathogenicityFirst_50_25_25 | S4_TargetabilityFirst_20_50_30    |         0.9286 |  0.000863 |
| S3_BioPathogenicityFirst_50_25_25 | S5_EqualWeights_33_33_33          |         1      |  0        |
| S4_TargetabilityFirst_20_50_30    | S5_EqualWeights_33_33_33          |         0.9286 |  0.000863 |

## 5. Mandatory Scientific Interpretation Guard
- `heuristic_composite_prioritization_index` (`0.30*Bio + 0.35*Targ + 0.35*Sel`) is a **heuristic prioritization index**, NOT a probability of RNAi gene silencing, fungal uptake, crop protection, or environmental safety.
- No experimental efficacy against *B. sorokiniana* is claimed for any candidate in this computational pipeline.
