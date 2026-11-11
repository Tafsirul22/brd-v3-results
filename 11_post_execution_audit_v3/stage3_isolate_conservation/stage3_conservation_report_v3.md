# Stage 3 — Rigorous Isolate-Aware Conservation & Threshold Sensitivity Report (BRD v3)

**Timestamp:** 2026-11-11 06:20:03 UTC

## Key Multi-Metric Findings
- **Artifact vs Divergence Separation:** In Stage 2, `WAI2406` and `WAI2411` showed genome-wide homopolymer frameshifts in >40% of housekeeping core genes. Consequently, `BRD_CsPKS1_CAND010` shows `75.0%` 21-mer preservation strictly due to indels in unpolished `WAI2411`, while maintaining **100.0% conservation across all 17 high-quality isolates (`min_21mer_hq17 = 100.0%`)**. In contrast, `BRD_CsSTE11_CAND001`–`004` show true biological SNP/indel divergence across `17/19` isolates, and `BRD_CsSp1` candidates have a genuine SNP/divergence in `LK93` (`min_21mer_hq17 = 76.67%–88.33%`).
- **Strictly Conserved Outliers:** `BRD_CsSTE11_CAND012`, `BRD_CsSTE11_CAND013`, `BRD_CsFUS3_CAND013`, and `BRD_CsPKS1_CAND006/007/008/012` achieve `100.0%` exact 21-mer and `100.0%` nucleotide identity across **19/19** isolates (including `WAI2406` and `WAI2411`).

## Key Candidates Reassessed (Section 3.3)
|   conservation_rank_v3 | candidate_id        | target   |   mean_nt_identity_all19_pct |   min_nt_identity_all19_pct |   min_nt_identity_hq17_pct |   mean_21mer_preservation_all19_pct |   min_21mer_preservation_all19_pct |   min_21mer_preservation_hq17_pct |   isolates_100pct_21mer_n |   efficient_sifi21_guides_100pct_hq17_pct |
|-----------------------:|:--------------------|:---------|-----------------------------:|----------------------------:|---------------------------:|------------------------------------:|-----------------------------------:|----------------------------------:|--------------------------:|------------------------------------------:|
|                      6 | BRD_CsSTE11_CAND012 | CsSTE11  |                      100     |                      100    |                     100    |                              100    |                             100    |                            100    |                        19 |                                    100    |
|                      7 | BRD_CsSTE11_CAND013 | CsSTE11  |                      100     |                      100    |                     100    |                              100    |                             100    |                            100    |                        19 |                                    100    |
|                     12 | BRD_CsFUS3_CAND003  | CsFUS3   |                       99.947 |                       99.5  |                     100    |                               98.95 |                              90    |                            100    |                        17 |                                    100    |
|                     23 | BRD_CsPKS1_CAND010  | CsPKS1   |                       99.942 |                       98.89 |                     100    |                               98.68 |                              75    |                            100    |                        18 |                                    100    |
|                     38 | BRD_CsPKS1_CAND013  | CsPKS1   |                       99.913 |                       98.9  |                      99.44 |                               99.18 |                              86.88 |                             86.88 |                        17 |                                     88.14 |
|                     50 | BRD_CsSp1_CAND008   | CsSp1    |                       99.942 |                       98.89 |                      98.89 |                               98.91 |                              79.38 |                             79.38 |                        18 |                                     77.78 |

## Threshold Sensitivity Scenarios
| scenario_id                             | evaluated_metric                  |   threshold_pct | mode         |   passing_candidates_n |   CsSTE11_pass |   CsPKS1_pass |   CsSp1_pass |   CsFUS3_pass | ste11_cand012_pass   | pks1_cand013_pass   | cssp1_cand008_pass   |
|:----------------------------------------|:----------------------------------|----------------:|:-------------|-----------------------:|---------------:|--------------:|-------------:|--------------:|:---------------------|:--------------------|:---------------------|
| Scenario_A1_Mean21mer_ge85pct           | mean_21mer_preservation_all19_pct |              85 | mean         |                     48 |             12 |            13 |           10 |            13 | True                 | True                | True                 |
| Scenario_A2_Mean21mer_ge90pct           | mean_21mer_preservation_all19_pct |              90 | mean         |                     46 |             10 |            13 |           10 |            13 | True                 | True                | True                 |
| Scenario_A3_Mean21mer_ge95pct           | mean_21mer_preservation_all19_pct |              95 | mean         |                     44 |              8 |            13 |           10 |            13 | True                 | True                | True                 |
| Scenario_B1_MinNtID_All19_ge85pct       | min_nt_identity_all19_pct         |              85 | all_19_min   |                     50 |             14 |            13 |           10 |            13 | True                 | True                | True                 |
| Scenario_B2_MinNtID_All19_ge90pct       | min_nt_identity_all19_pct         |              90 | all_19_min   |                     50 |             14 |            13 |           10 |            13 | True                 | True                | True                 |
| Scenario_B3_MinNtID_All19_ge95pct       | min_nt_identity_all19_pct         |              95 | all_19_min   |                     50 |             14 |            13 |           10 |            13 | True                 | True                | True                 |
| Scenario_C1_MinNtID_HQ17_ge90pct        | min_nt_identity_hq17_pct          |              90 | hq17_min     |                     50 |             14 |            13 |           10 |            13 | True                 | True                | True                 |
| Scenario_C2_MinNtID_HQ17_ge95pct        | min_nt_identity_hq17_pct          |              95 | hq17_min     |                     50 |             14 |            13 |           10 |            13 | True                 | True                | True                 |
| Scenario_C3_Min21mer_HQ17_ge95pct       | min_21mer_preservation_hq17_pct   |              95 | hq17_min     |                     29 |              6 |            12 |            0 |            11 | True                 | False               | False                |
| Scenario_D1_90pctOfIsolates_ge95pctNtID | isolates_ge95pct_nt_id_pct        |              90 | pct_isolates |                     50 |             14 |            13 |           10 |            13 | True                 | True                | True                 |

## Classification State Breakdown (All 950 Candidate-Isolate Pairs)
| target   |   Confirmed_Divergence |   Confirmed_Match |   Unresolved_Assembly_Polishing_Artifact |
|:---------|-----------------------:|------------------:|-----------------------------------------:|
| CsFUS3   |                      3 |               226 |                                       18 |
| CsPKS1   |                      4 |               233 |                                       10 |
| CsSTE11  |                    135 |               120 |                                       11 |
| CsSp1    |                     10 |               180 |                                        0 |
