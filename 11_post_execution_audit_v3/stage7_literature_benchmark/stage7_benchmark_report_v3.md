# Stage 7 — Retrospective Literature-Based Benchmark Analysis (BRD v3)

**Timestamp:** 2026-11-11 06:55:13 UTC

## 1. Benchmark Dataset Summary (`literature_benchmark_dataset_v3.tsv`)
| category                           |   records |
|:-----------------------------------|----------:|
| 3_Gene_Level_Only                  |        14 |
| 5_Ambiguous_or_Conflicting_Outcome |         5 |
| 4_Computational_Only               |         4 |
| 2_Reconstructed_Sequence           |         2 |

## 2. Qualitative Structured Benchmark Findings (Supporting & Refuting BRD v3 Design Principles)
1. **Why BRD v3 rejects single-metric target ranking (`PMC5943259`, `PMC13205770`, `PMC6465682`):** High expression (`Ssoah` FPKM >5,200; `SS1G_10167` FPKM ~24,000) and knockout essentiality (`Ssoah`) did NOT guarantee disease suppression under SIGS (`SS1G_06305` worsened infection by `129%`). Similarly, 3 of 4 `MoAP1/MoACT/MoSSADH/MoSOM1` asiRNAs failed to affect growth in *M. oryzae*, and 33% of `LUCID`-designed targets failed in-planta disease suppression despite mRNA knockdown.
2. **Why BRD v3 favors longer `180–550 nt` multi-guide dsRNAs & multi-gene options (`PMC7202221`, `PMC6859480`, `PMC13608320`):** Within `Werner et al. (2020)` (`n=8` *F. graminearum* constructs), length correlated positively with silencing (exploratory Spearman `rho=0.619, p=0.102`), while two tool-selected short `173-nt` constructs produced zero `FgAGO1` silencing. Likewise, dual/triple constructs outperformed single targets in *F. graminearum* (`CYP3RNA` 93% vs 78-82%) and Dothideomycete *A. alternata* (`AaGH10+AaSOD` SF > 1).
3. **Why BRD v3 requires explicit exogenous-uptake caution (`PMC8428832`, `PMC7391711`, `PMC4443973`):** `Leng et al. (2011)` proved intracellular hairpin RNAi works in *B. sorokiniana*, but **exogenous dsRNA uptake** (`Qiao 2021` uptake gate) has never been demonstrated in *B. sorokiniana* and must be tested empirically before assuming field SIGS viability.

## 3. Within-Study Exploratory Analysis (`Werner et al. 2020`, `n=8` Constructs)
| construct_group   | target   | design_mode   |   length_nt |   efficient_sirnas |   avg_silencing_pct |
|:------------------|:---------|:--------------|------------:|-------------------:|--------------------:|
| FgAGO1_Tool       | FgAGO1   | Tool          |         173 |                 13 |                   6 |
| FgAGO1_Manual     | FgAGO1   | Manual        |         658 |                 57 |                  32 |
| FgAGO2_Tool       | FgAGO2   | Tool          |         192 |                 12 |                  46 |
| FgAGO2_Manual     | FgAGO2   | Manual        |         871 |                 58 |                  49 |
| FgDCL1_Tool       | FgDCL1   | Tool          |         182 |                  5 |                  60 |
| FgDCL1_Manual     | FgDCL1   | Manual        |         912 |                 49 |                  77 |
| FgDCL2_Tool       | FgDCL2   | Tool          |         193 |                  9 |                  44 |
| FgDCL2_Manual     | FgDCL2   | Manual        |         870 |                 92 |                  63 |
- Length vs silencing: Spearman `rho = 0.6190` (`p = 1.0173e-01`)
- SiFi efficient siRNAs vs silencing: Spearman `rho = 0.2143` (`p = 6.1034e-01`)
*(Exploratory within-study summary only — n=8 from one detached-leaf barley experiment cannot calibrate predictive cross-species RNAi accuracy.)*
