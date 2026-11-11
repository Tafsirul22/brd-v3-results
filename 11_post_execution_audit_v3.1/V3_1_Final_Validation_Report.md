# Bipolaris-RNAi Designer (BRD) v3.1 — Final Validation Report

**Generated:** 2026-11-11 08:30:13 UTC  
**Branch / Release:** `v3.1` (`11_post_execution_audit_v3.1/`)  
**Scope Notice:** Computational target evaluation, orthology, conservation, off-target screening, and heuristic prioritization only. No wet-lab silencing, exogenous dsRNA uptake in *Bipolaris sorokiniana*, crop protection, or environmental biosafety is claimed.

---

## 1. Executive summary
BRD v3.1 upgrades the reproducible v3 pipeline across software audit, native Windows 11 execution (`Broccoli v2.0.1`, `Biopython 1.86` assembly QC), dataset verification (`pleosporales_odb12`, `fungi_odb12`, reference FASTA SHA-256s, exact 20/21 mismatch positions), dependency safety (guarded optional `ViennaRNA`), and clean-extraction verification (`hydrate_inputs_v3_1.py` + synchronized `hashes_v3_1.tsv`). Core biological conclusions and the 8-region clean candidate shortlist remain unchanged (`BRD_CsSTE11_CAND013`, `0.9418`, Pareto Tier 1), while `Broccoli v2.0.1` newly resolves true 1:1 orthology for `CsFUS3` (`OG_1315`) and `CsPKS1` (`OG_4251`) and empirical `compleasm run` limitations on native Windows are logged without synthetic scores.

## 2. v3 baseline findings
Stage 0 (`stage0_audit/V3_Baseline_Audit_Report.md`) preserved all original v1 (`10_final_reports/`), v2 (`11_post_execution_audit/`), and v3 (`11_post_execution_audit_v3/`) artifacts and confirmed all `10` independent audit issues (`AUDIT_01`–`AUDIT_10`), notably missing input genomes inside `11_post_execution_audit_v3.zip` (`AUDIT_05`), `hashes_v3.tsv` including 3 skipped Stage 2 files (`AUDIT_07`), and unexecuted Linux binaries (`AUDIT_01`–`04`).

## 3. v3.1 changes
See full auditable details in `V3_1_Change_Log.md`. Principal additions: `assembly_qc_biopython.py`, `stage4_compleasm_runner.py`, native compiled `tools/bin/broccoli.exe` (`v2.0.1`) + `stage5_broccoli_orthology.py`, `stage6_conservation_v31.py`, `stage7_offtarget_v31.py`, `stage8_to_11_v31.py`, `hydrate_inputs_v3_1.py`, `check_environment_v3_1.py`, Docker/Remote Linux templates, and `test_brd_v3_1.py`.

## 4. Tool investigation results
| Repository                                | Version                                                                           | License                                                            | Status                                                                                                                                                 |
|:------------------------------------------|:----------------------------------------------------------------------------------|:-------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------|
| https://github.com/biopython/biopython    | Installed 1.86 (upstream latest tag: biopython-188)                               | Biopython License Agreement / 3-Clause BSD (dual-licensed modules) | PASSED (Native Windows 11 Route A — all 19 isolates parsed & unit-tested; not QUAST)                                                                   |
| https://github.com/huangnengCSU/compleasm | v0.2.10 (commit a275e9f4, released 2026-10-07)                                    | Apache-2.0 (+ BUSCO MIT / CC-BY-ND-4.0 dataset license)            | BLOCKED on native Windows analysis (CLI/list/download PASSED rc=0; run BLOCKED rc=1 missing Win32 miniprot/hmmsearch; Docker/WSL2 absent)              |
| https://github.com/rderelle/Broccoli      | v2.0.1 Rust rewrite (commit e54756aa, released 2026-10-01); v1.4.0 Python audited | GPL-3.0-or-later                                                   | PASSED (Native Windows 11 Route A with -phylogenies nj; default -phylogenies bionj FastTree mode BLOCKED due to upstream MicrobesOnline anti-bot gate) |

Execution routes verified (`execution_environment/execution_routes_status.tsv`):
| route                     | availability                       | blocked_tools                                                                                                             |
|:--------------------------|:-----------------------------------|:--------------------------------------------------------------------------------------------------------------------------|
| Route_A_Native_Windows_11 | ACTIVE                             | compleasm run (missing Win32 miniprot/hmmsearch), Broccoli bionj (FastTree blocked), QUAST, BUSCO, OrthoFinder, RNAhybrid |
| Route_B_Docker_Desktop    | UNAVAILABLE (docker CLI missing)   | All container targets until Docker Desktop is installed                                                                   |
| Route_C_Remote_Linux_SSH  | UNCONFIGURED (no remote host/keys) | All remote Linux targets until SSH endpoint is supplied                                                                   |

## 5. Genome QC results (Biopython — Not QUAST)
All `19/19` *B. sorokiniana* assemblies parsed cleanly (`0` duplicate IDs, `0` duplicate contig sequences, `0` invalid nucleotide bases). Assemblies span `32.42–38.20 Mb`, N50 `130.08 kb–2.31 Mb`, L50 `7–73`. Two scaffolded assemblies (`ND90Pr` `3.4774%` N, `Sacc._Shoemaker` `3.4789%` N) triggered `ELEVATED_AMBIGUOUS_N_PCT` warnings (`assembly_qc_warnings_v3_1.tsv`). Intrinsic Biopython statistics **do not** detect misassemblies or local consensus homopolymers (`WAI2406`, `WAI2411`).

## 6. Compleasm results
| lineage            | orthodb_version   |   marker_count | verified_tar_md5                 | miniprot_version                                                 | hmmer_version                                                         |
|:-------------------|:------------------|---------------:|:---------------------------------|:-----------------------------------------------------------------|:----------------------------------------------------------------------|
| pleosporales_odb12 | OrthoDB_v12.1     |           6066 | 795a5af234fb9b7828b9927dfbebff4a | MISSING (Linux x64 binary only; unavailable natively on Windows) | MISSING (HMMER 3.4 POSIX-only build; unavailable natively on Windows) |
| fungi_odb12        | OrthoDB_v12.1     |           1122 | 79aba415d64ed9a7b1994452b99106cc | MISSING (Linux x64 binary only; unavailable natively on Windows) | MISSING (HMMER 3.4 POSIX-only build; unavailable natively on Windows) |

- **Empirical result:** `compleasm.py v0.2.10` CLI (`-h`, `list --remote`, `download`) succeeded on native Windows 11, downloading and verifying `pleosporales_odb12` (`6,066` markers) and `fungi_odb12` (`1,122` markers). Executing `compleasm.py run` on `ND90Pr` exited `rc=1` due to missing native Windows `miniprot` and `hmmsearch` (`compleasm_execution_log.txt`). All `19` rows in `compleasm_summary_v3_1.tsv` remain `NOT_RUN_DEPENDENCY_BLOCKED` (`NA` scores; never substituted or fabricated). Separate 50-marker tBLASTn core frame audit remains referenced in `compleasm_marker_details_v3_1.tsv`.

## 7. Broccoli results
- Compiled **Broccoli v2.0.1** (`GPL-3.0-or-later`, Rust `1.99.1` `x86_64-pc-windows-gnu`) natively to `tools/bin/broccoli.exe` and executed Steps 1–4 with `diamond.exe v2.2.8` using built-in `-phylogenies nj` (smoke run `4.07s`, full 36,714-protein 3-species run `255.98s`, `rc=0`).
- **Results:** `9,589` orthologous groups (`8,793` shared between `ND90Pr` and `Sacc._Shoemaker`; `5,710` across all 3 proteomes; `26,745` orthologous pairs; `0` chimeric proteins flagged). **Broccoli orthogroups are strictly labelled Broccoli (`OG_*`), never OrthoFinder.**

## 8. DIAMOND RBH comparison
| metric                                                             | analysis_a_diamond_rbh_v3   | analysis_b_broccoli_v2_nj   |
|:-------------------------------------------------------------------|:----------------------------|:----------------------------|
| Total orthogroups inferred                                         | 15015                       | 9589                        |
| Core orthogroups in BOTH B. sorokiniana proteomes (ND90Pr + Sacc.) | 7245                        | 8793                        |
| Three-species conserved orthogroups (+ Botrytis cinerea)           | 4430                        | 5710                        |
| Proteins assigned to any group (unique IDs)                        | 36714                       | 28051                       |
| Broccoli orthologous pairs (step 4)                                | NA (connected components)   | 26745                       |
| Broccoli ortholog pairs co-clustered in same DIAMOND RBH OG        | 22915                       | 22915/26745 (85.68%)        |
| Chimeric proteins flagged by Broccoli                              | nan                         | 0                           |

| target   | nd90pr_protein_id   | broccoli_og   | diamond_rbh_og   |
|:---------|:--------------------|:--------------|:-----------------|
| CsPKS1   | XP_007696213.1      | OG_4251       | RBH_OG000102     |
| CsFUS3   | XP_007697446.1      | OG_1315       | RBH_OG000002     |
| CsSTE11  | XP_007694574.1      | OG_1786       | RBH_OG002545     |
| CsSp1    | XP_007704210.1      | OG_9092       | RBH_OG008444     |

Whereas DIAMOND RBH graph clustering merged `CsFUS3` and `CsPKS1` into large kinase (`RBH_OG000002`) and PKS (`RBH_OG000102`) paralog superfamilies, Broccoli's tree/network partitioning isolated clean 1:1 orthologs (`OG_1315`, `OG_4251`). Concordance (`85.68%` = `22,915/26,745` pairs co-clustered) confirms broad agreement without constituting independent biological proof.

## 9. Conservation results
| target   |   isolates_complete_cds |   msa_min_nt_identity_vs_nd90pr_pct |   msa_ambiguous_n_bases_total |   candidate_level_min_nt_id_hq17_pct |   candidate_level_min_21mer_hq17_pct |
|:---------|------------------------:|------------------------------------:|------------------------------:|-------------------------------------:|-------------------------------------:|
| CsPKS1   |                      17 |                             99.8452 |                             0 |                                99.44 |                                86.88 |
| CsFUS3   |                      17 |                             99.9053 |                             0 |                                99.5  |                                90.87 |
| CsSTE11  |                      17 |                             99.5832 |                             0 |                                99    |                                78.89 |
| CsSp1    |                      19 |                             99.2593 |                             0 |                                98.89 |                                76.67 |

Full 50x19 identity/21-mer matrices and lowest-identity isolate labels (`LK93` reducing `CsSp1_CAND008` to `79.38%` HQ-17 21-mer preservation; `WAI2406`/`WAI2411` homopolymer indels isolated) are reported in `stage6_conservation/`.

## 10. Off-target results
Verified local FASTA SHA-256 checksums across all 18 reference panel entries (`stage7_offtarget/offtarget_reference_metadata.tsv`), resolved strand/coordinates across `1,214` exact 21-mer hit rows (`18/50` excluded, `32/50` pass), and resolved full 21-nt mismatch positions across `16,073` near-hits (`20/21`: `3,655`; `19/21`: `7,501`; `18/21`: `4,917`). Wheat, barley, and maize nuclear non-coding genomes (`>20 GB`) and unsequenced field microbes remain unscreened. **Zero exact matches never establish biological safety, and `7/8` shortlist leaders carry elevated `20/21` near-match burdens (`6–14` guides).**

## 11. Candidate-region results
`49/50` candidates match contiguous authoritative target CDS coordinates/substrings on `BRD_Target_Gene_Sequences.fasta` cleanly (`1` candidate, `BRD_CsSTE11_CAND009` from `dsRNAmax`, is a stitched `SYNTHETIC_KMER_CONSTRUCT`), and all `50` map to `12` independent non-redundant target regions (`8` clean passing). Sequence-derived descriptors, `ViennaRNA 2.8.0` structural metrics, blocked `RNAhybrid` status, and heuristic targetability scores remain cleanly separated (`stage8_candidate_design/candidate_master_v3_1.tsv`).

## 12. Biological target evidence
| gene_id   |   biological_relevance_score | off_target_concerns                                                                                                                                                                                                      | confidence_level                                                                                               |
|:----------|-----------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------|
| CsSTE11   |                       0.9921 | Region 04 (2226..2525) fails exact 21-mer screening in beneficial fungi; Regions 01/03 low near-match signal (3-4 20/21 hits); Region 02 has 14 20/21 near-hits                                                          | HIGH (direct B. sorokiniana KO + PHI curation + 1:1 Broccoli orthology)                                        |
| CsPKS1    |                       0.5796 | 0 exact non-target 21-mers across 4 candidate regions, but ALL 4 regions carry Elevated_20of21_Near_Match_Burden (6–12 unique 20/21 guides)                                                                              | MODERATE-HIGH for RNAi tractability, LOW for disease suppression (direct KO proves unaffected pathogenicity)   |
| CsSp1     |                       0.55   | 0 exact 21-mers across panel, 8 unique 20/21 near-matches; short 270-nt ORF leaves no alternate non-overlapping region                                                                                                   | MODERATE (strong direct wheat-leaf KO & valid 3B/4B/5B induction, but missing PHI/InterPro & SNP drop in LK93) |
| CsFUS3    |                       0.9921 | SEVERE BENEFICIAL-FUNGI & CROP HOMOLOGY: 13/13 candidates excluded by >=1 exact 21-mer match across Trichoderma/Metarhizium/Beauveria/Rice/Wheat/Maize; exhaustive CDS sweep shows only 5/1059 windows escape exact hits | HIGH biological relevance / EXCLUDED for selective dsRNA design across existing windows                        |

- `CsSTE11` leads (`0.9921`); `CsPKS1` preserves `PHI:2260` `unaffected pathogenicity`; `CsSp1` retains `3B/4B/5B` expression (`max = 139.61 RPKM`, invalid `1,438.44` excluded); `CsFUS3` retains beneficial-fungi homology exclusion.

## 13. Ranking results
`independent_region_ranking_v3_1.tsv` orders regions by `(passes_zero_exact_21mer_filter DESC, pareto_front_3d_regions ASC, heuristic_composite_prioritization_index DESC)` whereas pure composite ranking ignores Pareto tier dominance — explaining why Tier-1 `BRD_CsPKS1_CAND010` (`0.7762`) places `#5` ahead of Tier-2 `BRD_CsPKS1_CAND013` (`0.7841`, `#6`).

## 14. Literature benchmarking
Audited `25` curated benchmark studies (`14` gene-level only, `5` ambiguous/conflicting, `4` computational only, `2` reconstructed *P. pachyrhizi* EST sequences) plus `8` exploratory *F. graminearum* constructs (`stage11_literature_benchmark/literature_benchmark_v3_1.tsv`). *B. sorokiniana* still lacks any published exogenous dsRNA SIGS uptake assay.

## 15. Test results
Automated suite summary recorded in `V3_1_Test_Results.txt` and `V3_1_Test_Summary.md` (covering both workspace regression tests and independent clean-extraction execution).

## 16. Reproducibility assessment
Clean extractions can verify all packaged SHA-256 hashes immediately (`hashes_v3_1.tsv` synchronized 1:1 with ZIP entries), inspect environment readiness (`check_environment_v3_1.py`), and hydrate missing NCBI genomes with verified SHA-256 checksums (`hydrate_inputs_v3_1.py --primary-only` or all 19).

## 17. Remaining limitations
1. Linux-only external tools (`compleasm run` needing `miniprot`/`hmmsearch`, `QUAST`, `BUSCO`, `OrthoFinder`, `RNAhybrid`, `FastTree bionj`) remain blocked natively on Windows 11 without Docker or WSL2.
2. Computational rankings do not demonstrate exogenous dsRNA uptake, in vivo knockdown, or disease reduction in *B. sorokiniana*.
3. Whole nuclear non-coding genomes of wheat, barley, and maize remain outside local off-target coverage.
4. `7/8` clean independent-region leaders retain elevated `20/21` near-match burdens (`6–14` guides) alongside `CAND013` (`3`) and `CAND014` (`4`).

## 18. Final computational candidate shortlist
|   independent_region_rank_v3_1 |   composite_index_rank_within_filter_tier | candidate_id        | target   | region_id     |   candidate_length_nt |   biological_target_relevance_score |   dsrna_targetability_score |   sequence_selectivity_score |   heuristic_composite_prioritization_index |   pareto_front_3d_regions |   near_20of21_guides_total | selectivity_classification_v3     | rank_divergence_explanation                                                                                                                                    |
|-------------------------------:|------------------------------------------:|:--------------------|:---------|:--------------|----------------------:|------------------------------------:|----------------------------:|-----------------------------:|-------------------------------------------:|--------------------------:|---------------------------:|:----------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                              1 |                                         1 | BRD_CsSTE11_CAND013 | CsSTE11  | CsSTE11_REG01 |                   210 |                              0.9921 |                      0.8855 |                        0.955 |                                     0.9418 |                         1 |                          3 | Low_20of21_Near_Match_Signal      | Region Pareto-first and composite-index orders agree                                                                                                           |
|                              2 |                                         2 | BRD_CsSTE11_CAND014 | CsSTE11  | CsSTE11_REG03 |                   210 |                              0.9921 |                      0.8859 |                        0.94  |                                     0.9367 |                         1 |                          4 | Low_20of21_Near_Match_Signal      | Region Pareto-first and composite-index orders agree                                                                                                           |
|                              3 |                                         3 | BRD_CsSTE11_CAND010 | CsSTE11  | CsSTE11_REG02 |                   180 |                              0.9921 |                      0.8913 |                        0.79  |                                     0.8861 |                         1 |                         14 | Elevated_20of21_Near_Match_Burden | Region Pareto-first and composite-index orders agree                                                                                                           |
|                              4 |                                         4 | BRD_CsPKS1_CAND012  | CsPKS1   | CsPKS1_REG02  |                   210 |                              0.5796 |                      0.8879 |                        0.895 |                                     0.7979 |                         1 |                          7 | Elevated_20of21_Near_Match_Burden | Region Pareto-first and composite-index orders agree                                                                                                           |
|                              5 |                                         7 | BRD_CsPKS1_CAND010  | CsPKS1   | CsPKS1_REG04  |                   180 |                              0.5796 |                      0.9009 |                        0.82  |                                     0.7762 |                         1 |                         12 | Elevated_20of21_Near_Match_Burden | Pareto Tier 1 non-dominated trade-off prioritized ahead of Tier 2 higher weighted average (CAND010 composite 0.7762 Tier-1 vs CAND013 composite 0.7841 Tier-2) |
|                              6 |                                         5 | BRD_CsPKS1_CAND013  | CsPKS1   | CsPKS1_REG03  |                   180 |                              0.5796 |                      0.8335 |                        0.91  |                                     0.7841 |                         2 |                          6 | Elevated_20of21_Near_Match_Burden | Pareto Tier 1 non-dominated trade-off prioritized ahead of Tier 2 higher weighted average (CAND010 composite 0.7762 Tier-1 vs CAND013 composite 0.7841 Tier-2) |
|                              7 |                                         6 | BRD_CsPKS1_CAND006  | CsPKS1   | CsPKS1_REG01  |                   250 |                              0.5796 |                      0.8772 |                        0.865 |                                     0.7836 |                         2 |                          9 | Elevated_20of21_Near_Match_Burden | Region Pareto-first and composite-index orders agree                                                                                                           |
|                              8 |                                         8 | BRD_CsSp1_CAND008   | CsSp1    | CsSp1_REG01   |                   180 |                              0.55   |                      0.7524 |                        0.88  |                                     0.7363 |                         3 |                          8 | Elevated_20of21_Near_Match_Burden | Region Pareto-first and composite-index orders agree                                                                                                           |
