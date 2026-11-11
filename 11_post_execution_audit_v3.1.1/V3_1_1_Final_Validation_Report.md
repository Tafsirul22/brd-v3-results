# Bipolaris-RNAi Designer (BRD) v3.1.1 — Final Validation Report

**Generated:** 2026-11-11 09:28:31 UTC  
**OS / Python:** `Windows-11-10.0.26300-SP0` / Python `3.13.14`  
**Baseline Archive:** `11_post_execution_audit_v3.1.zip` (`8,089,873` bytes, SHA-256 `bd05ecc6fc8553ee99687020863e2008e1d1477e99f7e0e92299d09748b30066`)

> **Scientific Scope Boundary:** All findings in BRD v3.1.1 remain strictly computational. None of the analyses establish exogenous dsRNA uptake, transcript knockdown, fungicidal efficacy, crop protection, or environmental biosafety in *Bipolaris sorokiniana*.

---

## 1. What Was Fixed (`FIX_1`–`FIX_6`)
| fix_id   | area                   | defect_in_v31                                                                                                                                                                                    | target_resolution                                                                                                                                                       |
|:---------|:-----------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| FIX_1    | Optional Dependency    | Unconditional `import primer3` in `brd/integrations/lucid.py` (and `m07_dsRNA_design.py`, `m10_reporting.py`) crashes in clean env when `primer3-py` absent                                      | Guard optional imports, expose PRIMER3_AVAILABLE, raise actionable RuntimeError when Primer3 called without package, test both installed & absent states                |
| FIX_2    | Broccoli Proteomes     | Broccoli inputs (`Bsor_ND90Pr.faa`, `Bsor_Sacc._Shoemaker.faa`, `Bcin_ASM83294v1.faa`, `all_proteomes_combined_v3.faa`) excluded without hydration script                                        | Implement `hydrate_broccoli_proteomes_v3_1_1.py` fetching exact NCBI/GitHub sources, reconstructing identical files, verifying SHA-256 & Broccoli smoke test            |
| FIX_3    | Off-Target Panel       | Metadata & SHA-256 recorded for 18 panel refs (`17` non-targets + self), but no reference hydration script existed                                                                               | Implement `hydrate_offtarget_references_v3_1_1.py` fetching exact versioned NCBI accessions, verifying length & SHA-256, rebuilding DB & running smoke screen           |
| FIX_4    | Regression Testing     | v3.1 ZIP lacked `test_audit_v2.py` (`test_audit_v2.py` also hardcoded `OUT_DIR`), v2/v1 files, and complete clean-extraction test accounting                                                     | Restore historical tests with relative paths, include required non-huge v2/v1 fixtures, document core vs optional deps & passed vs blocked counts                       |
| FIX_5    | Manifest Completeness  | `hashes_v3_1.tsv` omitted `V3_1_Execution_Log.txt`, `V3_1_Test_Results.txt`, and `V3_1_Test_Summary.md`                                                                                          | Freeze all files & reports before hashing, include ALL release files in `hashes_v3_1_1.tsv`, and verify via standalone script detecting unlisted/missing files          |
| FIX_6    | Coordinate Consistency | `V3_1_Change_Log.md` claimed uniform verification of all 50 candidate coordinates on target FASTA, contradicting report (`49/50` contiguous CDS + `1` synthetic construct `BRD_CsSTE11_CAND009`) | Distinguish recorded coordinate assignments (`50`), contiguous CDS (`49`), stitched constructs (`1`: `CAND009`), and unresolved (`0`) consistently across all artifacts |

## 2. What Was Not Changed (Scientific Preservation)
All valid v2, v3, and v3.1 datasets and scientific results remain untouched (`11_post_execution_audit_v3.1.zip` unmodified; `11_post_execution_audit_v3/`, `11_post_execution_audit/`, `10_final_reports/`, `01_data_acquisition/` data preserved):
- Biopython `1.86` standalone assembly QC on `19/19` isolates (`32.42–38.20 Mb`, N50 `130.08 kb–2.31 Mb`; explicitly not QUAST).
- Broccoli `v2.0.1` (`9,589` orthogroups, `26,745` pairs, 1:1 `CsSTE11`/`CsPKS1`/`CsFUS3`/`CsSp1` OGs) and preserved DIAMOND RBH (`7,245` core *B. sorokiniana* OGs).
- Compleasm `v0.2.10` (`pleosporales_odb12` `6,066`, `fungi_odb12` `1,122`) status preserved honestly as `NOT_RUN_DEPENDENCY_BLOCKED` on Windows (`miniprot`/`hmmsearch` absent natively).
- Multi-isolate conservation matrices, candidate clustering (`12` independent regions, `8` clean), target evidence (`PHI:2260`, `CsSp1` `139.61 RPKM`, `CsFUS3` selectivity warning), and Pareto/sensitivity/bootstrap rankings.
- Preserved v3.1 off-target baseline figures (`1,214` exact-hit rows, `16,073` near-hit rows) labelled strictly as v3.1 baseline results (full BLASTn suite not re-executed).

## 3. Inputs Recovered & Verified
- **Broccoli Proteomes (`stage2_proteome_hydration/`):** `hydrate_broccoli_proteomes_v3_1_1.py` downloaded and verified all `3` raw proteomes (`Bsor_ND90Pr`: `12,210`, `Bsor_Sacc._Shoemaker`: `10,755`, `Bcin_ASM83294v1`: `13,749`) and deterministically reconstructed all `4` files (`36,714` combined sequences, `18,424,497` B) with **`4/4` exact historical SHA-256 matches** (`broccoli_proteome_validation.tsv`).
- **Off-Target Panel (`stage3_offtarget_hydration/`):** `hydrate_offtarget_references_v3_1_1.py` retrieved all **`18/18`** documented panel organisms (`33` FASTAs: `18` transcriptomes + `15` genomes) directly from NCBI Datasets v2 endpoints with **`33/33` exact byte/SHA-256 matches** (`0` failures) in `offtarget_reference_validation.tsv`.
- **19 *B. sorokiniana* Genomes:** Recoverable with exact SHA-256 via `hydrate_inputs_v3_1.py` (`genome_checksums_manifest_v3_1.tsv`).

## 4. Inputs Still Excluded from Compact Archive (Recoverable / Unscreened Boundaries)
- Multi-hundred-MB hydrated FASTAs (19 genomes ~`670` MB, 4 proteome FASTAs ~`36` MB, 18 off-target FASTAs ~`2.1` GB) and v3/v3.1 temporary intermediate caches (`all_proteomes_db_v3.dmnd` `18.9` MB, `diamond_all_vs_all_v3.tsv` `34.1` MB, Broccoli `dir_step1/*.fas`/`dir_step2/`) are intentionally omitted from the compact ZIP and enumerated in `excluded_large_files_v3_1_1.tsv`.
- Multi-GB nuclear non-coding genomes of wheat (`T. aestivum`), barley (`H. vulgare`), and maize (`Z. mays`) were never downloaded or screened in v3/v3.1 (RefSeq transcriptomes only).

## 5. Clean-Extraction Workflow Reproducibility
- **Broccoli (`stage2_proteome_hydration/broccoli_proteome_recovery_methods.md`):** Proteome inputs are **100% reproducibly recoverable** from a clean extraction. Executing `broccoli` itself additionally requires external `broccoli` (`v2.0.1`) and `diamond` (`v2.2.8`) binaries; when present, the `150`-protein/species smoke test passed (`4.31s`, `13` OGs, `17` pairs).
- **Off-Target Screening (`stage3_offtarget_hydration/offtarget_reference_recovery_methods.md`):** All `18` panel references are **100% reproducibly recoverable and SHA-256 verified** via `hydrate_offtarget_references_v3_1_1.py`; combined DB rebuild + `blastn-short` + dual-strand 21-mer smoke test confirmed `BRD_CsSTE11_CAND013` has `0` exact 21-mers while `BRD_CsFUS3_CAND001` matches exact beneficial-fungi hits (`v3.1` baseline counts `1,214` / `16,073` remain labeled historical v3.1 baseline results until full panel screening is explicitly rerun).

## 6. Candidate Coordinate Audit Reconciliation (`FIX_6`)
| metric                                                      |   count |   total | details                                                                                                   |
|:------------------------------------------------------------|--------:|--------:|:----------------------------------------------------------------------------------------------------------|
| candidates_with_recorded_coordinate_assignments             |      50 |      50 | All 50/50 carry 1-based inclusive CDS span [cds_start_1based, cds_end_1based], refseq_mrna_id, and strand |
| candidates_matching_contiguous_authoritative_cds_substrings |      49 |      50 | 49/50 match ND90Pr CDS substrings in BRD_Target_Gene_Sequences.fasta with 0 mismatches                    |
| candidates_stitched_or_synthetic_constructs                 |       1 |      50 | 1/50: BRD_CsSTE11_CAND009 (dsRNAmax, genomic_start=-1, 1 mismatches vs CsSTE11 1..300)                    |
| candidates_with_unresolved_coordinate_mapping               |       0 |      50 | 0/50 unresolved (all 50 explicitly accounted for without altering any candidate sequence)                 |

Every candidate (`50/50`) has 1-based inclusive CDS coordinates (`stage6_coordinates/candidate_coordinate_audit_v3_1_1.tsv`), comprising **`49/50`** exact contiguous substrings of `BRD_Target_Gene_Sequences.fasta` (`0` mismatches) plus **`1/50`** stitched synthetic construct (`BRD_CsSTE11_CAND009`, `dsRNAmax`, `is_contiguous=False`, `genomic_start_1based=-1`, `1` nt difference vs nominal `CsSTE11` `1..300`), with **`0/50`** unresolved mappings.

## 7. Regression Coverage, Test Inventory & Clean-Extraction Accounting
All `5` historical & current test suites (`53` tests total: `13` v2, `8` core, `9` v3, `11` v3.1, `12` v3.1.1) are packaged inside `11_post_execution_audit_v3.1.1.zip` and inventoried in `V3_1_1_Test_Inventory.tsv` (`48/53` runnable immediately from bare compact extraction; remaining cases depend transparently on hydrated FASTAs, optional `primer3-py`/`ViennaRNA`, or v3 `.dmnd`/`diamond_all_vs_all_v3.tsv` files — see `V3_1_1_Test_Summary.md` and `V3_1_1_Test_Results.txt`).

## 8. Release Integrity Manifest Completeness & Clean-Extraction Readiness (`FIX_5`)
`hashes_v3_1_1.tsv` enumerates **every packaged file in `11_post_execution_audit_v3.1.1.zip`** (including `V3_1_1_Execution_Log.txt`, `V3_1_1_Test_Results.txt`, `V3_1_1_Test_Summary.md`, reports, PDF, source, tests, and data), excluding only `hashes_v3_1_1.tsv` itself. Both disk files and the final packaged ZIP are strictly verified at build time via `verify_manifest_v3_1_1.py --zip-file 11_post_execution_audit_v3.1.1.zip` (`0` missing, `0` unlisted, `0` duplicate paths, `0` byte/SHA-256 mismatches, `zipfile.testzip() == None`), enabling clean-extraction verification upon unzip.

## 9. Scientific Rankings & Preserved Shortlist (Unchanged from v3.1)
**Did candidate rankings change? No.**

|   independent_region_rank_v3_1 |   composite_index_rank_overall_12regions |   pareto_front_3d_regions | region_id     | candidate_id        | target   |   heuristic_composite_prioritization_index |   near_20of21_guides_total |
|-------------------------------:|-----------------------------------------:|--------------------------:|:--------------|:--------------------|:---------|-------------------------------------------:|---------------------------:|
|                              1 |                                        1 |                         1 | CsSTE11_REG01 | BRD_CsSTE11_CAND013 | CsSTE11  |                                     0.9418 |                          3 |
|                              2 |                                        2 |                         1 | CsSTE11_REG03 | BRD_CsSTE11_CAND014 | CsSTE11  |                                     0.9367 |                          4 |
|                              3 |                                        3 |                         1 | CsSTE11_REG02 | BRD_CsSTE11_CAND010 | CsSTE11  |                                     0.8861 |                         14 |
|                              4 |                                        4 |                         1 | CsPKS1_REG02  | BRD_CsPKS1_CAND012  | CsPKS1   |                                     0.7979 |                          7 |
|                              5 |                                        7 |                         1 | CsPKS1_REG04  | BRD_CsPKS1_CAND010  | CsPKS1   |                                     0.7762 |                         12 |
|                              6 |                                        5 |                         2 | CsPKS1_REG03  | BRD_CsPKS1_CAND013  | CsPKS1   |                                     0.7841 |                          6 |
|                              7 |                                        6 |                         2 | CsPKS1_REG01  | BRD_CsPKS1_CAND006  | CsPKS1   |                                     0.7836 |                          9 |
|                              8 |                                        8 |                         3 | CsSp1_REG01   | BRD_CsSp1_CAND008   | CsSp1    |                                     0.7363 |                          8 |

## 10. Remaining Limitations
1. Purely computational prioritisation without wet-lab silencing, exogenous dsRNA uptake, or crop protection data.
2. Linux-only binaries (`compleasm run`, `QUAST`, `BUSCO`, `OrthoFinder`, `RNAhybrid`, `FastTree BioNJ`) remain blocked natively on Windows 11 without Docker or WSL2.
3. Whole nuclear non-coding genomes of wheat, barley, and maize remain unscreened, and `7/8` clean shortlist leaders carry `Elevated_20of21_Near_Match_Burden` (`6–14` unique `20/21` guides alongside `3–4` in `CAND013`/`CAND014`).
