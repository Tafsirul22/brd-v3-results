# BRD v3.1 to v3.1.1 Focused Remediation Change Log

**Release Date:** 2026-11-11 09:28:31 UTC  
**Authoritative v3.1 Baseline:** `11_post_execution_audit_v3.1.zip` (SHA-256 `bd05ecc6fc8553ee99687020863e2008e1d1477e99f7e0e92299d09748b30066`)

## 1. Code Fixes
- **`brd/integrations/lucid.py` (`FIX_1`):** Guarded `import primer3`, `import RNA`, and `import sifi21_plus.free_energy` at import time; exposed `PRIMER3_AVAILABLE`, `VIENNARNA_AVAILABLE`, `SIFI21_AVAILABLE`, and `is_primer3_available()`; isolated `design_amplification_primers(...)` and added `design_primers: bool = True` to `design_lucid_windows(...)` so invoking Primer3 without `primer3-py` raises an actionable `RuntimeError` while unrelated RNAi modules and primer-free window scoring remain operational.
- **`brd/modules/m07_dsRNA_design.py` & `brd/modules/m10_reporting.py` (`FIX_1`):** Removed unused top-level `import primer3` from `m07_dsRNA_design.py` and guarded optional imports (`RNA`, `primer3`, `snakemake`) inside `_collect_software_versions()` in `m10_reporting.py`.
- **`11_post_execution_audit/tests/test_audit_v2.py` (`FIX_4` Action 10):** Replaced hardcoded developer machine `OUT_DIR` path (pre-edit SHA-256 `77d27274ad381715506add838ea0bbdc50670796c5889846543731f2442e24f9`) with relative `Path(__file__).resolve().parents[2]` plus `_resolve_relative_manifest_fasta(...)` resolving extracted repo / hydrated Stage 3 paths only; all 13 v2 assertions preserved untouched.
- **New v3.1.1 scripts:** `stage0_inspect_v31.py`, `hydrate_broccoli_proteomes_v3_1_1.py`, `hydrate_offtarget_references_v3_1_1.py`, `stage6_coordinates_audit_v311.py`, `check_environment_v3_1_1.py`, `release_layout_v3_1_1.py`, `verify_manifest_v3_1_1.py`, `stage8_to_9_release_v311.py`, and `tests/test_brd_v3_1_1.py`.

## 2. Packaging Fixes
- **Restored historical test assets inside `11_post_execution_audit_v3.1.1.zip` (`FIX_4`):** Added `11_post_execution_audit/tests/test_audit_v2.py`, required v2 deliverables (`11_post_execution_audit/`), `10_final_reports/BRD_{dsRNA_Candidates,Target_Gene_Summary}.tsv`, `02_genome_qc/assembly_quality_summary.tsv`, and `08_off_target_screening/blast_databases/database_manifest.tsv` alongside existing v3, v3.1, and core tests.
- **1:1 ZIP vs Manifest bijection (`FIX_5`):** Unified `release_layout_v3_1_1.py` across packaging and manifest verification so `hashes_v3_1_1.tsv` covers `100%` of packaged files (including `V3_1_1_Execution_Log.txt`, `V3_1_1_Test_Results.txt`, `V3_1_1_Test_Summary.md`, and `V3_1_1_Final_Validation_Report.*`) with zero unlisted files and only `hashes_v3_1_1.tsv` itself excluded from self-hashing.

## 3. Dependency Fixes
- **Mandatory vs optional split:** Separated core runtime requirements in `requirements.txt` from optional extras (`ViennaRNA>=2.7.0`, `primer3-py>=2.0`, `snakemake>=9.0`) in `requirements-optional.txt` and `environment.yml`.
- **Dual-state testing:** Tested `primer3-py` installed (`2.4.0` live call) and absent (isolated `sys.meta_path` subprocess without `primer3` and without optional `RNA`/`sifi21_plus`).

## 4. Reproducibility Improvements
- **`hydrate_broccoli_proteomes_v3_1_1.py` (`FIX_2`):** Downloads pinned sources (`GCF_000338995.1`, `GCA_013416765.1`, `Brunxi/LUCID@37910bd` Ensembl `B. cinerea`), reconstructs `Bsor_ND90Pr.faa`, `Bsor_Sacc._Shoemaker.faa`, `Bcin_ASM83294v1.faa`, and `all_proteomes_combined_v3.faa` with `100%` byte-for-byte SHA-256 identity (`4/4` validated), and executes a `150`-protein/species Broccoli smoke test when `broccoli`/`diamond` binaries exist.
- **`hydrate_offtarget_references_v3_1_1.py` (`FIX_3`):** Downloads and validates all `18` panel organisms (`33` FASTAs: `18` transcriptomes + `15` genomes; `33/33` SHA-256 PASS, `0` failed), rebuilds combined FASTA + smoke BLAST DB, and runs dual-strand exact 21-mer screening smoke comparison.

## 5. Documentation Corrections
- **`FIX_6` Coordinate Audit Reconciliation:** Corrected the ambiguous v3.1 change-log wording (which failed to separate the stitched synthetic construct from exact CDS matches) without altering any candidate sequence. Explicit breakdown (`stage6_coordinates/candidate_coordinate_summary_v3_1_1.tsv`):
  - Candidates with recorded coordinate assignments: **`50/50`**
  - Candidates matching contiguous authoritative CDS substrings: **`49/50`**
  - Candidates that are stitched or synthetic constructs: **`1/50`** (`BRD_CsSTE11_CAND009` from `dsRNAmax`, `exact_cds_substring=False`, `is_stitched_or_synthetic=True`, `genomic_start_1based=-1`)
  - Candidates with unresolved coordinate mapping: **`0/50`**

## 6. Scientific Output Changes
- **None.** All v3.1 candidate sequences, `49/50` + `1` coordinate assignments, conservation metrics, off-target baseline hits (`1,214` exact, `16,073` near-hits), Broccoli `v2.0.1` (`9,589` OGs) & DIAMOND RBH (`7,245` core OGs) outputs, biological relevance scores, `8`-candidate clean independent-region shortlist, Pareto tiers, bootstrap results, and historical v2/v3 reports remain **100% unchanged**.
