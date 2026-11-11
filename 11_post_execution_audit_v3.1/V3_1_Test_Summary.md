# BRD v3.1 Test Verification Summary

**Timestamp:** 2026-11-11 08:31:45 UTC  
**PDF Engine:** Edge_headless (rendered `V3_1_Final_Validation_Report.pdf`)

## 1. Workspace Full Regression Run
- **Command:** `pytest 11_post_execution_audit_v3.1/tests/test_brd_v3_1.py 11_post_execution_audit_v3/tests/test_brd_v3.py tests/test_brd.py 11_post_execution_audit/tests/test_audit_v2.py -v`
- **Result:** **`41 passed in 46.48s`** (`Exit Code: 0`, `0` failed, `0` skipped)
  - `11_post_execution_audit_v3.1/tests/test_brd_v3_1.py`: `11 passed`
  - `11_post_execution_audit_v3/tests/test_brd_v3.py`: `9 passed`
  - `tests/test_brd.py`: `8 passed`
  - `11_post_execution_audit/tests/test_audit_v2.py`: `13 passed`

## 2. Independent Clean-Environment Extraction Validation (`11_post_execution_audit_v3.1.zip`)
- **ZIP Verification:** `zipfile.ZipFile('11_post_execution_audit_v3.1.zip').testzip()` -> `None` across `340` entries (`7.72 MB` compressed release archive; `0` CRC errors).
- **Clean-Extraction Pytest (`cwd=%TEMP%\brd_v31_clean_extract`):**
  - Command: `pytest 11_post_execution_audit_v3.1/tests/test_brd_v3_1.py tests/test_brd.py -v`
  - Result: **`19 passed in 5.14s`** (`Exit Code: 0`, `0` failed, `0` skipped; includes `100%` SHA-256 verification of `hashes_v3_1.tsv`).
- **Optional `ViennaRNA` Import Guard Check:** Passed cleanly when `RNA` and `sifi21_plus` are masked as unavailable (`OPTIONAL_VIENNARNA_ABSENT_IMPORT_GUARD_OK`).
- **Live Genome Hydration Check (`--primary-only`):** `hydrate_inputs_v3_1.py --primary-only` downloaded `GCF_000338995.1_Cocsa1_genomic.fna` from NCBI Datasets v2 (`total_requested=1, verified=1, downloaded=1, failed=[]`) and `analyze_fasta_assembly` returned `OK total_length_bp=34409167 n50_bp=1789485 l50=7 warnings=1`.

## 3. Honest Tool Execution & Blocker Inventory
- **Executed natively on Windows 11:** `Biopython 1.86` assembly QC (`19/19`), `Broccoli v2.0.1` (`-phylogenies nj` + `diamond v2.2.8` across `36,714` proteins, `rc=0`), `compleasm.py v0.2.10` CLI & OrthoDB `v12.1` downloads (`pleosporales_odb12`, `fungi_odb12`).
- **Blocked (NOT_RUN_DEPENDENCY_BLOCKED — never fabricated or skipped silently):** `compleasm run` (requires Linux/POSIX `miniprot` and `hmmsearch`), `QUAST`, `BUSCO`, `OrthoFinder`, `RNAhybrid`, and `FastTree` BioNJ mode.
