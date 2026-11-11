# BRD v3.1.1 Test Verification Summary

**Timestamp:** 2026-11-11 09:28:31 UTC  
**OS:** `Windows-11-10.0.26300-SP0`  
**Python:** `3.13.14` (`C:\Users\ASUS BANGLADESH\Desktop\New folder (2)\Output\environments\brd_venv\Scripts\python.exe`)  
**Key Package Versions:** `biopython=1.86`, `pandas=2.3.3`, `numpy=2.4.0`, `scipy=1.16.3`, `pytest=9.1.1`, `ViennaRNA=2.8.0`, `primer3-py=2.4.0`  
**PDF Engine:** `Edge_headless (300041 B)`

## 1. Full Workspace Regression Suite
- **Command:** `pytest 11_post_execution_audit_v3.1.1/tests/test_brd_v3_1_1.py 11_post_execution_audit_v3.1/tests/test_brd_v3_1.py 11_post_execution_audit_v3/tests/test_brd_v3.py tests/test_brd.py 11_post_execution_audit/tests/test_audit_v2.py -v`
- **Observed Outcome:** **`53 passed in 51.14s`** (`Exit Code: 0`, `0` failed, `0` skipped, `0` blocked in workspace)
  - `11_post_execution_audit_v3.1.1/tests/test_brd_v3_1_1.py` (`12` tests, includes isolated Primer3-absent, live Primer3-installed, isolated optional-all-absent, Fix 2/3/5/6 tests)
  - `11_post_execution_audit_v3.1/tests/test_brd_v3_1.py` (`11` tests)
  - `11_post_execution_audit_v3/tests/test_brd_v3.py` (`9` tests)
  - `tests/test_brd.py` (`8` tests)
  - `11_post_execution_audit/tests/test_audit_v2.py` (`13` tests, relative paths verified)

## 2. Clean-Extraction Test Coverage Matrix (`V3_1_1_Test_Inventory.tsv`)
- **Immediately runnable from compact ZIP without optional deps or large downloads:** `48/53` tests (`11/12` v3.1.1 + `11/11` v3.1 + `8/8` core + `12/13` v2 + `6/9` v3).
- **Conditional / Blocked cases in bare compact clean extraction:**
  - `test_brd_v3_1_1.py::test_02b_primer3_installed_live_design` requires optional `primer3-py>=2.0` (`SKIPPED` cleanly if omitted).
  - `test_brd_v3.py::test_stage4_expanded_panel_and_selectivity_invariants` requires optional `ViennaRNA>=2.7.0`.
  - `test_audit_v2.py::test_11_both_orientations_screened_and_all_exact_21mers_reported` requires beneficial-fungi references (`hydrate_offtarget_references_v3_1_1.py`).
  - `test_brd_v3.py::test_stage1_authoritative_invariants` requires all `19` genomes (`hydrate_inputs_v3_1.py`).
  - `test_brd_v3.py::test_stage8_report_18_sections_and_hashes` remains **BLOCKED in compact ZIP by design** because v3's `hashes_v3.tsv` hashes excluded multi-MB `.dmnd` / `diamond_all_vs_all_v3.tsv` caches (resolved cleanly for v3.1.1 in `hashes_v3_1_1.tsv`).
