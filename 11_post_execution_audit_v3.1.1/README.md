# Bipolaris-RNAi Designer (BRD) v3.1.1

**Release Date:** 2026-11-11 09:28:31 UTC  
**Base Release:** `11_post_execution_audit_v3.1.zip` (`8,089,873` bytes, SHA-256 `bd05ecc6fc8553ee99687020863e2008e1d1477e99f7e0e92299d09748b30066`)

## 1. Overview & Scientific Boundary
BRD v3.1.1 is a focused reproducibility, dependency, testing, manifest-integrity, and coordinate-documentation remediation of BRD v3.1. All valid v2, v3, and v3.1 scientific outputs and candidate rankings are preserved unmodified. All findings remain strictly computational and do **not** constitute wet-lab proof of dsRNA uptake, silencing efficacy, crop protection, or biological safety.

## 2. Clean Installation (Mandatory vs Optional Dependencies)
```powershell
# Mandatory core Python environment (imports brd cleanly without primer3-py or ViennaRNA)
pip install -r 11_post_execution_audit_v3.1.1/requirements.txt

# Optional extras (PCR primer design in lucid.py + ViennaRNA RNA folding / v3 stage4 test)
pip install -r 11_post_execution_audit_v3.1.1/requirements-optional.txt
```

## 3. Input Hydration & Reproducibility Workflow
```powershell
# A. Hydrate primary B. sorokiniana ND90Pr genome (or omit --primary-only for all 19 genomes)
python 11_post_execution_audit_v3.1/src/hydrate_inputs_v3_1.py --primary-only

# B. Reconstruct exact Broccoli v2.0.1 proteomes (36,714 sequences, 4 files SHA-256 verified + smoke test)
python 11_post_execution_audit_v3.1.1/src/hydrate_broccoli_proteomes_v3_1_1.py --smoke-test

# C. Hydrate & SHA-256 verify all 18 off-target panel organisms (32 FASTAs) + rebuild DB + smoke test
python 11_post_execution_audit_v3.1.1/src/hydrate_offtarget_references_v3_1_1.py

# D. Verify release integrity (0 missing, 0 unlisted, 0 size/SHA256 mismatches)
python 11_post_execution_audit_v3.1.1/src/verify_manifest_v3_1_1.py
```

## 4. Regression Testing
```powershell
# Immediate clean-extraction suites (without external multi-hundred-MB downloads):
pytest 11_post_execution_audit_v3.1.1/tests/test_brd_v3_1_1.py 11_post_execution_audit_v3.1/tests/test_brd_v3_1.py tests/test_brd.py -v

# Full regression after hydrating genomes & off-target references (+ optional deps installed):
pytest 11_post_execution_audit_v3.1.1/tests/test_brd_v3_1_1.py 11_post_execution_audit_v3.1/tests/test_brd_v3_1.py tests/test_brd.py 11_post_execution_audit/tests/test_audit_v2.py -v
```
