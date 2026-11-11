# Stage 0 — Supplied v3.1 Archive Inspection & Baseline Snapshot (BRD v3.1.1)

**Inspected:** 2026-11-11 09:28:29 UTC  
**Authoritative Archive:** `11_post_execution_audit_v3.1.zip` (`8,089,873` bytes, SHA-256 `bd05ecc6fc8553ee99687020863e2008e1d1477e99f7e0e92299d09748b30066`, `testzip() = None` across `340` entries)

## Comparison vs Workspace Repository
- `IDENTICAL` entries in extracted archive vs repo: `333`
- `MODIFIED_IN_REPO` (since archive creation — expected only for Stage 1 fixes applied on `v3.1.1`): `7`
- `OMITTED_FROM_V3_1_ZIP` (`input_proteomes/`, `dir_step1/*.fas`, `dir_step2/`, caches): `9`

## Preserved v3.1 Baseline Shortlist (8 Clean Independent Regions)
| Rank | Candidate | Target | Region | Score |
|---:|---|---|---|---:|
| 1 | BRD_CsSTE11_CAND013 | CsSTE11 | CsSTE11_REG01 | 0.9418 |
| 2 | BRD_CsSTE11_CAND014 | CsSTE11 | CsSTE11_REG03 | 0.9367 |
| 3 | BRD_CsSTE11_CAND010 | CsSTE11 | CsSTE11_REG02 | 0.8861 |
| 4 | BRD_CsPKS1_CAND012 | CsPKS1 | CsPKS1_REG02 | 0.7979 |
| 5 | BRD_CsPKS1_CAND010 | CsPKS1 | CsPKS1_REG04 | 0.7762 |
| 6 | BRD_CsPKS1_CAND013 | CsPKS1 | CsPKS1_REG03 | 0.7841 |
| 7 | BRD_CsPKS1_CAND006 | CsPKS1 | CsPKS1_REG01 | 0.7836 |
| 8 | BRD_CsSp1_CAND008 | CsSp1 | CsSp1_REG01 | 0.7363 |
