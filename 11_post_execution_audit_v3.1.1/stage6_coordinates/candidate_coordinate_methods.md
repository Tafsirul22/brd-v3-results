# Candidate Coordinate Audit & Reconciliation Methods (BRD v3.1.1)

**Audited:** 2026-11-11 09:12:49 UTC

## 1. Reconciliation of v3.1 Defect (FIX 6)
The v3.1 `V3_1_Final_Validation_Report.md` accurately reported `49/50` contiguous CDS substring matches plus `1` synthetic construct (`BRD_CsSTE11_CAND009`), whereas `V3_1_Change_Log.md` ambiguously stated `Verified 50/50 CDS coordinates on target FASTA`. All v3.1.1 tables, change logs, and validation reports now consistently separate all four coordinate categories without modifying any candidate sequence.

## 2. Explicit 4-Category Audit Counts
| Category | Count | Details |
|---|---:|---|
| 1. Candidates with recorded coordinate assignments | **50/50** | 1-based inclusive target CDS span `[cds_start_1based, cds_end_1based]` + RefSeq mRNA ID & strand |
| 2. Contiguous authoritative CDS substring matches | **49/50** | Exact contiguous substring on `BRD_Target_Gene_Sequences.fasta` (`0` mismatches, valid genomic coordinates) |
| 3. Stitched or synthetic constructs | **1/50** | `BRD_CsSTE11_CAND009` (`dsRNAmax`), `is_contiguous=False`, `is_stitched_or_synthetic=True`, `genomic_start_1based=-1`, `1` nt difference vs nominal `CsSTE11` `1..300` CDS slice |
| 4. Unresolved coordinate mappings | **0/50** | Every candidate is cleanly classified (`49 + 1 + 0 = 50`) |

## 3. Coordinate Convention & Sequence Domain
- **Convention:** `1-based_inclusive` (corresponding Python slice `[cds_start_1based - 1 : cds_end_1based]`).
- **Domains recorded:** Both spliced target **CDS** coordinates on RefSeq mRNAs (`XM_007698023.1`, `XM_007699256.1`, `XM_007696384.1`, `XM_007706020.1`) and reference **genomic scaffold** coordinates on `GCF_000338995.1` (`NW_006911896.1 (+)`, `NW_006911898.1 (-)`, `NW_006911895.1 (+)`, `NW_006911909.1 (+)`).
