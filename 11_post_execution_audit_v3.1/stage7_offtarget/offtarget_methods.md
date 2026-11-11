# Stage 7 — Expanded Off-Target Screening & Mismatch Resolution Methodology (BRD v3.1)

**Executed:** 2026-11-11 08:23:44 UTC

## Scope, Thresholds & Explicit Caveats
1. **Taxonomic panel (`offtarget_reference_metadata.tsv`):** `17` non-target organisms (`4` crops: wheat, barley, rice, maize; `9` beneficial fungi: *T. atroviride*, *T. virens*, *T. reesei*, *T. harzianum*, *B. bassiana*, *M. anisopliae*, *M. robertsii*, *R. irregularis*, *R. clarus*; `4` rhizosphere/soil bacteria: *B. subtilis*, *B. velezensis*, *P. fluorescens*, *P. putida*) plus masked self *B. sorokiniana* `ND90Pr` control. Verified local FASTA SHA-256 checksums are recorded for every present transcript & genome file.
2. **Exact & near-hit resolution:**
   - **Exact 21-nt hard threshold preserved (`passes_zero_exact_21mer_filter`):** `32/50` pass; `18/50` excluded (including all `13/13` `CsFUS3` candidates). Alignment strand (`+`/`-`) and 1-based start offsets are resolved directly against reference sequences in `offtarget_exact_hits_v3_1.tsv`.
   - **Near matches (`offtarget_near_hits_v3_1.tsv`):** Enriches all `16,073` BLASTn-short hits (`20/21`, `19/21`, `18/21`) with full 21-nt subject window extraction, mismatch count (`total_21nt_mismatches`), and **exact 1-based mismatch positions (`1..21`)** inside the guide.
3. **Near-match burden visibility (`AUDIT_08`):** `29/50` candidates (`7/8` clean shortlist representatives!) carry `Elevated_20of21_Near_Match_Burden` (`3–14` unique non-target `20/21` guides across leading representatives). Passing exact 21-mer filtering NEVER equals biological safety.
4. **Coverage boundaries:** Whole nuclear genomes of wheat (`T. aestivum`), barley (`H. vulgare`), and maize (`Z. mays`) (`>20 GB`) were excluded locally (RefSeq transcriptomes screened), and unsequenced soil/phyllosphere biota remain uncovered.
