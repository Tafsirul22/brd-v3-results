# Stage 6 — Multi-Isolate Conservation Methodology (BRD v3.1)

**Executed:** 2026-11-11 08:22:41 UTC

## Audit Principles Enforced
1. **Per-gene and per-candidate separation:** Target-level MAFFT CDS alignments (`conservation_alignment_files/*_cds_aligned.fasta`) and candidate-level semi-global alignments (`candidate_by_isolate_alignments_v3.tsv`, `950` candidate-by-isolate pairs) are evaluated independently across all `19` *B. sorokiniana* isolates.
2. **Minimum identity & 21-mer window preservation:** Reporting enforces `min_nt_identity_hq17_pct`, `min_21mer_preservation_hq17_pct`, `lowest_identity_isolate_hq17`, and `lowest_identity_isolate_all19` alongside means so localized SNPs (e.g., `LK93` in `CsSp1`, dropping 21-mer preservation to `79.38%` in `BRD_CsSp1_CAND008` despite `98.89%` nucleotide identity) are never masked.
3. **Assembly polishing artifact separation:** `WAI2406` (`+137` excess internal stops in 50 core markers) and `WAI2411` (`+182` excess stops) are preserved in `ALL-19` summaries and cleanly separated from the `17` high-quality (`HQ-17`) assemblies.
4. **21-nt sliding window assessment:** Every possible 21-nt guide (`candidate_length_nt - 20`) within all `50` candidates (`9,301` guide windows total) is tracked across `ALL-19` and `HQ-17`.
