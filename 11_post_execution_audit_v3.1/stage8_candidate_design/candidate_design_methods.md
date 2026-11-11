# Stage 8 — Independent Candidate-Region Design & Targetability Separation (BRD v3.1)

**Audited:** 2026-11-11 08:30:11 UTC

- **Coordinate & Sequence Audit:** All `50/50` candidate sequences match `candidate_sequences_v3_1.fasta`, exact length, and SHA-256. Of the `50` candidates, `49/50` map 1:1 to contiguous CDS substrings on `BRD_Target_Gene_Sequences.fasta` (`49/49` verified identically), while `1` candidate (`BRD_CsSTE11_CAND009`, `dsRNAmax`) is a multi-isolate stitched `SYNTHETIC_KMER_CONSTRUCT` (`exact_cds_substring=False`, `genomic_start_1based=-1`).
- **Redundancy Reduction:** The `50` overlapping windows cluster cleanly into `12` independent target regions (`>=25%` CDS overlap or `>=25%` shared 21-mers); `8` independent regions retain clean zero-exact-21mer candidates.
- **Strict Separation in `candidate_master_v3_1.tsv`:**
  1. **Sequence-derived metrics** (`seq_*`): GC%, homopolymer run length, dinucleotide entropy, DUST repeat score.
  2. **RNA structural predictions** (`struct_viennarna_*`): `ViennaRNA 2.8.0` whole-candidate self-fold MFE and unpaired accessibility probability (`RNAhybrid` reported transparently as `NOT_RUN_DEPENDENCY_BLOCKED`).
  3. **Heuristic targetability index** (`heuristic_dsrna_targetability_score`): Never interpreted as proof of experimental RNAi efficacy.
