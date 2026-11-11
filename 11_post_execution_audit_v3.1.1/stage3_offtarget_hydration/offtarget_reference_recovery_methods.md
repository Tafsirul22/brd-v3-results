# Off-Target Reference Hydration & Verification Methods (BRD v3.1.1)

**Verified:** 2026-11-11 09:21:09 UTC

## 1. Panel Coverage & Authoritative Checksums
`hydrate_offtarget_references_v3_1_1.py` hydrates all **18 references** (`17` non-target crops, beneficial fungi, and soil bacteria + `1` self *B. sorokiniana* `ND90Pr` control; `33` FASTA files: `18` transcriptomes + `15` genomes, excluding unscreened multi-GB wheat/barley/maize nuclear genomes) via official **NCBI Datasets v2 REST endpoints**.

## 2. Integrity Verification Protocol
For every retrieved organism archive:
1. Confirms `assembly_data_report.jsonl` accession (`GCF_*` / `GCA_*`) and NCBI `tax_id` match the documented baseline (`fail-fast` on mismatch).
2. Validates FASTA readability, non-zero nucleotide records, IUPAC nucleotide character set, exact sequence/scaffold count, exact byte size, and exact **SHA-256 checksum** against `11_post_execution_audit_v3.1/stage7_offtarget/offtarget_reference_metadata.tsv`.
3. Never silently substitutes newer assembly versions or accepts SHA-256 mismatches (`failed=0`, `validated_files=33/33`).

## 3. Screening Reproducibility Status
- **Clean-extraction hydration & checksum verification:** Verified across retrieved files (`failures=0`).
- **Combined DB rebuild & screening smoke test:** `BUILT_OK (smoke_representative_panel.fna); blastn_smoke_hits=102`; smoke test verified `BRD_CsSTE11_CAND013` has `0` exact 21-mer hits while `BRD_CsFUS3_CAND001` reproduces expected exact beneficial-fungi hits (`26` hits observed).
- **Historical baseline preservation:** Existing v3.1 screening outputs (`1,214` exact-hit rows, `16,073` near-hit rows) remain labelled as preserved v3.1 baseline results unless full panel BLAST+ screening is explicitly rerun.
