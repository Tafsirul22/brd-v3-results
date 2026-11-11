# Stage 3 — Biopython Assembly Quality Control Methodology (BRD v3.1)

**Executed:** 2026-11-11 08:15:12 UTC  
**Execution Route:** Route A — Native Windows 11 (`Biopython SeqIO`)  
**Tool Replacement Notice:** This module is a dedicated **Biopython** Python assembly statistics and FASTA integrity analyzer. **QUAST was NOT executed** (`Linux/WSL2` absent) and these outputs are never reported as QUAST results.

## Metrics & Integrity Checks Computed
- **Per-assembly size & count metrics:** `total_length_bp`, `contig_count`, strict length thresholds `contigs_gt_500bp`, `contigs_gt_1kb`, `contigs_gt_5kb`, `contigs_gt_10kb`, `largest_contig_bp`, `smallest_contig_bp`, `mean_contig_bp`, `median_contig_bp`.
- **Contiguity metrics:** Exact `N50` and `L50` computed via `compute_n50_l50()` (validated on synthetic unit test assemblies).
- **Base composition & ambiguity:** `gc_percent` (total span), `gc_percent_ungapped_acgt` (`G+C / (A+C+G+T)` excluding `N` gaps), `n_base_count`, `other_iupac_ambiguous_count`, `ambiguous_nucleotide_count`, `ambiguous_nucleotide_percent`.
- **Integrity alerts (`assembly_qc_warnings_v3_1.tsv`):** Detects missing files, empty FASTAs, duplicate sequence IDs, MD5-identical duplicate contig records, invalid non-IUPAC characters, assemblies `<25 Mb`, and scaffold gap fractions `>1.0%`.

## Summary Snapshot Across All 19 Isolates
- **Available & parsed:** `19/19`
- **Warnings raised:** `2` (`ELEVATED_AMBIGUOUS_N_PCT`)

## Explicit Scientific Limitations (Not QUAST)
Biopython assembly metrics measure only intrinsic FASTA composition and contiguity. They **do not** independently establish:
1. Reference-based structural misassemblies, translocations, or relocations.
2. Biological accuracy of local consensus bases (e.g., long-read homopolymer indels in `WAI2406` / `WAI2411`).
3. Chromosome-scale completeness or absence of chimeric scaffolding.
4. Complete gene space representation or annotation completeness.
