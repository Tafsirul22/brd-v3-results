# Bipolaris-RNAi Designer v3 (`11_post_execution_audit_v3`)

Reproducible computational genomics, isolate-aware dsRNA design, expanded off-target selectivity screening, and literature benchmarking upgrade for *Bipolaris sorokiniana*.

## Key Scope Boundaries
- Purely computational research release; **no wet-lab validation, fungal dsRNA uptake assay, or crop protection assay in *B. sorokiniana* is claimed**.
- Windows native tools executed: Python 3.13.14, BLAST+ 2.17.0+, DIAMOND 2.2.8, NCBI Datasets 18.41.0, MAFFT v7.526, ViennaRNA 2.8.0.
- Linux-only tools (`QUAST`, `BUSCO`, `OrthoFinder`, `RNAhybrid`) were **not available** (no WSL2 Linux distro installed) and are recorded honestly as `NOT_RUN_DEPENDENCY_BLOCKED` in `stage2_genome_qc_comparative/unavailable_failed_or_incomplete_analyses_v3.tsv`.

## Quick Verification & Clean Execution
```bash
# Set working directory to repository root (parent of 11_post_execution_audit_v3)
pytest 11_post_execution_audit_v3/tests/test_brd_v3.py -v
```
All code inside `brd/` and `11_post_execution_audit_v3/src/` uses relative paths via `11_post_execution_audit_v3/src/paths_v3.py` (override root via `BRD_OUTPUT_DIR`).

## Principal Directory Layout
- `authoritative_data/` — Immutable Stage 1 target/isolate/candidate/guide tables & FASTAs.
- `stage1_reconciliation/` — Baseline audit and discrepancy reconciliation log.
- `stage2_genome_qc_comparative/` — Native QC, repaired DIAMOND RBH orthology, 50-marker core frame audit, and tool blocker log.
- `stage3_isolate_conservation/` — Isolate-aware multi-metric conservation & HQ-17 vs artifact separation.
- `stage4_off_target_expanded/` — Expanded 17-non-target + self panel manifest, coverage/gaps, dual-strand 21-mer/near/long hits.
- `stage5_redundancy_design/` — Independent region clustering & exhaustive 180-nt CDS sweep.
- `stage6_multidim_ranking/` — Strictly separated 3D scores, Pareto fronts, weight sensitivity, bootstrap stability, 50 Evidence Cards.
- `stage7_literature_benchmark/` — Open-access primary JATS XMLs, verified benchmark dataset (`literature_benchmark_dataset_v3.tsv`), and qualitative analysis.
- `V3_Validation_and_Execution_Report.md` (`.html`, `.pdf`) — Complete 18-section scientific validation report.
- `hashes_v3.tsv` — SHA256 verification manifest of deliverables.
