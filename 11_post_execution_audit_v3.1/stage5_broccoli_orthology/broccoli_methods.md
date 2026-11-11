# Stage 5 — Broccoli v2.0.1 Orthology Methodology & Comparison (BRD v3.1)

**Executed:** 2026-11-11 08:16:19 UTC  
**Execution Route:** Route A — Native Windows 11 (`broccoli.exe v2.0.1` + `diamond.exe v2.2.8`)  
**Naming Guard:** These results are **Broccoli v2.0.1** (`-phylogenies nj`) orthogroups and are never labelled OrthoFinder orthogroups. Original `DIAMOND RBH` outputs in `11_post_execution_audit_v3/` are preserved completely untouched.

## Methodology
1. **Proteome validation (`broccoli_proteome_qc_v3_1.tsv`):** Split `all_proteomes_combined_v3.faa` cleanly into `3` single-species proteomes (`Bsor_ND90Pr` = `12,210`, `Bsor_Sacc._Shoemaker` = `10,755`, `Bcin_ASM83294v1` = `13,749`; `0` duplicate IDs, `0` empty sequences).
2. **Search strategy (`-combined_search` omitted):** Because proteomes span `10,755–13,749` proteins and pairwise search on `3` species takes seconds, default per-pair DIAMOND searches were executed cleanly to avoid combined-database e-value approximation assumptions.
3. **Tree inference (`-phylogenies nj`):** Used Broccoli v2's built-in Rust neighbor-joining tree builder because FastTree `MicrobesOnline` upstream downloads were gated by an anti-bot challenge.
4. **Comparison (`broccoli_vs_diamond_rbh_comparison_v3_1.tsv`):** `22915/26745` (`85.68%`) of Broccoli orthologous pairs co-cluster within the same DIAMOND RBH group. Concordance reflects shared sequence-similarity signal (not independent biological proof).
