# Bipolaris-RNAi Designer (BRD) v3.1

**Updated:** 2026-11-11 08:30:13 UTC

## 1. Project Overview
Bipolaris-RNAi Designer v3.1 is a controlled reproducibility, tool-validation, and provenance upgrade of the computational RNAi candidate prioritisation pipeline (`11_post_execution_audit_v3.1/`) targeting *Bipolaris sorokiniana* (syn. *Cochliobolus sativus*).

## 2. Research Objectives
1. Replace undocumented Windows gaps with transparent GitHub tool investigations (`biopython`, `compleasm`, `Broccoli`) and empirical Route A/B/C execution logs.
2. Execute dedicated Biopython assembly QC across all 19 *B. sorokiniana* isolates with synthetic N50/L50 tests.
3. Evaluate Compleasm (`v0.2.10`) with verified fungal OrthoDB v12 lineages (`pleosporales_odb12`, `fungi_odb12`) natively on Windows without fabricating scores.
4. Build and run **Broccoli v2.0.1** natively on Windows 11 with DIAMOND `v2.2.8` and compare against preserved DIAMOND RBH results.
5. Fix clean-extraction packaging, SHA-256 manifest synchronization, input genome hydration (`GCF_000338995.1_Cocsa1_genomic.fna`), optional `ViennaRNA` guards, and clean-environment testing.

## 3. Data Sources
- **19 *B. sorokiniana* genome assemblies** (`authoritative_data/auth_isolates.tsv`, `genome_checksums_manifest_v3_1.tsv`).
- **4 target genes** (`CsSTE11`, `CsPKS1`, `CsSp1`, `CsFUS3`) verified against RefSeq/GenBank/UniProt/PHI-base v4.19/DEG/InterPro/GSE130517.
- **18 expanded off-target reference entries** (`17` non-target crops, beneficial fungi, soil bacteria + masked self `ND90Pr`, hashed in `stage7_offtarget/offtarget_reference_metadata.tsv`).

## 4. Target-Gene Selection
See `stage9_biological_targets/target_evidence_table_v3_1.tsv` (`CsSTE11`, `CsPKS1`, `CsSp1`, `CsFUS3`) preserving all PHI:2260, GSE130517 `3B/4B/5B`, and `CsFUS3` selectivity caveats.

## 5. Genome QC Methodology
See `stage3_assembly_qc/assembly_qc_methods.md` (`src/assembly_qc_biopython.py`). Strictly labelled Biopython QC (not QUAST).

## 6. Compleasm Methodology
See `stage4_compleasm/compleasm_methods.md` (`src/stage4_compleasm_runner.py`). CLI/downloads verified (`rc=0`); full genome execution (`run`) recorded as `NOT_RUN_DEPENDENCY_BLOCKED` on Windows without `miniprot`/`hmmsearch`, Docker, or WSL2.

## 7. Broccoli Methodology
See `stage5_broccoli_orthology/broccoli_methods.md` (`src/stage5_broccoli_orthology.py`). Compiled `broccoli.exe v2.0.1` via Rust `1.99.1` GNU toolchain and ran all 4 steps natively (`-phylogenies nj` + `diamond.exe v2.2.8`) across 36,714 proteins (`9,589` OGs, `26,745` orthologous pairs).

## 8. DIAMOND RBH Methodology
Original v3 `DIAMOND v2.2.8` RBH results remain preserved unmodified under `11_post_execution_audit_v3/stage2_genome_qc_comparative/orthology_diamond_rbh/` (`7,245` core *B. sorokiniana* OGs).

## 9. Isolate-Conservation Methodology
See `stage6_conservation/conservation_methods.md`. Tracks ALL-19 vs HQ-17 (`WAI2406`/`WAI2411` separated), minimum identity, coverage, gaps, ambiguous N bases, and 21-nt windows across `conservation_alignment_files/`.

## 10. Off-Target Methodology
See `stage7_offtarget/offtarget_methods.md`. Enforces 0 exact non-target/paralog 21-mer filter (`32/50` pass) while reporting exact 20/21, 19/21, and 18/21 mismatch positions and elevated near-match burdens.

## 11. Candidate-Design Methodology
See `stage8_candidate_design/candidate_design_methods.md`. Clusters 50 overlapping windows into 12 independent target regions (`8` clean passing) with explicit separation of sequence, ViennaRNA structure, and targetability heuristic scores.

## 12. Ranking Methodology
See `stage10_ranking/ranking_methods.md`. Separates Biological Relevance (`0.30`), Targetability (`0.35`), and Selectivity (`0.35`), explaining Pareto-first `independent_region_rank_v3_1` vs pure composite ranking.

## 13. Literature Benchmark
See `stage11_literature_benchmark/literature_benchmark_audit.md` (`25` curated studies + `8` exploratory *F. graminearum* constructs).

## 14. Software Dependencies
See `requirements.txt`, `environment.yml`, and `V3_1_Environment_Report.md`.

## 15. Installation Instructions
See `tool_investigation/installation_instructions.md` and `execution_environment/V3_1_Execution_Strategy.md`.

## 16. Execution Instructions
```bash
# 1. Hydrate missing primary reference genome (or all 19 with no flag) in clean extraction:
python 11_post_execution_audit_v3.1/src/hydrate_inputs_v3_1.py --primary-only
# 2. Validate environment & run Biopython QC / Broccoli / Stages 6-11:
python 11_post_execution_audit_v3.1/src/check_environment_v3_1.py
python 11_post_execution_audit_v3.1/src/assembly_qc_biopython.py
```

## 17. Test Results
See `V3_1_Test_Results.txt` and `V3_1_Test_Summary.md`.

## 18. Limitations
- Computational research only; exogenous dsRNA uptake, transcript silencing, crop protection, and environmental biosafety in *B. sorokiniana* are untested.
- Linux-only tools (`compleasm run`, `QUAST`, `BUSCO`, `OrthoFinder`, `RNAhybrid`, `FastTree bionj`) did not execute natively on Windows 11.
- Multi-GB nuclear genomes of wheat, barley, and maize remain outside local off-target coverage (RefSeq transcriptomes screened).

## 19. Reproducibility Instructions
Extract `11_post_execution_audit_v3.1.zip` into a clean directory, install `11_post_execution_audit_v3.1/requirements.txt`, optionally run `hydrate_inputs_v3_1.py --primary-only`, and run `pytest 11_post_execution_audit_v3.1/tests/test_brd_v3_1.py tests/test_brd.py -v`.

## 20. Citation and License Information
Third-party license notices for Biopython, Compleasm, BUSCO/OrthoDB datasets (`CC-BY-ND-4.0`), Broccoli (`GPL-3.0`), NCBI E-utilities, Europe PMC, and OpenAlex are archived in `.licenses/` and documented in `tool_investigation/tool_licenses.md`.
