# Bipolaris-RNAi Designer (BRD) v3 — Validation and Execution Report

**Generated:** 2026-11-11 07:01:26 UTC
**Execution Mode:** Full 8-Stage Computational Research Release (`11_post_execution_audit_v3`)
**Mandatory Scientific Disclosure:** All scores, Pareto tiers, and candidate shortlists in this report are **computational prioritizations** across genomic/transcriptomic datasets. Nothing in this report claims demonstrated fungal dsRNA uptake, in vivo gene knockdown, crop protection against *Bipolaris sorokiniana*, or environmental biosafety.

---

## 1. Research objective and computational scope
Bipolaris-RNAi Designer v3 upgrades the prior v1 (`10_final_reports/`) and v2 (`11_post_execution_audit/`) computational RNAi candidate prioritisation pipeline targeting **4 *Bipolaris sorokiniana* (syn. *Cochliobolus sativus*) pathogenicity/biology genes** (`CsFUS3`, `CsSTE11`, `CsSp1`, `CsPKS1`) across **19 public *B. sorokiniana* genome assemblies**.

Key methodological upgrades in v3:
- Reconciles all v1/v2 numerical and filename discrepancies (`reconciliation_log_v3.tsv`) and locks **1 authoritative dataset** (`authoritative_data/`).
- Fixes the M04 stale `.dmnd` cache bug in `brd/integrations/diamond.py`, recovering real core-proteome orthology and distinguishing annotation gaps from gene loss.
- Separates unpolished long-read homopolymer frameshift artifacts (`WAI2406`, `WAI2411`) from genuine biological polymorphisms across all 50 dsRNA candidates.
- Expands non-target screening from 7 transcriptomes to **17 non-target organisms** (4 crops, 9 beneficial fungi with whole genomes, 4 soil/rhizosphere bacteria) plus self-*B. sorokiniana* paralogs/genome on **both strands**.
- Clusters the 50 overlapping M07 windows into **12 non-redundant independent target regions**.
- Enforces strict **3-dimensional separation** (`biological_target_relevance_score`, `dsrna_targetability_score`, `sequence_selectivity_score`), Pareto optimality, weight sensitivity, isolate bootstrapping, and primary-source literature benchmarking.

---

## 2. Data sources and accession provenance
### Authoritative Target Loci (`auth_targets.tsv`)
| target   | locus_tag         | refseq_protein_id   | refseq_mrna_id   | uniprot_accession   |   cds_length_nt | verification_status                                                                           |
|:---------|:------------------|:--------------------|:-----------------|:--------------------|----------------:|:----------------------------------------------------------------------------------------------|
| CsPKS1   | COCSADRAFT_157080 | XP_007696213.1      | XM_007698023.1   | M2RQ22              |            6462 | Confirmed (Primer + Full RefSeq/GenBank Verified)                                             |
| CsFUS3   | COCSADRAFT_24574  | XP_007697446.1      | XM_007699256.1   | M2TCC8              |            1059 | Confirmed (Primer + Full RefSeq/GenBank Verified)                                             |
| CsSTE11  | COCSADRAFT_105370 | XP_007694574.1      | XM_007696384.1   | M2TIY7              |            2880 | Confirmed (Primer + Full RefSeq/GenBank Verified)                                             |
| CsSp1    | COCSADRAFT_164106 | XP_007704210.1      | XM_007706020.1   | M2RXZ0              |             270 | Confirmed via 5/5 Primers + Reconstructed 270-nt RefSeq CDS (No standalone nuccore accession) |

### Reconciliation of Prior v1/v2 Discrepancies (`reconciliation_log_v3.tsv`)
| discrepancy_id                      | entity                                              | earlier_value                                                                                                                                                                                                                                                        | later_v2_value                                                                                                                                                                                                                                                                                                                                                                          | verified_root_cause                                                                                                                                                                                                                                                                                                                                                                        | v3_resolution                                                                                                                                             |
|:------------------------------------|:----------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------|
| REC_01_STE11_CAND012_COMPOSITE      | BRD_CsSTE11_CAND012                                 | v1=0.7832; audit(2).zip=0.8744                                                                                                                                                                                                                                       | 0.917                                                                                                                                                                                                                                                                                                                                                                                   | In audit(2).zip, CsSTE11 used M05 target_evidence_score = 0.6*0.90 + 0.4*min(1, 87.63/150) = 0.7737 (weight 0.20 -> composite 0.8744). In audit_v2.zip, M06 biological rubric raised CsSTE11 biological_target_score to 0.9866 and injected it directly into M09 candidate composite (+0.20*(0.9866 - 0.7737) = +0.0426 -> 0.9170).                                                        | RESOLVED: Separated gene-level biological relevance completely from candidate-level dsRNA targetability & selectivity in Stage 6 (no gene-score leakage). |
| REC_02_STE11_CAND013_COMPOSITE      | BRD_CsSTE11_CAND013                                 | v1=0.7321; audit(2).zip=0.8633                                                                                                                                                                                                                                       | 0.9059                                                                                                                                                                                                                                                                                                                                                                                  | Exact same +0.0426 shift as CAND012 (0.8633 -> 0.9059) when CsSTE11 target_evidence_score changed from 0.7737 to 0.9866 at weight 0.20.                                                                                                                                                                                                                                                    | RESOLVED: Same Stage 6 separation; additionally clustered with CAND012 into single 5' target region in Stage 5.                                           |
| REC_03_CSFUS3_BIO_SCORE             | CsFUS3 biological_target_score                      | audit(2).zip=0.9034                                                                                                                                                                                                                                                  | 0.9915                                                                                                                                                                                                                                                                                                                                                                                  | audit(2).zip used hardcoded bio_eval_meta=0.98 & included GSE130517 RPKM/150 (0.40*0.98 + 0.25*0.9760 + 0.20*0.6014 + 0.15*1.00 = 0.9034). audit_v2.zip replaced this with the 4-criterion expression-excluded rubric (0.35*1.0 + 0.20*1.0 + 0.25*0.9660 + 0.20*1.0 = 0.9915).                                                                                                             | VERIFIED & RETAINED AS TRANSPARENT RUBRIC DIMENSION (with distinct KO pathogenicity, PHI-base/DEG/InterPro, conservation & valid expression context).     |
| REC_04_CSSTE11_BIO_SCORE            | CsSTE11 biological_target_score                     | audit(2).zip=0.8837                                                                                                                                                                                                                                                  | 0.9866                                                                                                                                                                                                                                                                                                                                                                                  | audit(2).zip computed 0.40*0.96 + 0.25*0.9320 + 0.20*0.5842 + 0.15*1.00 = 0.8837 (penalised by 1-read 2B_III RPKM 87.63/150). audit_v2.zip computed 0.35*1.0 + 0.20*1.0 + 0.25*0.9462 + 0.20*1.0 = 0.9866.                                                                                                                                                                                 | VERIFIED & RESOLVED: Invalid 2B_III expression denominator removed; biological evaluation transparent.                                                    |
| REC_05_CSSP1_BIO_SCORE              | CsSp1 biological_target_score                       | audit(2).zip=0.8904                                                                                                                                                                                                                                                  | 0.6993                                                                                                                                                                                                                                                                                                                                                                                  | audit(2).zip gave 0.8904 (0.40*0.86 + 0.25*0.9956 + 0.20*1.00 + 0.15*0.65), inflated by expression_score=1.00 from 1-read 1438.44 RPKM AND even caused a stale rank-text contradiction (#2 in TSV vs #3 in summary prose). audit_v2.zip excluded invalid expression and scored comparative DB support (0/3 PHI-base/DEG/InterPro = 0.0) and standalone verification (2/4 = 0.5) -> 0.6993. | VERIFIED & RESOLVED cleanly.                                                                                                                              |
| REC_06_CSSP1_EXPRESSION_SPIKE       | CsSp1 GSE130517 RPKM                                | v1/infection_expression.tsv=1438.44 (functional_evidence_summary=1438.44)                                                                                                                                                                                            | valid 3B/4B/5B max=139.61 (sum raw reads=75); 3A_I RPKM=1438.44 from raw_reads=1                                                                                                                                                                                                                                                                                                        | Direct inspection of GSE130517 matrix confirmed Bipso1_10075 (322 nt) had RPKM_3A_I=1438.44 from exactly 1 raw read in shallow sample 3A_I. Valid colonized root samples 3B/4B/5B give max=139.61 RPKM (75 total raw reads).                                                                                                                                                               | RESOLVED: Full raw matrix provenance verified; RPKM/150 excluded permanently; valid 3B/4B/5B metrics retained strictly as qualitative context.            |
| REC_07_TEST_PATHS_AND_FILENAME      | test_audit_v2.py & archives                         | audit(3).zip test: ['OUT_DIR = Path(r"C:\\Users\\ASUS BANGLADESH\\Desktop\\New folder (2)\\Output")', 'near = pd.read_csv(AUDIT_DIR / "off_target_near_matches_v2.tsv", sep="\\t")'] | archive near=['11_post_execution_audit/off_target_near_guide_matches_v2.tsv'] | audit_v2.zip test: ['OUT_DIR = Path(r"C:\\Users\\ASUS BANGLADESH\\Desktop\\New folder (2)\\Output")', '# Note: off_target_near_guide_matches_v2.tsv is the exact filename written by run_post_execution_audit.py', 'near = pd.read_csv(AUDIT_DIR / "off_target_near_guide_matches_v2.tsv", sep="\\t")'] | archive near=['11_post_execution_audit/off_target_near_guide_matches_v2.tsv'] | audit(3).zip packaged test_audit_v2.py while test_11 still queried 'off_target_near_matches_v2.tsv', whereas run_post_execution_audit.py wrote 'off_target_near_guide_matches_v2.tsv' at root (plus BLAST intermediates under intermediate/). Also OUT_DIR used an absolute Windows drive path.                                                                                            | RESOLVED: Filename fixed (passed 13/13 in _v2.zip) and ALL v3 code/tests + brd/config.py refactored to relative/env-resolved portable paths.              |
| REC_08_M04_CORE_ACCESSORY_CACHE_BUG | 04_comparative_genomics/core_accessory & diamond.py | v1 core orthogroups=0/20288                                                                                                                                                                                                                                          | Bsor_Sacc._Shoemaker absent from all_proteomes_db.dmnd subjects (0/249,075 subject hits)                                                                                                                                                                                                                                                                                                | NEW v3 DISCOVERY: brd/integrations/diamond.py cached build_diamond_db on `if not dmnd.exists()`, silently skipping rebuild after Bsor_Sacc._Shoemaker was added to all_proteomes_combined.faa.                                                                                                                                                                                             | REPAIRED in Stage 2: Check SHA256/mtime before reusing .dmnd and rerun clean DIAMOND RBH across annotated proteomes.                                      |

### Pipeline Inputs & Software Versions (`software_and_reference_versions_v3.tsv`)
- **Isolate panel:** `19` *B. sorokiniana* assemblies (`GCF_000338995.1` RefSeq `ND90Pr` + 18 GenBank assemblies; `auth_isolates.tsv`).
- **Functional databases:** PHI-base v4.19, DEG essential genes, InterPro (`05_functional_evidence/`), NCBI GEO `GSE130517` barley-root RNA-seq.
| component                           | version_or_status                                                | notes                                                                         |
|:------------------------------------|:-----------------------------------------------------------------|:------------------------------------------------------------------------------|
| Python                              | 3.13.14                                                          | Windows-11-10.0.26300-SP0                                                     |
| biopython                           | 1.86                                                             | Sequence parsing & semi-global alignments                                     |
| pandas / numpy / scipy              | 2.3.3 / 2.4.0 / 1.16.3                                           | Statistical tables, Spearman, bootstrap                                       |
| ViennaRNA                           | 2.8.0                                                            | Duplexfold MFE & whole-candidate RNA folding                                  |
| matplotlib / networkx               | 3.10.8 / 3.6.1                                                   | Publication figures & RBH/overlap clustering                                  |
| BLAST+ (blastn/tblastn/makeblastdb) | blastn: 2.17.0+                                                  | Off-target blastn/blastn-short & core marker tblastn                          |
| DIAMOND                             | diamond version 2.2.8                                            | Repaired all-vs-all RBH orthology clustering                                  |
| NCBI datasets CLI                   | datasets version: 18.41.0                                        | Expanded reference genome & transcriptome retrieval                           |
| MAFFT                               | Active code page: 65001                                          | M04/M07 baseline multiple sequence alignments                                 |
| dsRNAmax                            | 2026/11/11 13:01:27 dsRNAmax - dsRNA maximizer (Version: 1.1.14) | Included in M07 baseline provenance                                           |
| QUAST / BUSCO / OrthoFinder         | NOT_RUN_DEPENDENCY_BLOCKED                                       | Linux/WSL2/Conda dependencies absent on Windows host; zero outputs fabricated |
| RNAhybrid                           | NOT_RUN_DEPENDENCY_BLOCKED                                       | Linux-only Bioconda binary; ViennaRNA 2.8.0 used as native comparator         |

---

## 3. Genome QC and assembly limitations
All 19 assemblies span `32.43–38.24 Mb` and `49.24–50.25%` GC (`assembly_metrics_19_isolates_v3.tsv`), confirming conspecificity and ruling out large-scale genome contamination. However, two structural assembly limitations materially impact downstream inference:
1. **Scaffold fragmentation & N-gaps:** `ND90Pr` (`154` scaffolds, `3.48%` N-bases = `~1.19 Mb`) and `Sacc._Shoemaker` (`GCA_013416765.1`, `3.48%` N-bases) contain extensive scaffolding gaps, partially truncating flanking or edge alignments.
2. **Unpolished long-read homopolymer indel frameshifts (`WAI2406`, `WAI2411`):** Testing 50 conserved single-copy core fungal markers via `tblastn` establishes an exon-intron boundary HSP stop baseline of `52–54` stops (`26–27` markers) across `ND90Pr` and all `16` remaining high-quality assemblies (`HQ-17`). By contrast, `WAI2406` (`190` stops, `+137` excess over `ND90Pr`) and `WAI2411` (`235` stops, `+182` excess) exhibit systematic frameshift disruptions across conserved housekeeping loci despite 50/50 marker detection (`isolates_artifact_vs_quality_v3.tsv`). Therefore, indels inside `WAI2406`/`WAI2411` cannot be treated as true biological gene loss.
3. **External Linux QC tool blockage:** `QUAST` and `BUSCO` did **not** run natively on Windows (no WSL2 Linux distro present) and are never substituted or faked (`Section 16`).

---

## 4. Comparative genomics results
1. **Stale `.dmnd` Cache Repair (`brd/integrations/diamond.py`):** M04 previously reported `0/20,288` core orthogroups shared with `Bsor_Sacc._Shoemaker` because `build_diamond_db` skipped rebuilding an outdated `all_proteomes_db.dmnd`. After repairing the DB rebuild and running clean `DIAMOND v2.2.8` all-vs-all RBH clustering across `Bsor_ND90Pr` (`12,210` proteins), `Bsor_Sacc._Shoemaker` (`10,755` proteins), and `Botrytis cinerea ASM83294v1`, **`7,245` core orthogroups** are cleanly recovered across both *B. sorokiniana* isolates.
2. **Annotation Gap vs Gene Loss (`annotation_comparability_audit_v3.tsv`):**
| metric                                                              |   count | interpretation                                                                                                                                                   |
|:--------------------------------------------------------------------|--------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ND90Pr total annotated proteins                                     |   12210 | NCBI RefSeq GCF_000338995.1 (JGI Cocsa1 model)                                                                                                                   |
| Sacc._Shoemaker total annotated proteins                            |   10755 | NCBI GenBank GCA_013416765.1 (1,455 fewer protein models than ND90Pr)                                                                                            |
| v1/v2 stale-cache core orthogroups (ND90Pr + Sacc.)                 |       0 | BUG: old all_proteomes_db.dmnd omitted Sacc._Shoemaker subjects                                                                                                  |
| v3 repaired Core orthogroups (ND90Pr & Sacc._Shoemaker present)     |    7245 | Genuine conserved core proteome shared across both annotated B. sorokiniana genomes                                                                              |
| Orthogroups annotated only in ND90Pr (Sacc._Shoemaker proteome = 0) |    2162 | Of 2385 such ND90Pr proteins, 2263 (94.9%) match Sacc._Shoemaker raw genome (>=70% aa tblastn) -> primarily unannotated gene models/gaps, not genuine gene loss! |
Of `2,385` `ND90Pr`-only proteins lacking a `Sacc._Shoemaker` proteome annotation, **`tblastn` recovers `2,263` (`94.9%`) directly in the raw `Sacc._Shoemaker` genome assembly (`>=70%` aa identity)**, proving the discrepancy reflects annotation incompleteness rather than biological deletion.
3. **OrthoFinder status:** Linux `OrthoFinder` (`mcl`/`fastme`) is blocked on Windows (`NOT_RUN_DEPENDENCY_BLOCKED`); v3 explicitly labels its orthology table as `DIAMOND_v2.2.8_AllVsAll_RBH`.

---

## 5. Conservation across isolates
Evaluating conservation using **minimum identity and exact 21-mer preservation across `HQ-17` isolates** separately from `WAI2406`/`WAI2411` prevents single-isolate indels or high means from masking divergence (`candidate_conservation_multimetric_v3.tsv`):
- **Total candidates assessed:** `50` across 4 targets; `50/50` achieve complete detection across the 19 assemblies (`2/19` isolates flagged `Unresolved_Assembly_Polishing_Artifact` for frameshift indels).
- **Mean vs Minimum disparity:** Several candidates maintain high mean identity (`>98%`) while dropping in `min_21mer_preservation_hq17_pct` or `min_21mer_preservation_all19_pct` due to single-nucleotide polymorphisms or `WAI2406/WAI2411` indels.

### Conservation Threshold Sensitivity (`conservation_threshold_sensitivity_v3.tsv`)
| scenario_id                             | evaluated_metric                  |   threshold_pct |   passing_candidates_n |   CsSTE11_pass |   CsPKS1_pass |   CsSp1_pass |   CsFUS3_pass |
|:----------------------------------------|:----------------------------------|----------------:|-----------------------:|---------------:|--------------:|-------------:|--------------:|
| Scenario_A1_Mean21mer_ge85pct           | mean_21mer_preservation_all19_pct |              85 |                     48 |             12 |            13 |           10 |            13 |
| Scenario_A2_Mean21mer_ge90pct           | mean_21mer_preservation_all19_pct |              90 |                     46 |             10 |            13 |           10 |            13 |
| Scenario_A3_Mean21mer_ge95pct           | mean_21mer_preservation_all19_pct |              95 |                     44 |              8 |            13 |           10 |            13 |
| Scenario_B1_MinNtID_All19_ge85pct       | min_nt_identity_all19_pct         |              85 |                     50 |             14 |            13 |           10 |            13 |
| Scenario_B2_MinNtID_All19_ge90pct       | min_nt_identity_all19_pct         |              90 |                     50 |             14 |            13 |           10 |            13 |
| Scenario_B3_MinNtID_All19_ge95pct       | min_nt_identity_all19_pct         |              95 |                     50 |             14 |            13 |           10 |            13 |
| Scenario_C1_MinNtID_HQ17_ge90pct        | min_nt_identity_hq17_pct          |              90 |                     50 |             14 |            13 |           10 |            13 |
| Scenario_C2_MinNtID_HQ17_ge95pct        | min_nt_identity_hq17_pct          |              95 |                     50 |             14 |            13 |           10 |            13 |
| Scenario_C3_Min21mer_HQ17_ge95pct       | min_21mer_preservation_hq17_pct   |              95 |                     29 |              6 |            12 |            0 |            11 |
| Scenario_D1_90pctOfIsolates_ge95pctNtID | isolates_ge95pct_nt_id_pct        |              90 |                     50 |             14 |            13 |           10 |            13 |

---

## 6. Off-target screening coverage and findings
The v3 screen expands from M08's 7 transcriptomes to **18 reference entries** (`17` non-target organisms + `1` self *B. sorokiniana* `ND90Pr` control with target loci masked; `expanded_reference_manifest_v3.tsv`, `reference_panel_coverage_and_gaps_v3.tsv`):
- **Crops (`4`):** Wheat (`T. aestivum`), Barley (`H. vulgare`), Rice (`O. sativa` `GCF_034140825.1` tx+genome), Maize (`Z. mays` `GCF_902167145.1` tx).
- **Beneficial Fungi (`9`):** *Trichoderma atroviride*, *T. virens*, *T. reesei*, *T. harzianum* (`GCF_003025095.1`), *Beauveria bassiana*, *Metarhizium anisopliae* (`GCF_039654215.1`), *M. robertsii* (`GCA_000591435.1`), *Rhizophagus irregularis*, *R. clarus* (`GCA_003203555.1`), all including **whole-genome dual-strand 21-mer scans**.
- **Rhizosphere/Soil Bacteria (`4`):** *Bacillus subtilis* (`GCF_000009045.1`), *B. velezensis* (`GCF_000508125.2`), *Pseudomonas fluorescens* (`GCF_900215245.1`), *P. putida* (`GCF_000412675.1`).
- **Passing zero-exact-21mer hard filter:** **`32/50`** (`18/50` excluded by `>=1` exact dual-strand 21-mer match across the expanded panel!).

### Expanded Selectivity Classification Summary (`candidate_selectivity_expanded_v3.tsv`)
| selectivity_classification_v3          |   candidates_n |
|:---------------------------------------|---------------:|
| Elevated_20of21_Near_Match_Burden      |             29 |
| Exact_21mer_NonTarget_Or_Paralog_Match |             18 |
| Low_20of21_Near_Match_Signal           |              3 |

---

## 7. Candidate redundancy and independent target regions
Many of the original 50 M07 candidates were overlapping sliding/multi-tool windows (`SiFi21`, `dsRNAbuilder`, `dsRNAmax`) on identical CDS loci. Clustering windows sharing `>=25%` CDS overlap or `>=25%` 21-mers reduces **50 candidates into `12` non-redundant independent target regions** (`independent_target_regions_v3.tsv`):
| region_id     | target   |   cds_span_start_nt |   cds_span_end_nt |   region_span_nt |   n_redundant_candidates |   n_passing_clean_hq85 | design_tools_present                              |
|:--------------|:---------|--------------------:|------------------:|-----------------:|-------------------------:|-----------------------:|:--------------------------------------------------|
| CsFUS3_REG01  | CsFUS3   |                   1 |               300 |              300 |                        5 |                      0 | SIREN;dsRNAmax                                    |
| CsFUS3_REG02  | CsFUS3   |                 326 |               710 |              385 |                        7 |                      0 | BRD_IsolateAware;LUCID_sifi21_plus                |
| CsFUS3_REG03  | CsFUS3   |                 676 |               875 |              200 |                        1 |                      0 | LUCID_sifi21_plus                                 |
| CsPKS1_REG01  | CsPKS1   |                   1 |               300 |              300 |                        5 |                      5 | SIREN;dsRNAmax                                    |
| CsPKS1_REG02  | CsPKS1   |                4526 |              4735 |              210 |                        1 |                      1 | BRD_IsolateAware                                  |
| CsPKS1_REG03  | CsPKS1   |                5651 |              5830 |              180 |                        1 |                      1 | BRD_IsolateAware                                  |
| CsPKS1_REG04  | CsPKS1   |                5951 |              6210 |              260 |                        6 |                      6 | BRD_IsolateAware;LUCID_sifi21_plus                |
| CsSTE11_REG01 | CsSTE11  |                   1 |               300 |              300 |                        8 |                      6 | BRD_IsolateAware;LUCID_sifi21_plus;SIREN;dsRNAmax |
| CsSTE11_REG02 | CsSTE11  |                1076 |              1285 |              210 |                        2 |                      2 | BRD_IsolateAware                                  |
| CsSTE11_REG03 | CsSTE11  |                1401 |              1610 |              210 |                        1 |                      1 | BRD_IsolateAware                                  |
| CsSTE11_REG04 | CsSTE11  |                2226 |              2525 |              300 |                        3 |                      0 | LUCID_sifi21_plus                                 |
| CsSp1_REG01   | CsSp1    |                   1 |               255 |              255 |                       10 |                     10 | BRD_IsolateAware;LUCID_sifi21_plus;SIREN          |

- **Exhaustive 180-nt Full-CDS Sweep (`exhaustive_cds_180nt_sweep_v3.tsv`):** Scanning `25-nt` steps across all 4 ORFs shows only **`5`** clean 180-nt windows on `CsFUS3` (`1,059 nt`) escape exact 21-mer homology with beneficial fungi, proving that MAPK kinase conservation spans virtually the entire `CsFUS3` gene (rather than being an M07 sampling quirk). Outside M07 windows, `112` clean high-conservation exploratory windows (`HQ-17 >= 95%`, `0` exact 21-mers across non-targets) were identified on `CsSTE11`/`CsPKS1`.

---

## 8. Biological target evidence
Biological target relevance (`biological_target_relevance_score`) is scored strictly at the gene level with context-aware PHI-base evaluation and valid `3B/4B/5B` RNA-seq reporting excluded from composite equations (`biological_target_relevance_v3.tsv`):
| target   |   biological_target_relevance_score |   subscore_direct_ko_pathogenicity |   subscore_phi_context_aware |   subscore_sequence_traceability |   gse130517_valid_3b4b5b_mean_rpkm_separate | phi_base_context_interpretation                                                                                                                                     |
|:---------|------------------------------------:|-----------------------------------:|-----------------------------:|---------------------------------:|--------------------------------------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CsFUS3   |                              0.9921 |                               1    |                         1    |                              1   |                                       62.53 | Direct/close PHI loss_of_pathogenicity (PHI:142__PHI:2366__PHI:2398__PHI:125577); 10/10 plant-fungal homologs show virulence reduction                              |
| CsSTE11  |                              0.9921 |                               1    |                         1    |                              1   |                                       22.12 | Direct/close PHI loss_of_pathogenicity (PHI:4778, PHI:2397__PHI:125579); 9/10 plant-fungal homologs show virulence reduction; note M. robertsii ortholog unaffected |
| CsPKS1   |                              0.5796 |                               0.25 |                         0.25 |                              1   |                                       13.19 | Direct B. sorokiniana entry PHI:2260 records 'unaffected pathogenicity' (secondary DHN melanin phenotype; mixed toxin roles only in distant homologs)               |
| CsSp1    |                              0.55   |                               1    |                         0    |                              0.5 |                                       79.3  | Absent from PHI-base v4.19 (effect relies on Zhang et al. 2022 primary paper, not PHI curation)                                                                     |

- Crucial biological context: **`CsPKS1`** ranks **`4th`** (`0.5796`) in target relevance because its direct *B. sorokiniana* PHI-base entry `PHI:2260` (`Pks1`) explicitly documents **`unaffected pathogenicity`** (affecting DHN spore melanization, not host virulence!).

---

## 9. Candidate targetability
Candidate targetability (`dsrna_targetability_score`) isolates sequence-intrinsic properties (`min_nt_identity_hq17_pct`, `min_21mer_preservation_hq17_pct`, `guides_pct_preserved_hq17`, `efficient_sifi21_fraction`, `guide_acc_mean`, and sequence entropy/complexity) with zero biological-target leakage and zero tool-agreement weighting.

---

## 10. Selectivity assessment
All off-target counts (`exact_21_crops_total`, `exact_21_beneficial_fungi_tx`, `exact_21_beneficial_fungi_gen`, bacteria, self-paralog/genome, `20/21` near matches, `>=30 nt` alignments) feed `sequence_selectivity_score` (`heuristic_sequence_selectivity_index`). Zero exact 21-mer hits across covered panels does **not** prove environmental safety because wheat, barley, and maize whole non-coding genomes (~21 GB combined) and unsequenced soil/phyllosphere biota remain outside local coverage (`reference_panel_coverage_and_gaps_v3.tsv`). Conversely, insect 21-nt feeding comparators (e.g., *Diabrotica* Bachman et al. 2013) are explicitly distinguished as cross-kingdom comparators rather than fungal silencing thresholds.

---

## 11. Pareto and sensitivity analysis
- **3D Pareto Optimal Tiers:** Evaluated both across all 50 candidates (`pareto_front_3d_all50`) and across the **`12` independent target region representatives** (`pareto_front_3d_regions`). Among passing independent regions (`n=8`), Pareto Tier-1 captures candidates trading off higher `CsSTE11`/`CsSp1` biological relevance against higher `CsPKS1`/`CsSTE11` sequence selectivity.
- **Weight Sensitivity (`stage6_weight_sensitivity_spearman_v3.tsv`):**
| scenario_a                        | scenario_b                        |   spearman_rho |   p_value |
|:----------------------------------|:----------------------------------|---------------:|----------:|
| S1_Balanced_30_35_35              | S2_SelectivityFirst_20_30_50      |         1      |  0        |
| S1_Balanced_30_35_35              | S3_BioPathogenicityFirst_50_25_25 |         1      |  0        |
| S1_Balanced_30_35_35              | S4_TargetabilityFirst_20_50_30    |         0.9286 |  0.000863 |
| S1_Balanced_30_35_35              | S5_EqualWeights_33_33_33          |         1      |  0        |
| S2_SelectivityFirst_20_30_50      | S3_BioPathogenicityFirst_50_25_25 |         1      |  0        |
| S2_SelectivityFirst_20_30_50      | S4_TargetabilityFirst_20_50_30    |         0.9286 |  0.000863 |
| S2_SelectivityFirst_20_30_50      | S5_EqualWeights_33_33_33          |         1      |  0        |
| S3_BioPathogenicityFirst_50_25_25 | S4_TargetabilityFirst_20_50_30    |         0.9286 |  0.000863 |
| S3_BioPathogenicityFirst_50_25_25 | S5_EqualWeights_33_33_33          |         1      |  0        |
| S4_TargetabilityFirst_20_50_30    | S5_EqualWeights_33_33_33          |         0.9286 |  0.000863 |
- **1,000-Iteration HQ-17 Isolate Bootstrap (`stage6_bootstrap_stability_v3.tsv`):**
| candidate_id        | target   | region_id     |   observed_clean_region_rank |   bootstrap_mean_rank |   bootstrap_median_rank |   bootstrap_rank_sd |   bootstrap_rank_ci95_low |   bootstrap_rank_ci95_high |   prob_top1_pct |   prob_top3_pct |   prob_top5_pct |
|:--------------------|:---------|:--------------|-----------------------------:|----------------------:|------------------------:|--------------------:|--------------------------:|---------------------------:|----------------:|----------------:|----------------:|
| BRD_CsSTE11_CAND013 | CsSTE11  | CsSTE11_REG01 |                            1 |                     1 |                       1 |                   0 |                         1 |                          1 |             100 |             100 |             100 |
| BRD_CsSTE11_CAND014 | CsSTE11  | CsSTE11_REG03 |                            2 |                     2 |                       2 |                   0 |                         2 |                          2 |               0 |             100 |             100 |
| BRD_CsSTE11_CAND010 | CsSTE11  | CsSTE11_REG02 |                            3 |                     3 |                       3 |                   0 |                         3 |                          3 |               0 |             100 |             100 |
| BRD_CsPKS1_CAND012  | CsPKS1   | CsPKS1_REG02  |                            4 |                     4 |                       4 |                   0 |                         4 |                          4 |               0 |               0 |             100 |
| BRD_CsPKS1_CAND010  | CsPKS1   | CsPKS1_REG04  |                            5 |                     7 |                       7 |                   0 |                         7 |                          7 |               0 |               0 |               0 |
| BRD_CsPKS1_CAND013  | CsPKS1   | CsPKS1_REG03  |                            6 |                     5 |                       5 |                   0 |                         5 |                          5 |               0 |               0 |             100 |
| BRD_CsPKS1_CAND006  | CsPKS1   | CsPKS1_REG01  |                            7 |                     6 |                       6 |                   0 |                         6 |                          6 |               0 |               0 |               0 |
| BRD_CsSp1_CAND008   | CsSp1    | CsSp1_REG01   |                            8 |                     8 |                       8 |                   0 |                         8 |                          8 |               0 |               0 |               0 |

---

## 12. Literature benchmarking
| benchmark_category                 |   study_records |
|:-----------------------------------|----------------:|
| 3_Gene_Level_Only                  |              14 |
| 5_Ambiguous_or_Conflicting_Outcome |               5 |
| 4_Computational_Only               |               4 |
| 2_Reconstructed_Sequence           |               2 |

Across `25` curated primary records spanning *Fusarium graminearum*, *Botrytis cinerea*, *Sclerotinia sclerotiorum*, *Phakopsora pachyrhizi*, *Magnaporthe oryzae*, *Alternaria alternata*, *Aspergillus flavus*, *Rhizoctonia solani*, and *Bipolaris sorokiniana* (`literature_benchmark_dataset_v3.tsv`), primary full-text XML analysis shows:
1. **Delivery & uptake dependency (`Qiao et al. 2021, PMC8428832`; `Wytinck et al. 2020, PMC7391711`):** High dsRNA uptake species (*B. cinerea*, *S. sclerotiorum*, *A. niger*, *R. solani*, *V. dahliae*) contrast sharply with zero/weak uptake species (*Colletotrichum gloeosporioides*, *Trichoderma virens*, *P. infestans*). *B. sorokiniana* possesses intracellular RNAi machinery (`Leng et al. 2011`), yet **exogenous dsRNA uptake has not been experimentally demonstrated in *B. sorokiniana***.
2. **Construct length vs SiFi/tool length (`Werner et al. 2020, PMC7202221`):** Manual long constructs (`658–912 nt`) outperformed short tool-designed constructs (`173–193 nt`) on *F. graminearum* `AGO1/2` & `DCL1/2`, including failing (`No silencing`) in two short tool constructs (`AGO1-DCL2`, `AGO1-AGO2`).
3. **Dataset scope limit:** Because primary datasets span heterogeneous delivery routes (SIGS vs HIGS vs liquid culture), variable concentrations, and limited verified full construct sequences (`stage7_benchmark_report_v3.md`), Stage 7 deliberately reports a **structured qualitative + exploratory benchmark** rather than claiming predictive statistical accuracy.

---

## 13. Final prioritized candidate shortlist
The redundancy-reduced computational shortlist selects **1 highest-ranking representative per independent Stage 5 target region** (`8` clean representatives passing the zero-exact-21mer hard filter out of `12` independent regions; full 50-candidate Evidence Cards in `stage6_multidim_ranking/candidate_evidence_cards_v3.md`):

|   independent_region_rank_v3 | candidate_id        | target   | region_id     | cds_coords   |   candidate_length_nt |   biological_target_relevance_score |   dsrna_targetability_score |   sequence_selectivity_score |   heuristic_composite_prioritization_index |   pareto_front_3d_regions |   min_nt_identity_hq17_pct |   min_21mer_preservation_hq17_pct | selectivity_classification_v3     |
|-----------------------------:|:--------------------|:---------|:--------------|:-------------|----------------------:|------------------------------------:|----------------------------:|-----------------------------:|-------------------------------------------:|--------------------------:|---------------------------:|----------------------------------:|:----------------------------------|
|                            1 | BRD_CsSTE11_CAND013 | CsSTE11  | CsSTE11_REG01 | 1..210       |                   210 |                              0.9921 |                      0.8855 |                        0.955 |                                     0.9418 |                         1 |                     100    |                            100    | Low_20of21_Near_Match_Signal      |
|                            2 | BRD_CsSTE11_CAND014 | CsSTE11  | CsSTE11_REG03 | 1401..1610   |                   210 |                              0.9921 |                      0.8859 |                        0.94  |                                     0.9367 |                         1 |                     100    |                            100    | Low_20of21_Near_Match_Signal      |
|                            3 | BRD_CsSTE11_CAND010 | CsSTE11  | CsSTE11_REG02 | 1101..1280   |                   180 |                              0.9921 |                      0.8913 |                        0.79  |                                     0.8861 |                         1 |                     100    |                            100    | Elevated_20of21_Near_Match_Burden |
|                            4 | BRD_CsPKS1_CAND012  | CsPKS1   | CsPKS1_REG02  | 4526..4735   |                   210 |                              0.5796 |                      0.8879 |                        0.895 |                                     0.7979 |                         1 |                     100    |                            100    | Elevated_20of21_Near_Match_Burden |
|                            5 | BRD_CsPKS1_CAND010  | CsPKS1   | CsPKS1_REG04  | 6001..6180   |                   180 |                              0.5796 |                      0.9009 |                        0.82  |                                     0.7762 |                         1 |                     100    |                            100    | Elevated_20of21_Near_Match_Burden |
|                            6 | BRD_CsPKS1_CAND013  | CsPKS1   | CsPKS1_REG03  | 5651..5830   |                   180 |                              0.5796 |                      0.8335 |                        0.91  |                                     0.7841 |                         2 |                      99.44 |                             86.88 | Elevated_20of21_Near_Match_Burden |
|                            7 | BRD_CsPKS1_CAND006  | CsPKS1   | CsPKS1_REG01  | 5..254       |                   250 |                              0.5796 |                      0.8772 |                        0.865 |                                     0.7836 |                         2 |                     100    |                            100    | Elevated_20of21_Near_Match_Burden |
|                            8 | BRD_CsSp1_CAND008   | CsSp1    | CsSp1_REG01   | 1..180       |                   180 |                              0.55   |                      0.7524 |                        0.88  |                                     0.7363 |                         3 |                      98.89 |                             79.38 | Elevated_20of21_Near_Match_Burden |

---

## 14. Reasons each candidate ranks where it does
All `12` independent target region representatives (`8` passing + `4` excluded for exact 21-mer matches):
- **Rank 1 — `BRD_CsSTE11_CAND013` (`CsSTE11`, `CsSTE11_REG01`, n=8 overlapping windows, CLEAN_0_EXACT_21MER):** Bio=0.9921, Targ=0.8855 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.955 (crops 21-mer=0, ben.fungi tx/gen=0/0, 20/21=3, long=0). Heuristic index=`0.9418` (3D Region Pareto Tier `1`).
- **Rank 2 — `BRD_CsSTE11_CAND014` (`CsSTE11`, `CsSTE11_REG03`, n=1 overlapping windows, CLEAN_0_EXACT_21MER):** Bio=0.9921, Targ=0.8859 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.94 (crops 21-mer=0, ben.fungi tx/gen=0/0, 20/21=4, long=0). Heuristic index=`0.9367` (3D Region Pareto Tier `1`).
- **Rank 3 — `BRD_CsSTE11_CAND010` (`CsSTE11`, `CsSTE11_REG02`, n=2 overlapping windows, CLEAN_0_EXACT_21MER):** Bio=0.9921, Targ=0.8913 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.79 (crops 21-mer=0, ben.fungi tx/gen=0/0, 20/21=14, long=0). Heuristic index=`0.8861` (3D Region Pareto Tier `1`).
- **Rank 4 — `BRD_CsPKS1_CAND012` (`CsPKS1`, `CsPKS1_REG02`, n=1 overlapping windows, CLEAN_0_EXACT_21MER):** Bio=0.5796, Targ=0.8879 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.895 (crops 21-mer=0, ben.fungi tx/gen=0/0, 20/21=7, long=0). Heuristic index=`0.7979` (3D Region Pareto Tier `1`).
- **Rank 5 — `BRD_CsPKS1_CAND010` (`CsPKS1`, `CsPKS1_REG04`, n=6 overlapping windows, CLEAN_0_EXACT_21MER):** Bio=0.5796, Targ=0.9009 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.82 (crops 21-mer=0, ben.fungi tx/gen=0/0, 20/21=12, long=0). Heuristic index=`0.7762` (3D Region Pareto Tier `1`).
- **Rank 6 — `BRD_CsPKS1_CAND013` (`CsPKS1`, `CsPKS1_REG03`, n=1 overlapping windows, CLEAN_0_EXACT_21MER):** Bio=0.5796, Targ=0.8335 (HQ17 min nt=99.44%, min 21-mer=86.88%), Sel=0.91 (crops 21-mer=0, ben.fungi tx/gen=0/0, 20/21=6, long=0). Heuristic index=`0.7841` (3D Region Pareto Tier `2`).
- **Rank 7 — `BRD_CsPKS1_CAND006` (`CsPKS1`, `CsPKS1_REG01`, n=5 overlapping windows, CLEAN_0_EXACT_21MER):** Bio=0.5796, Targ=0.8772 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.865 (crops 21-mer=0, ben.fungi tx/gen=0/0, 20/21=9, long=0). Heuristic index=`0.7836` (3D Region Pareto Tier `2`).
- **Rank 8 — `BRD_CsSp1_CAND008` (`CsSp1`, `CsSp1_REG01`, n=10 overlapping windows, CLEAN_0_EXACT_21MER):** Bio=0.55, Targ=0.7524 (HQ17 min nt=98.89%, min 21-mer=79.38%), Sel=0.88 (crops 21-mer=0, ben.fungi tx/gen=0/0, 20/21=8, long=0). Heuristic index=`0.7363` (3D Region Pareto Tier `3`).
- **Rank 9 — `BRD_CsFUS3_CAND001` (`CsFUS3`, `CsFUS3_REG02`, n=7 overlapping windows, EXCLUDED_EXACT_21MER):** Bio=0.9921, Targ=0.8958 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.1 (crops 21-mer=0, ben.fungi tx/gen=4/4, 20/21=35, long=14). Heuristic index=`0.6462` (3D Region Pareto Tier `1`).
- **Rank 10 — `BRD_CsFUS3_CAND003` (`CsFUS3`, `CsFUS3_REG03`, n=1 overlapping windows, EXCLUDED_EXACT_21MER):** Bio=0.9921, Targ=0.8836 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.09 (crops 21-mer=1, ben.fungi tx/gen=0/0, 20/21=14, long=4). Heuristic index=`0.6384` (3D Region Pareto Tier `2`).
- **Rank 11 — `BRD_CsFUS3_CAND009` (`CsFUS3`, `CsFUS3_REG01`, n=5 overlapping windows, EXCLUDED_EXACT_21MER):** Bio=0.9921, Targ=0.858 (HQ17 min nt=100.0%, min 21-mer=100.0%), Sel=0.1 (crops 21-mer=0, ben.fungi tx/gen=5/5, 20/21=66, long=11). Heuristic index=`0.6329` (3D Region Pareto Tier `2`).
- **Rank 12 — `BRD_CsSTE11_CAND004` (`CsSTE11`, `CsSTE11_REG04`, n=3 overlapping windows, EXCLUDED_EXACT_21MER):** Bio=0.9921, Targ=0.7789 (HQ17 min nt=99.0%, min 21-mer=78.89%), Sel=0.1 (crops 21-mer=0, ben.fungi tx/gen=2/2, 20/21=24, long=7). Heuristic index=`0.6052` (3D Region Pareto Tier `3`).

---

## 15. Candidate-specific weaknesses and uncertainty
- **`BRD_CsSTE11_CAND013` (`CsSTE11_REG01`):** 3 unique 20/21-mer near-matches
- **`BRD_CsSTE11_CAND014` (`CsSTE11_REG03`):** 4 unique 20/21-mer near-matches
- **`BRD_CsSTE11_CAND010` (`CsSTE11_REG02`):** 14 unique 20/21-mer near-matches
- **`BRD_CsPKS1_CAND012` (`CsPKS1_REG02`):** 7 unique 20/21-mer near-matches; Target biology: B. sorokiniana PHI:2260 Pks1 is 'unaffected pathogenicity' (melanin only)
- **`BRD_CsPKS1_CAND010` (`CsPKS1_REG04`):** 12 unique 20/21-mer near-matches; Target biology: B. sorokiniana PHI:2260 Pks1 is 'unaffected pathogenicity' (melanin only)
- **`BRD_CsPKS1_CAND013` (`CsPKS1_REG03`):** 6 unique 20/21-mer near-matches; Target biology: B. sorokiniana PHI:2260 Pks1 is 'unaffected pathogenicity' (melanin only)
- **`BRD_CsPKS1_CAND006` (`CsPKS1_REG01`):** 9 unique 20/21-mer near-matches; Target biology: B. sorokiniana PHI:2260 Pks1 is 'unaffected pathogenicity' (melanin only)
- **`BRD_CsSp1_CAND008` (`CsSp1_REG01`):** 8 unique 20/21-mer near-matches; HQ-17 min nt identity dips to 98.89%; Target traceability/length: 270-nt short CDS without standalone nuccore accession

---

## 16. Failed, unavailable, or incomplete analyses
Honest inventory (`unavailable_failed_or_incomplete_analyses_v3.tsv`):
| tool        | status                     | blocker_or_alternative                                                                                    | score_fabrication_guard                                |
|:------------|:---------------------------|:----------------------------------------------------------------------------------------------------------|:-------------------------------------------------------|
| QUAST       | NOT_RUN_DEPENDENCY_BLOCKED | Requires Linux/WSL2 + Perl/minimap2/GlimmerHMM; Windows host has no WSL2 distro.                          | Zero BUSCO/QUAST/OrthoFinder synthetic outputs claimed |
| BUSCO       | NOT_RUN_DEPENDENCY_BLOCKED | Requires Linux/macOS HMMER/MetaEuk/Augustus + fungi_odb10/pleosporales_odb10.                             | Zero BUSCO/QUAST/OrthoFinder synthetic outputs claimed |
| OrthoFinder | NOT_RUN_DEPENDENCY_BLOCKED | Requires Linux MCL/FastME/OrthoFinder runtime; mcl & fastme absent natively on Windows.                   | Zero BUSCO/QUAST/OrthoFinder synthetic outputs claimed |
| RNAhybrid   | NOT_RUN_DEPENDENCY_BLOCKED | Linux/macOS Bioconda binary unavailable natively on Windows; ViennaRNA 2.8.0 RNA.duplexfold used cleanly. | Zero BUSCO/QUAST/OrthoFinder synthetic outputs claimed |

Additional incomplete / missing elements documented without substitution:
- `WAI2406` (`GCA_008452705.1`) and `WAI2411` (`GCA_008452715.1`) contain unpolished homopolymer frameshift artifacts (`+137` / `+182` stops above baseline), blocking true coding polymorphism resolution in those 2 assemblies.
- Wheat (`T. aestivum`), barley (`H. vulgare`), and maize (`Z. mays`) whole non-coding nuclear genomes (`>20 GB`) were excluded from local storage; off-target screening covers their complete RefSeq transcriptomes.
- Unsequenced rhizosphere/phyllosphere microorganisms and wet-lab dsRNA uptake in *B. sorokiniana* remain empirically untested.

---

## 17. Reproducibility instructions
All v3 modules resolve relative to repository root (`11_post_execution_audit_v3/src/paths_v3.py`, optional `BRD_OUTPUT_DIR` env override) and contain zero hardcoded machine paths:
```bash
# 1. Activate environment containing packages listed in 11_post_execution_audit_v3/requirements.txt
# 2. Verify Stages 1-8 end-to-end & execute full test suite:
pytest 11_post_execution_audit_v3/tests/test_brd_v3.py tests/test_brd.py 11_post_execution_audit/tests/test_audit_v2.py -v
# 3. Verify SHA256 integrity manifest:
# Check 11_post_execution_audit_v3/hashes_v3.tsv
```
(See `11_post_execution_audit_v3/README.md` for optional Linux/WSL2 installation of `QUAST`, `BUSCO`, and `OrthoFinder`.)

---

## 18. Explicit distinction between computational predictions and experimentally demonstrated effects
| Layer | Status in *Bipolaris sorokiniana* | Source / Basis |
|---|---|---|
| **Intracellular fungal RNAi pathway** | Demonstrated (`CsPKS1` / `GFP` / `ToxA` hairpin RNAi) | `Leng et al. (2011)` (`doi:10.1111/j.1364-3703.2010.00666.x`) |
| **Target gene virulence phenotypes (`CsFUS3`, `CsSTE11`, `CsSp1`, `CsPKS1`)** | Demonstrated via gene knockouts / hairpin RNAi | `Leng & Zhong (2015) PMC4443973`; `Zhang et al. (2022) PMC8743017`; `Leng et al. (2011) PHI:2260` |
| **Exogenous dsRNA uptake & topical SIGS protection in *B. sorokiniana*** | **NOT DEMONSTRATED (untested computational hypothesis)** | Extrapolated cautiously from `Wang et al. (2016)`, `Koch et al. (2016)`, `McLoughlin et al. (2018)`, `Qiao et al. (2021)` |
| **Candidate targetability, conservation & selectivity rankings (`11_post_execution_audit_v3/`)** | **100% Computational Prediction Only** | Multi-isolate `tblastn`/semi-global alignments, `blastn`/21-mer scans, ViennaRNA folding, heuristic Pareto indices |

**Conclusion:** `BRD_CsSTE11_CAND013` (Region `CsSTE11_REG01`) and its companion non-redundant representatives (`BRD_CsSTE11_CAND014 [CsSTE11_REG03], BRD_CsSTE11_CAND010 [CsSTE11_REG02], BRD_CsPKS1_CAND012 [CsPKS1_REG02], BRD_CsPKS1_CAND010 [CsPKS1_REG04]`) are **prioritized computational candidates** for future laboratory validation — **not** proven fungicides, silencing reagents, or biosafe treatments.
