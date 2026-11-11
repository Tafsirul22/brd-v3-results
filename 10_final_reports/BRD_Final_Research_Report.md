# Bipolaris-RNAi Designer (BRD) — Final Research Report

## 1. Title and Project Description

**Bipolaris-RNAi Designer (BRD) v1.0.0** is a modular, reproducible bioinformatics pipeline engineered to identify, verify, and computationally prioritise dsRNA candidates targeting pathogenicity-associated genes in *Bipolaris sorokiniana* (teleomorph *Cochliobolus sativus*). Generated autonomously on 2026-11-10 23:41 UTC.

## 2. Research Objectives

1. Verify local and public *B. sorokiniana* genome assemblies and annotations.
2. Resolve and sequence-verify the four target genes (`CsPKS1`, `CsFUS3`, `CsSTE11`, and `CsSp1`).
3. Measure within-species nucleotide and 21-mer conservation across all qualifying public isolates.
4. Integrate species-specific genetics, PHI-base v4.19, InterPro, UniProt, and GEO infection-stage RNA-seq expression (`GSE130517`).
5. Design candidate dsRNAs using LUCID (`sifi21_plus` + ViennaRNA), SIREN, dsRNAmax, and BRD isolate-aware sliding windows.
6. Screen candidates and 21-mer guides against hexaploid wheat (`IWGSC v2.1`), barley (`MorexV3`), and five beneficial fungal species.
7. Rank candidates with Pareto optimality, 1,000-iteration isolate bootstrap CIs, and weight sensitivity testing.

## 3. Biological Rationale

Spot blotch, common root rot, and black point caused by *B. sorokiniana* threaten wheat and barley production, particularly in warm humid regions such as South Asia (Al-Sadi, 2021; Roy et al., 2023). While transformation-based hairpin RNAi silenced endogenous `CsPKS1` in *B. sorokiniana* (Leng et al., 2011), demonstrating functional intracellular RNAi machinery, **no exogenous dsRNA uptake or spray-induced gene silencing (SIGS) experiment has been published for this species**. BRD frames exogenous dsRNA uptake as a testable experimental hypothesis rather than an established fact (Qiao et al., 2021).

## 4. Input Dataset Inventory

Recursive inspection of `C:\Users\ASUS BANGLADESH\Desktop\New folder (2)` identified 17 inventory records (local NCBI Datasets archive `Bipolaris sorokiniana.zip` for `GCF_000338995.1` ND90Pr verified against all 11 MD5 checksums, plus 3 project planning/review `.docx` documents). Original input files remained strictly unmodified.

## 5. Genome Accession and Isolate Information

Querying NCBI TaxID 45130 identified 20 assembly records representing **19 distinct biological isolates** after deduplicating paired GenBank/RefSeq `GCA_000338995.1` vs `GCF_000338995.1` (ND90Pr). Isolates span North America (`ND90Pr`, `ND93-1` T2T), China (`LK93`, `Yt-6`, Gansu), Australia (`WAI2406`, `WAI2411`, `BRIP10943a`, `BRIP27492a`), and India (`BS112` plus 9 recent 2026 assemblies).

## 6. Genome Quality Assessment

All 19 biological isolate assemblies passed FASTA validity, genome size (32.42–38.20 Mb), and GC% (49.24%–50.29%) checks. Contig/scaffold N50 ranged from 130.1 kb to 2.31 Mb. Note: Official BUSCO v6.1.0 and QUAST v5.3.0 require Linux/WSL2 (no WSL2 distro installed); BUSCO completeness is explicitly recorded as `NOT_RUN_LINUX_ONLY_NO_WSL_DISTRO` and never fabricated.

## 7. Annotation Harmonisation

RefSeq `GCF_000338995.1` (`ND90Pr`, JGI/NCBI annotation) passed structural GFF3 validation with 0 structural errors. Out of 12210 evaluated CDS sequences, 12202 translated with 100% identity to `protein.faa`. Re-annotation via Linux-only Funannotate/BRAKER3 was not run; unannotated isolate genomes were interrogated directly via splice-aware `tblastn` + `blastn` exon reconstruction in M04.

## 8. Comparative Genomics

DIAMOND all-vs-all RBH orthogroup clustering grouped annotated proteomes, and splice-aware genomic locus mapping reconstructed all 4 targets across the 19-isolate panel, followed by MAFFT v7.526 alignment.

| target   | orthogroup_id   |   isolates_evaluated |   isolates_with_gene |   isolates_complete_cds |   ref_cds_length_nt |   mean_pairwise_nt_identity_pct |   min_pairwise_nt_identity_pct |   mean_pairwise_aa_identity_pct |   snp_sites |   indel_columns |   conserved_blocks_ge21nt_count |   conserved_blocks_ge150nt_count |   longest_100pct_conserved_block_nt |   mean_21mer_isolate_conservation_pct |   pct_100pct_conserved_21mers | orthology_method                                                                                |
|:---------|:----------------|---------------------:|---------------------:|------------------------:|--------------------:|--------------------------------:|-------------------------------:|--------------------------------:|------------:|----------------:|--------------------------------:|---------------------------------:|------------------------------------:|--------------------------------------:|------------------------------:|:------------------------------------------------------------------------------------------------|
| CsPKS1   | BRD_OG000098    |                   19 |                   19 |                      17 |                6462 |                          99.958 |                         99.845 |                          96.465 |          22 |              23 |                              24 |                               14 |                                 852 |                                98.973 |                         89.77 | DIAMOND_all_vs_all_RBH_graph_clustering (OrthoFinder Linux/MCL unavailable natively on Windows) |
| CsFUS3   | BRD_OG000010    |                   19 |                   19 |                      17 |                1059 |                          99.984 |                         99.905 |                          94.43  |           3 |               5 |                               6 |                                3 |                                 364 |                                99.205 |                         88.35 | DIAMOND_all_vs_all_RBH_graph_clustering (OrthoFinder Linux/MCL unavailable natively on Windows) |
| CsSTE11  | BRD_OG000001    |                   19 |                   19 |                      17 |                2880 |                          99.657 |                         99.583 |                          95.785 |          19 |              16 |                              25 |                                6 |                                 273 |                                92.602 |                         78.5  | DIAMOND_all_vs_all_RBH_graph_clustering (OrthoFinder Linux/MCL unavailable natively on Windows) |
| CsSp1    | BRD_OG013358    |                   19 |                   19 |                      19 |                 270 |                          99.959 |                         99.259 |                          99.875 |           2 |               0 |                               3 |                                0 |                                 124 |                                99.116 |                         83.2  | DIAMOND_all_vs_all_RBH_graph_clustering (OrthoFinder Linux/MCL unavailable natively on Windows) |

## 9. Target Gene Identification

All four targets achieved **Confirmed identity** in `GCF_000338995.1`:
- `CsPKS1`: `COCSADRAFT_157080` (`XP_007696213.1`, 6,462 nt CDS, 100.0% match to `HQ830033.1`).
- `CsFUS3`: `COCSADRAFT_24574` (`XP_007697446.1`, RefSeq symbol `CHK1`, 1,059 nt CDS, confirmed by `CsFUS3-F1`/`CsFUS3-com-F/R` primers).
- `CsSTE11`: `COCSADRAFT_105370` (`XP_007694574.1`, RefSeq symbol `STE11`, 2,880 nt CDS, confirmed by `CsSTE11-F1/F4/F5/F6` primers).
- `CsSp1`: **Resolved** via 100% matches to 5 independent primers from Zhang et al. (2022, PMC8743017 Table 2) as `COCSADRAFT_164106` (`XP_007704210.1`, `XM_007706020.1`, 270 nt CDS, 89 aa).

## 10. Accession Verification

Complete 24-column accession verification including `HQ830033.1`, `KF769134.1`–`KF769141.1`, paired `GCA/GCF_000338995.1`, and `CsSp1` primer proof is saved in `10_final_reports/BRD_Target_Gene_Accession_Verification.tsv`.

## 11. Functional Evidence

Direct *B. sorokiniana* experimental genetics support `CsFUS3` and `CsSTE11` (Leng & Zhong, 2015), `CsSp1` (Zhang et al., 2022), and `CsPKS1` hairpin silencing (Leng et al., 2011). PHI-base v4.19 DIAMOND searches, InterPro/UniProt annotations, and GEO `GSE130517` barley-infection RNA-seq confirmed in-planta transcription of all 4 targets and identified 14 core RNAi machinery homologs (DCL, AGO, RdRP) in `ND90Pr`.

## 12. Target Prioritisation

Using provisional research weights (25% conservation, 25% off-target safety, 20% functional evidence, 15% expression, 15% designability), `CsFUS3` and `CsSTE11` lead on combined pathogenicity evidence, length flexibility, and expression, while `CsSp1` is constrained by its 270 nt CDS length and `CsPKS1` serves as the endogenous hairpin-RNAi benchmark.

## 13. dsRNA Design Methodology

Four complementary computational workflows were executed: (1) LUCID baseline (`sifi21_plus` free-energy/S1-S19 asymmetry + ViennaRNA 2.8.0 unpaired accessibility + Primer3, 200–400 nt), (2) SIREN v0.1.6 (`siren_prefilter.py` + ViennaRNA `duplexfold` + `siren_designVIII.py`), (3) `dsRNAmax.exe` v1.1.14 multi-isolate pan-CDS k-mer maximisation, and (4) BRD isolate-aware sliding windows across 180–500 nt.

## 14. Candidate Sequences and Characteristics

M07 generated **50 unique dsRNA candidates** (lengths 180–300 nt; GC 47.0%–58.9%) and 2876 unique 21-mer guide sequences.

## 15. Off-Target Screening

All 50 candidates and their 21-mer guides were screened via `blastn-short` (`word_size=7`), exact 21-mer matching, and ViennaRNA `RNA.duplexfold` against 7 RefSeq non-target transcriptomes: Hexaploid Wheat (`GCF_018294505.1`), Barley (`GCF_904849725.1`), *Trichoderma atroviride*, *T. virens*, *T. reesei*, *Beauveria bassiana*, and *Rhizophagus irregularis*. Absence of hits in these databases indicates no detected sequence matches under stated thresholds, **not** proven biological safety.

## 16. Candidate Ranking

Of 50 candidates evaluated in M09, **39 passed all computational filters** (0 exact 21-mer host/beneficial/paralog hits, QC pass, >=85% conservation), including **1 High computational priority**, **26 Moderate computational priority**, and **14 Pareto rank-1 non-dominated** candidates.

**Top 5 Ranked Candidates:**

|   final_rank | candidate_id        | target   | tool_name         |   candidate_length_nt |   conservation_pct |   off_target_safety_score |   composite_score | computational_priority_category   |
|-------------:|:--------------------|:---------|:------------------|----------------------:|-------------------:|--------------------------:|------------------:|:----------------------------------|
|            1 | BRD_CsSTE11_CAND012 | CsSTE11  | BRD_IsolateAware  |                   180 |             100    |                      0.55 |            0.7832 | High computational priority       |
|            2 | BRD_CsSp1_CAND008   | CsSp1    | BRD_IsolateAware  |                   180 |              98.91 |                      0.35 |            0.7458 | Moderate computational priority   |
|            3 | BRD_CsSp1_CAND001   | CsSp1    | LUCID_sifi21_plus |                   200 |              98.77 |                      0.35 |            0.7452 | Moderate computational priority   |
|            4 | BRD_CsSp1_CAND009   | CsSp1    | BRD_IsolateAware  |                   210 |              98.84 |                      0.35 |            0.7448 | Moderate computational priority   |
|            5 | BRD_CsSp1_CAND005   | CsSp1    | SIREN             |                   200 |              98.77 |                      0.35 |            0.7437 | Moderate computational priority   |

## 17. Uncertainty and Sensitivity Analysis

1,000-iteration bootstrap resampling across the 19 biological isolates quantified 95% confidence intervals for 21-mer conservation and rank stability. Five weight-perturbation scenarios (`baseline`, `safety_first`, `conservation_first`, `evidence_first`, `equal_weights`) evaluated sensitivity of the top shortlist.

## 18. Limitations

1. **Exogenous dsRNA uptake is unproven in *B. sorokiniana*** (Qiao et al., 2021).
2. **Linux-only tools (QUAST, BUSCO, OrthoFinder MCL, Funannotate/BRAKER3, RNAhybrid)** could not execute natively on Windows 11 without an installed WSL2 Linux distro; transparent Python/BLAST+/DIAMOND/ViennaRNA equivalents were run where scientifically valid, and BUSCO scores remain `NA`.
3. **RNA-seq (`GSE130517`)** represents *H. vulgare* root infection in `ND90Pr`, whereas `PRJNA743515` (wheat root rot) raw FASTQ quantification remains future work.
4. **Off-target screening** reflects the 7 RefSeq transcriptomes used and does not cover unannotated non-coding RNAs or untested soil microbiomes.

## 19. Recommended Candidates for Experimental Validation

Top computationally passing candidates per target gene selected for sequential bench evaluation:

- `CsSTE11`: `BRD_CsSTE11_CAND012` (BRD_IsolateAware, 180 nt, cons=100.0%, safety=0.55, score=0.7832, rank #1, High computational priority)
- `CsFUS3`: `BRD_CsFUS3_CAND013` (BRD_IsolateAware, 210 nt, cons=100.0%, safety=0.45, score=0.7281, rank #13, Moderate computational priority)
- `CsSp1`: `BRD_CsSp1_CAND008` (BRD_IsolateAware, 180 nt, cons=98.91%, safety=0.35, score=0.7458, rank #2, Moderate computational priority)
- `CsPKS1`: `BRD_CsPKS1_CAND013` (BRD_IsolateAware, 180 nt, cons=99.18%, safety=0.55, score=0.6545, rank #26, Moderate computational priority)

Sequential validation roadmap: (1) fluorescently labelled dsRNA uptake assays in *B. sorokiniana* conidia/hyphae (essential prerequisite, as exogenous dsRNA uptake is unproven in this species), (2) *in vitro* transcript knockdown qRT-PCR, and (3) detached wheat/barley leaf SIGS assays with scrambled dsRNA and water controls.

## 20. References

- Al-Sadi, A. M. (2021). *Frontiers in Cellular and Infection Microbiology*, 11, 584899. https://doi.org/10.3389/fcimb.2021.584899
- Appels, R., et al. (2018). *Science*, 361(6403), eaar7191. https://doi.org/10.1126/science.aar7191
- Fletcher, S. J., et al. (2025). dsRNAmax. *NAR Genomics and Bioinformatics*, 7(2), lqaf064. https://doi.org/10.1093/nargab/lqaf064
- Jiménez-Castro, L., et al. (2026). LUCID. *Plant Biotechnology Journal*, 24(6), 3597–3615. https://doi.org/10.1111/pbi.70526
- Koch, A., et al. (2016). *PLOS Pathogens*, 12(10), e1005901. https://doi.org/10.1371/journal.ppat.1005901
- Leng, Y., et al. (2011). *Molecular Plant Pathology*, 12(3), 289–298. https://doi.org/10.1111/j.1364-3703.2010.00666.x
- Leng, Y., & Zhong, S. (2015). *PLOS ONE*, 10(5), e0128291. https://doi.org/10.1371/journal.pone.0128291
- Leng, Y., et al. (2024). *PhytoFrontiers*, 4(4), 247–250. https://doi.org/10.1094/PHYTOFR-08-23-0108-A
- Piombo, E., et al. (2021). *Microbiology Spectrum*, 9(2), e0109921. https://doi.org/10.1128/Spectrum.01099-21
- Qiao, L., et al. (2021). *Plant Biotechnology Journal*, 19(9), 1756–1768. https://doi.org/10.1111/pbi.13589
- Roy, C., et al. (2023). *Frontiers in Plant Science*, 14, 1098648. https://doi.org/10.3389/fpls.2023.1098648
- Urban, M., et al. (2022). PHI-base. *Nucleic Acids Research*, 50(D1), D837–D847. https://doi.org/10.1093/nar/gkab1037
- Vargas Mejía, P., & Vega-Arreguín, J. C. (2025/2026). SIREN. https://doi.org/10.1101/2025.05.26.656188
- Zhang, W., et al. (2022). CsSp1. *Molecular Plant Pathology*, 23(2), 218–236. https://doi.org/10.1111/mpp.13155
- Zhang, W., et al. (2023). LK93 genome. *MPMI*, 36(7), 452–456. https://doi.org/10.1094/MPMI-09-22-0196-A

## 21. Software Versions and Parameters

| software               | version                                                                                              | repository_or_source                                | license                          | execution_status                                                                                   |
|:-----------------------|:-----------------------------------------------------------------------------------------------------|:----------------------------------------------------|:---------------------------------|:---------------------------------------------------------------------------------------------------|
| Python                 | 3.13.14                                                                                              | https://www.python.org/                             | PSF                              | Native Windows 11                                                                                  |
| Biopython              | 1.86                                                                                                 | https://biopython.org/                              | Biopython                        | Native Windows 11                                                                                  |
| ViennaRNA (Python RNA) | 2.8.0                                                                                                | https://www.tbi.univie.ac.at/RNA/                   | ViennaRNA                        | Native Windows 11                                                                                  |
| Snakemake              | 9.27.1                                                                                               | https://snakemake.readthedocs.io/                   | MIT                              | Native Windows 11                                                                                  |
| pandas                 | 2.3.3                                                                                                | https://pandas.pydata.org/                          | BSD-3-Clause                     | Native Windows 11                                                                                  |
| NumPy                  | 2.4.0                                                                                                | https://numpy.org/                                  | BSD-3-Clause                     | Native Windows 11                                                                                  |
| SciPy                  | 1.16.3                                                                                               | https://scipy.org/                                  | BSD-3-Clause                     | Native Windows 11                                                                                  |
| Matplotlib             | 3.10.8                                                                                               | https://matplotlib.org/                             | PSF-based                        | Native Windows 11                                                                                  |
| primer3-py             | 2.4.0                                                                                                | https://github.com/libnano/primer3-py               | GPL-2.0                          | Native Windows 11                                                                                  |
| reportlab              | 4.4.7                                                                                                | https://www.reportlab.com/                          | BSD                              | Native Windows 11                                                                                  |
| PyYAML                 | 6.0.3                                                                                                | https://pyyaml.org/                                 | MIT                              | Native Windows 11                                                                                  |
| NCBI BLAST+            | blastn: 2.17.0+ | Package: blast 2.17.0, build Jul  1 2025 08:57:20                                  | https://blast.ncbi.nlm.nih.gov/                     | Public Domain                    | Native Windows 11                                                                                  |
| NCBI Datasets CLI      | datasets version: 18.41.0                                                                            | https://www.ncbi.nlm.nih.gov/datasets/              | Public Domain                    | Native Windows 11                                                                                  |
| DIAMOND                | diamond version 2.2.8                                                                                | https://github.com/bbuchfink/diamond                | GPL-3.0                          | Native Windows 11                                                                                  |
| MAFFT                  | Preparing environment to run MAFFT on Windows. | This may take a while, if real-time scanning by ant | https://mafft.cbrc.jp/alignment/software/           | BSD-3-Clause                     | Native Windows 11                                                                                  |
| dsRNAmax               | 2026/11/11 05:41:33 dsRNAmax - dsRNA maximizer (Version: 1.1.14) | -biasHeader string                | https://github.com/sfletc/dsRNAmax                  | BSD-3-Clause                     | Native Windows 11 (commit 674b5b5)                                                                 |
| SIREN (siren-rnai)     | 0.1.6 (commit 0f5a5a9)                                                                               | https://github.com/pablovargasmejia/SIREN           | GPL-3.0                          | Native Windows Python modules (prefilter + designVIII); RNAhybrid replaced by ViennaRNA duplexfold |
| LUCID (sifi21_plus)    | 0.1.0 (commit 37910bd)                                                                               | https://github.com/Brunxi/LUCID                     | No explicit LICENSE file in repo | sifi21_plus free-energy/efficiency reproduced natively with ViennaRNA 2.8.0                        |
| tensorflow/lucid       | v0.3.10 (rejected)                                                                                   | https://github.com/tensorflow/lucid                 | Apache-2.0                       | EXCLUDED — neural network interpretability library, not RNAi design                                |
| evanpeikon/BLAST       | unreleased (rejected)                                                                                | https://github.com/evanpeikon/BLAST                 | MIT                              | EXCLUDED — educational script; official NCBI BLAST+ 2.17.0+ used                                   |
| metashot/busco         | 1.0.0 (rejected)                                                                                     | https://github.com/metashot/busco                   | GPL-3.0                          | EXCLUDED — deprecated Nextflow wrapper                                                             |
| QUAST                  | 5.3.0 (not installed)                                                                                | https://github.com/ablab/quast                      | GPL-2.0                          | Linux-only (WSL2 distro absent); Python QUAST-compatible metrics computed                          |
| BUSCO                  | 6.1.0 (not installed)                                                                                | https://busco.ezlab.org/                            | MIT                              | Linux-only (WSL2 distro absent); scores marked NA (not fabricated)                                 |
| OrthoFinder            | 3.2 (not installed)                                                                                  | https://github.com/OrthoFinder/OrthoFinder          | GPL-3.0                          | Linux MCL/FastME unavailable natively; DIAMOND RBH graph clustering executed                       |
| Funannotate / BRAKER3  | 1.9.1 / 3.0.8 (not installed)                                                                        | https://github.com/nextgenusfs/funannotate          | BSD-2 / Artistic                 | Linux-only; NCBI RefSeq/GenBank annotations reused + tblastn/blastn exon alignment                 |
| RNAhybrid              | 2.1.2 (not installed)                                                                                | https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid | GPL-2.0                          | Linux/macOS Bioconda only; ViennaRNA 2.8.0 RNA.duplexfold executed                                 |

## 22. Data Provenance

Complete SHA-256 provenance ledger saved at `C:\Users\ASUS BANGLADESH\Desktop\New folder (2)\Output\00_project_metadata\data_provenance.tsv` and `C:\Users\ASUS BANGLADESH\Desktop\New folder (2)\Output\00_project_metadata\provenance.json`.