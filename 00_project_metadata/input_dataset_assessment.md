# Input Dataset Assessment — Bipolaris-RNAi Designer (BRD)

**Input Directory:** `C:\Users\ASUS BANGLADESH\Desktop\New folder (2)`  
**Assessment Date:** 2026-11-11  
**Inventory File:** `00_project_metadata/input_file_inventory.tsv`

## 1. Summary of Discovered Input Files

A recursive inspection of the input directory identified four input items (excluding the designated `Output/` directory):

1. **`Bipolaris sorokiniana.zip`** (48,920,317 bytes; NCBI Datasets Data Package)
   - **Organism:** *Bipolaris sorokiniana* ND90Pr (NCBI TaxID: 665912; synonym *Cochliobolus sativus* ND90Pr)
   - **Assembly:** `GCF_000338995.1` (`Cocsa1`, paired GenBank assembly `GCA_000338995.1`), RefSeq Reference Genome, Scaffold-level (154 scaffolds, 504 contigs, 34,409,167 bp, contig N50 = 243,393 bp, scaffold N50 = 1,789,485 bp, GC = 50%).
   - **Annotation:** Annotation submitted by US DOE Joint Genome Institute (JGI), propagated via NCBI Eukaryotic Annotation Propagation Pipeline (2020-02-06); 12,304 total genes (12,210 protein-coding, 92 non-coding, 2 pseudogenes).
   - **BioProject / BioSample:** `PRJNA53923` / `SAMN02981285` (SRA `SRS1566654`, WGS `AEIN00000000.1`).
   - **Contained sequence & annotation files:**
     - `GCF_000338995.1_Cocsa1_genomic.fna` (34,857,722 bytes) — Genomic nucleotide FASTA (M01, M02, M03, M06)
     - `genomic.gff` (29,296,335 bytes) — GFF3 structural & functional annotation (M03, M06)
     - `genomic.gtf` (38,616,812 bytes) — GTF annotation (M03)
     - `genomic.gbff` (62,660,172 bytes) — GenBank flat file (M01, M03, M05, M06)
     - `cds_from_genomic.fna` (19,813,523 bytes) — Annotated CDS nucleotide FASTA (M03, M04, M06, M07)
     - `rna.fna` (20,359,143 bytes) — Transcript/RNA nucleotide FASTA (M03, M07)
     - `protein.faa` (6,582,504 bytes) — Annotated protein FASTA (M03, M04, M05, M06)
     - `assembly_data_report.jsonl`, `sequence_report.jsonl`, `data_summary.tsv`, `dataset_catalog.json`, `md5sum.txt` — Metadata & integrity checksums (M01, M02).
   - **Suitability:** Complete, integrity-verified primary reference genome package suitable for M01–M07 without requiring *de novo* re-annotation of ND90Pr.

2. **`BRD_Complete_Pipeline_Technical_Report_v1.0 (1).docx`** (361,743 bytes)
   - Technical engineering specification for BRD v1.0, defining Snakemake architecture, manifest schemas, LUCID (`https://github.com/Brunxi/LUCID`), SIREN (`https://github.com/pablovargasmejia/SIREN`), and dsRNAmax (`https://github.com/sfletc/dsRNAmax`) integration roles, Pareto/bootstrap ranking, and prediction-only scientific boundaries.

3. **`BRD_Literature_Review.docx`** (1,642,956 bytes)
   - Critical literature synthesis covering RNAi machinery in fungi, SIGS vs transformation-based hairpin RNAi (Leng et al., 2011), *B. sorokiniana* molecular/genomic resources (Aggarwal et al., 2019, 2022; Zhang et al., 2022, 2023; Leng et al., 2024; Basak et al., 2024; Shukla et al., 2025), off-target considerations in hexaploid wheat and beneficial fungi, and 82 verified references.

4. **`Bsorokiniana_four_target_comparison.docx`** (45,010 bytes)
   - Provisional literature-based comparison of the four target genes (`CsPKS1`, `CsFUS3`, `CsSTE11`, `CsSp1`), establishing `CsFUS3` and `CsSTE11` as distinct joint-lead MAPK candidates, `CsSp1` as a wheat-colonisation effector comparator, and `CsPKS1` as the endogenous hairpin-RNAi benchmark.

## 2. Dataset Gaps & Expansion Strategy

- **Single local isolate:** The local zip contains only one isolate (`ND90Pr`, `GCF_000338995.1`). Multi-isolate comparative genomics (M04) and bootstrap uncertainty analysis (M09) require expanding the panel using public *B. sorokiniana* assemblies from NCBI (e.g., chromosome-level `LK93` `GCA_024722535.1`, `ND93-1` T2T `GCA_033006315.1`, Australian isolates `WAI2406`, `BRIP10943a`, `BRIP27492a`, `WAI2411`, Indian isolate `BS112` `GCA_004329375.1`, and `Yt-6` `GCA_037974985.1`, plus additional NCBI isolates).
- **Host & non-target references:** No wheat, barley, or beneficial fungi FASTA files are present in the input directory; M08 must acquire versioned reference transcript/CDS/genome datasets from NCBI/EnsemblPlants with recorded accessions and checksums.
- **Immutability:** All original input files remain unmodified; `Bipolaris sorokiniana.zip` will be extracted into `Output/01_data_acquisition/` with MD5 verification against `md5sum.txt`.
