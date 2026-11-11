configfile: "config/project_configuration.yaml"

PY = "environments/brd_venv/Scripts/python.exe"

rule all:
    input:
        "10_final_reports/BRD_Final_Research_Report.html",
        "10_final_reports/BRD_Final_Research_Report.md",
        "10_final_reports/BRD_Final_Research_Report.pdf",
        "10_final_reports/BRD_Execution_Summary.txt",
        "10_final_reports/BRD_dsRNA_Candidate_Ranking.tsv",
        "10_final_reports/BRD_Off_Target_Analysis.tsv"

rule m01_acquisition:
    output:
        "01_data_acquisition/isolate_metadata.tsv",
        "01_data_acquisition/genome_metadata.tsv",
        "01_data_acquisition/download_manifest.tsv"
    shell:
        '"{PY}" BRD_pipeline.py --run m01'

rule m02_genome_qc:
    input:
        "01_data_acquisition/isolate_metadata.tsv"
    output:
        "02_genome_qc/assembly_quality_summary.tsv",
        "02_genome_qc/quast_summary.tsv",
        "02_genome_qc/busco_summary.tsv"
    shell:
        '"{PY}" BRD_pipeline.py --run m02'

rule m03_annotation:
    input:
        "01_data_acquisition/isolate_metadata.tsv",
        "02_genome_qc/assembly_quality_summary.tsv"
    output:
        "03_annotation/harmonised/gene_id_mapping.tsv",
        "03_annotation/cds/genome_cds.fasta",
        "03_annotation/proteins/genome_proteins.fasta",
        "03_annotation/transcripts/genome_transcripts.fasta"
    shell:
        '"{PY}" BRD_pipeline.py --run m03'

rule m04_comparative_genomics:
    input:
        "01_data_acquisition/isolate_metadata.tsv",
        "03_annotation/proteins/genome_proteins.fasta"
    output:
        "04_comparative_genomics/orthofinder/Orthogroups.tsv",
        "04_comparative_genomics/nucleotide_conservation/target_nucleotide_conservation.tsv"
    shell:
        '"{PY}" BRD_pipeline.py --run m04'

rule m05_functional_evidence:
    input:
        "03_annotation/proteins/genome_proteins.fasta",
        "01_data_acquisition/download_manifest.tsv"
    output:
        "05_functional_evidence/functional_evidence_summary/functional_evidence_summary.tsv"
    shell:
        '"{PY}" BRD_pipeline.py --run m05'

rule m06_target_prioritisation:
    input:
        "04_comparative_genomics/nucleotide_conservation/target_nucleotide_conservation.tsv",
        "05_functional_evidence/functional_evidence_summary/functional_evidence_summary.tsv"
    output:
        "06_target_prioritisation/target_scoring/target_prioritisation.tsv",
        "10_final_reports/BRD_Target_Gene_Accession_Verification.tsv"
    shell:
        '"{PY}" BRD_pipeline.py --run m06'

rule m07_dsrna_design:
    input:
        "06_target_prioritisation/target_scoring/target_prioritisation.tsv",
        "03_annotation/transcripts/genome_transcripts.fasta"
    output:
        "07_dsRNA_design/candidate_sequences/dsRNA_candidates.fasta",
        "07_dsRNA_design/candidate_sequences/dsRNA_candidates.tsv",
        "07_dsRNA_design/guide_sequences/guide_sequences.fasta"
    shell:
        '"{PY}" BRD_pipeline.py --run m07'

rule m08_off_target:
    input:
        "07_dsRNA_design/candidate_sequences/dsRNA_candidates.fasta",
        "07_dsRNA_design/guide_sequences/guide_sequences.fasta"
    output:
        "08_off_target_screening/summaries/off_target_summary.tsv",
        "10_final_reports/BRD_Off_Target_Analysis.tsv"
    shell:
        '"{PY}" BRD_pipeline.py --run m08'

rule m09_ranking:
    input:
        "07_dsRNA_design/candidate_sequences/dsRNA_candidates.tsv",
        "08_off_target_screening/summaries/off_target_summary.tsv",
        "05_functional_evidence/functional_evidence_summary/functional_evidence_summary.tsv"
    output:
        "09_candidate_ranking/integrated_scores/final_candidate_ranking.tsv",
        "10_final_reports/BRD_dsRNA_Candidate_Ranking.tsv"
    shell:
        '"{PY}" BRD_pipeline.py --run m09'

rule m10_reporting:
    input:
        "09_candidate_ranking/integrated_scores/final_candidate_ranking.tsv",
        "06_target_prioritisation/target_scoring/target_prioritisation.tsv"
    output:
        "10_final_reports/BRD_Final_Research_Report.html",
        "10_final_reports/BRD_Final_Research_Report.md",
        "10_final_reports/BRD_Final_Research_Report.pdf",
        "10_final_reports/BRD_Execution_Summary.txt"
    shell:
        '"{PY}" BRD_pipeline.py --run m10'
