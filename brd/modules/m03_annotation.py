"""Module M03 — Annotation Harmonisation & Validation."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, List

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from brd.config import BRDConfig
from brd.integrations.funannotate import check_braker, check_funannotate
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker
from brd.sequence.cds_validation import validate_cds_translation
from brd.sequence.gff_parser import validate_gff3
from brd.utils.report_utils import df_to_html_table, wrap_html_report, write_text


def run_m03(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    logger = get_logger("BRD.M03", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    m01_iso = cfg.output_dir / "01_data_acquisition" / "isolate_metadata.tsv"
    df_iso = pd.read_csv(m01_iso, sep="\t")
    m03_dir = cfg.module_dir("03_annotation")
    orig_dir = m03_dir / "original_annotations"
    harm_dir = m03_dir / "harmonised"
    gff_out_dir = m03_dir / "gff3"
    cds_out_dir = m03_dir / "cds"
    tx_out_dir = m03_dir / "transcripts"
    prot_out_dir = m03_dir / "proteins"
    qc_dir = m03_dir / "annotation_qc"
    for d in (orig_dir, harm_dir, gff_out_dir, cds_out_dir, tx_out_dir, prot_out_dir, qc_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Inspect primary reference GCF_000338995.1 (ND90Pr)
    ref_row = df_iso[df_iso["assembly_accession"] == "GCF_000338995.1"].iloc[0]
    ref_fna = Path(str(ref_row["genome_fasta"]))
    ref_gff = Path(str(ref_row["gff_path"]))
    ref_pkg = ref_gff.parent
    ref_cds = ref_pkg / "cds_from_genomic.fna"
    ref_rna = ref_pkg / "rna.fna"
    ref_prot = ref_pkg / "protein.faa"

    # Preserve original annotations
    shutil.copy2(ref_gff, orig_dir / "GCF_000338995.1_genomic.gff")
    seq_lengths = {r.id: len(r.seq) for r in SeqIO.parse(str(ref_fna), "fasta")}

    logger.info("Validating primary reference GFF3: %s", ref_gff)
    gff_info = validate_gff3(ref_gff, seq_lengths)
    prot_dict = {r.id: str(r.seq) for r in SeqIO.parse(str(ref_prot), "fasta")}
    prot_recs = list(SeqIO.parse(str(ref_prot), "fasta"))
    cds_recs = list(SeqIO.parse(str(ref_cds), "fasta"))
    rna_recs = list(SeqIO.parse(str(ref_rna), "fasta"))

    # Parse CDS headers and validate translation against protein.faa
    mapping_rows: List[dict] = []
    harmonised_cds: List[SeqRecord] = []
    harmonised_prot: List[SeqRecord] = []
    seen_prot_ids = set()

    cds_val_stats = {"total": 0, "complete": 0, "exact_prot_match": 0, "internal_stops": 0, "mod3_fail": 0}
    import re

    for idx, rec in enumerate(cds_recs, 1):
        desc = rec.description
        m_prot = re.search(r"\[protein_id=([^\]]+)\]", desc)
        m_locus = re.search(r"\[locus_tag=([^\]]+)\]", desc)
        m_gene = re.search(r"\[gene=([^\]]+)\]", desc)
        m_pname = re.search(r"\[protein=([^\]]+)\]", desc)
        m_loc = re.search(r"\[location=([^\]]+)\]", desc)

        pid = m_prot.group(1) if m_prot else rec.id
        locus = m_locus.group(1) if m_locus else ""
        gene_sym = m_gene.group(1) if m_gene else ""
        pname = m_pname.group(1) if m_pname else ""
        loc_str = m_loc.group(1) if m_loc else ""
        exp_prot = prot_dict.get(pid)
        val = validate_cds_translation(str(rec.seq), expected_protein=exp_prot)

        cds_val_stats["total"] += 1
        if val["is_complete_cds"]:
            cds_val_stats["complete"] += 1
        if val["matches_reference_protein"]:
            cds_val_stats["exact_prot_match"] += 1
        if val["internal_stop_count"] > 0:
            cds_val_stats["internal_stops"] += 1
        if not val["mod3_ok"]:
            cds_val_stats["mod3_fail"] += 1

        std_id = f"BRD_ND90PR_{locus or f'G{idx:05d}'}"
        mapping_rows.append(
            {
                "brd_gene_id": std_id,
                "assembly_accession": "GCF_000338995.1",
                "isolate_id": "ND90Pr",
                "original_locus_tag": locus,
                "gene_symbol": gene_sym,
                "protein_id": pid,
                "cds_header_id": rec.id,
                "genomic_location": loc_str,
                "product": pname,
                "cds_length_nt": val["cds_length_nt"],
                "protein_length_aa": val["translated_aa_len"],
                "is_complete_cds": val["is_complete_cds"],
                "matches_reference_protein": val["matches_reference_protein"],
                "internal_stop_count": val["internal_stop_count"],
            }
        )
        harmonised_cds.append(SeqRecord(rec.seq, id=std_id, description=f"protein_id={pid} locus_tag={locus} gene={gene_sym}"))
        if pid in prot_dict and pid not in seen_prot_ids:
            seen_prot_ids.add(pid)
            harmonised_prot.append(SeqRecord(Seq(prot_dict[pid]), id=std_id, description=f"protein_id={pid} locus_tag={locus} gene={gene_sym}"))

    # Export harmonised files
    harm_gff_p = harm_dir / "harmonised_annotation.gff3"
    shutil.copy2(ref_gff, harm_gff_p)
    shutil.copy2(ref_gff, gff_out_dir / "GCF_000338995.1_ND90Pr.gff3")

    genes_fa_p = harm_dir / "genome_genes.fasta"
    cds_fa_p = cds_out_dir / "genome_cds.fasta"
    tx_fa_p = tx_out_dir / "genome_transcripts.fasta"
    prot_fa_p = prot_out_dir / "genome_proteins.fasta"
    shutil.copy2(ref_cds, cds_fa_p)
    shutil.copy2(ref_cds, genes_fa_p)
    shutil.copy2(ref_rna, tx_fa_p)
    shutil.copy2(ref_prot, prot_fa_p)
    SeqIO.write(harmonised_cds, str(harm_dir / "harmonised_cds.fasta"), "fasta")
    SeqIO.write(harmonised_prot, str(harm_dir / "harmonised_proteins.fasta"), "fasta")

    # Copy any secondary annotated genome (e.g., GCA_013416765.1)
    annot_status_rows = []
    funa_st = check_funannotate()
    brak_st = check_braker()
    for _, r in df_iso.iterrows():
        acc = str(r["assembly_accession"])
        has_gff = bool(r["gff_path"]) and str(r["gff_path"]) != "nan" and Path(str(r["gff_path"])).exists()
        has_prot = bool(r["protein_path"]) and str(r["protein_path"]) != "nan" and Path(str(r["protein_path"])).exists()
        if has_gff and Path(str(r["gff_path"])) != ref_gff:
            shutil.copy2(str(r["gff_path"]), gff_out_dir / f"{acc}.gff3")
        if has_prot and Path(str(r["protein_path"])) != ref_prot:
            shutil.copy2(str(r["protein_path"]), prot_out_dir / f"{acc}_protein.faa")

        if acc == "GCF_000338995.1":
            st = "REUSED_PRIMARY_REFSEQ_JGI"
            note = f"GFF errors={gff_info['error_count']}, orphan_mRNAs={gff_info['orphan_mrnas']}, CDS exact protein match={cds_val_stats['exact_prot_match']}/{cds_val_stats['total']}"
        elif has_gff and has_prot:
            st = "REUSED_GENBANK_ANNOTATION"
            note = "Public GenBank structural annotation retained; proteomes included in M04 comparative genomics."
        else:
            st = "UNANNOTATED_GENOME_LINUX_TOOL_REQUIRED_FOR_FULL_GFF"
            note = f"Funannotate/BRAKER3 require Linux/WSL2 ({funa_st['status']}); M04 maps target loci directly via tblastn+blastn exon alignment."
        annot_status_rows.append(
            {
                "isolate_id": r["isolate_id"],
                "assembly_accession": acc,
                "has_gff3": has_gff,
                "has_protein_faa": has_prot,
                "annotation_strategy": st,
                "notes": note,
            }
        )

    map_df = pd.DataFrame(mapping_rows)
    map_p = harm_dir / "gene_id_mapping.tsv"
    map_df.to_csv(map_p, sep="\t", index=False)

    val_df = pd.DataFrame(annot_status_rows)
    val_p = qc_dir / "annotation_validation.tsv"
    val_df.to_csv(val_p, sep="\t", index=False)

    summary_metrics = pd.DataFrame(
        [
            {"metric": "GFF3 genes parsed", "value": len(gff_info["genes"])},
            {"metric": "GFF3 mRNAs parsed", "value": len(gff_info["mrnas"])},
            {"metric": "GFF3 orphan mRNAs", "value": gff_info["orphan_mrnas"]},
            {"metric": "GFF3 structural errors", "value": gff_info["error_count"]},
            {"metric": "Total CDS evaluated", "value": cds_val_stats["total"]},
            {"metric": "Complete CDS (ATG..stop, mod3, 0 internal stops)", "value": cds_val_stats["complete"]},
            {"metric": "CDS translations 100% matching protein.faa", "value": cds_val_stats["exact_prot_match"]},
            {"metric": "CDS with internal stops", "value": cds_val_stats["internal_stops"]},
            {"metric": "Funannotate status", "value": funa_st["status"]},
            {"metric": "BRAKER3 status", "value": brak_st["status"]},
        ]
    )
    summary_metrics.to_csv(qc_dir / "nd90pr_cds_translation_qc.tsv", sep="\t", index=False)

    rep_html = qc_dir / "annotation_qc_report.html"
    body = (
        f"<h2>Primary Reference Annotation QC (GCF_000338995.1 / ND90Pr)</h2>"
        f"{df_to_html_table(summary_metrics)}"
        f"<h2>Isolate Panel Annotation Inventory</h2>"
        f"{df_to_html_table(val_df)}"
    )
    write_text(rep_html, wrap_html_report("BRD M03 — Annotation Harmonisation & Validation", "GCF_000338995.1 & Isolate Panel", body))

    prov.record("M03", "Annotation Harmonisation", str(ref_gff), map_p, "Biopython+GFF3-QC", "1.86", "GCF_000338995.1", f"{len(map_df)} CDS validated")
    logger.info("M03 complete: %d CDS validated (%d exact protein matches).", cds_val_stats["total"], cds_val_stats["exact_prot_match"])
    return {
        "harmonised_gff3": harm_gff_p,
        "gene_id_mapping": map_p,
        "annotation_validation": val_p,
        "genome_cds": cds_fa_p,
        "genome_proteins": prot_fa_p,
        "genome_transcripts": tx_fa_p,
        "report": rep_html,
    }
