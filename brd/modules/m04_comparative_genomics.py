"""Module M04 — Comparative Genomics, Orthology, and Nucleotide Conservation Analysis."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from brd.config import BRDConfig
from brd.integrations.blast import make_blast_db, run_blast
from brd.integrations.diamond import build_diamond_db, run_diamond_blastp
from brd.integrations.orthofinder import run_comparative_orthogroups
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker
from brd.sequence.cds_validation import validate_cds_translation
from brd.sequence.conservation import analyze_msa_conservation
from brd.sequence.sequence_alignment import run_mafft
from brd.utils.report_utils import df_to_html_table, wrap_html_report, write_text


def _extract_reference_target_models(cfg: BRDConfig, gff_path: Path, genome_fna: Path, cds_fna: Path, prot_faa: Path) -> Dict[str, dict]:
    """Extract verified ND90Pr reference gene models, exon coordinates, genomic span, CDS, and protein."""
    genome_dict = {r.id: r.seq for r in SeqIO.parse(str(genome_fna), "fasta")}
    cds_by_prot = {}
    for rec in SeqIO.parse(str(cds_fna), "fasta"):
        m = re.search(r"\[protein_id=([^\]]+)\]", rec.description)
        if m:
            cds_by_prot[m.group(1)] = rec
    prot_dict = {r.id: r for r in SeqIO.parse(str(prot_faa), "fasta")}

    # Parse CDS exon intervals from GFF3 per protein_id
    exons_by_prot: Dict[str, List[Tuple[str, int, int, str, int]]] = {}
    with open(gff_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            m = re.search(r"Name=([^;]+)", parts[8])
            if not m:
                continue
            pid = m.group(1)
            exons_by_prot.setdefault(pid, []).append((parts[0], int(parts[3]), int(parts[4]), parts[6], int(parts[7]) if parts[7].isdigit() else 0))

    models: Dict[str, dict] = {}
    for tname, tcfg in cfg.raw["targets"].items():
        pid = tcfg["nd90pr_protein_id"]
        locus = tcfg["nd90pr_locus_tag"]
        cds_rec = cds_by_prot[pid]
        prot_rec = prot_dict[pid]
        exons = sorted(exons_by_prot[pid], key=lambda x: x[1])
        seqid = exons[0][0]
        strand = exons[0][3]
        g_start = min(e[1] for e in exons)
        g_end = max(e[2] for e in exons)
        genomic_seq = genome_dict[seqid][g_start - 1 : g_end]
        if strand == "-":
            genomic_seq = genomic_seq.reverse_complement()
            exons_ordered = sorted(exons, key=lambda x: x[1], reverse=True)
        else:
            exons_ordered = exons

        exon_seqs = []
        exon_rel_intervals = []
        for scf, es, ee, st, ph in exons_ordered:
            es_seq = genome_dict[scf][es - 1 : ee]
            if st == "-":
                es_seq = es_seq.reverse_complement()
                rel_s = g_end - ee
                rel_e = g_end - es + 1
            else:
                rel_s = es - g_start
                rel_e = ee - g_start + 1
            exon_seqs.append(str(es_seq).upper())
            exon_rel_intervals.append((rel_s, rel_e))

        models[tname] = {
            "target": tname,
            "locus_tag": locus,
            "protein_id": pid,
            "gene_id": tcfg["nd90pr_gene_id"],
            "scaffold": seqid,
            "start": g_start,
            "end": g_end,
            "strand": strand,
            "cds_seq": str(cds_rec.seq).upper(),
            "cds_id": cds_rec.id,
            "protein_seq": str(prot_rec.seq),
            "genomic_seq": str(genomic_seq).upper(),
            "exon_seqs": exon_seqs,
            "exon_rel_intervals": exon_rel_intervals,
        }
    return models


def _splice_align_target_in_genome(
    model: dict,
    iso_id: str,
    iso_acc: str,
    iso_fna: Path,
    blast_db_prefix: Path,
    tblastn_bin: str,
    blastn_bin: str,
    work_dir: Path,
    log_path: Path,
) -> dict:
    """Locate target gene in an isolate genome via tblastn + exon-level nucleotide alignment."""
    tname = model["target"]
    if iso_acc == "GCF_000338995.1":
        return {
            "target": tname,
            "isolate_id": iso_id,
            "assembly_accession": iso_acc,
            "status": "PRESENT_COMPLETE_CDS" if tname != "CsPKS1" else "PRESENT_COMPLETE_CDS_REFSEQ_PARTIAL_FLAG",
            "classification": "Present (Primary Reference)",
            "scaffold": model["scaffold"],
            "start": model["start"],
            "end": model["end"],
            "strand": model["strand"],
            "cds_length_nt": len(model["cds_seq"]),
            "protein_length_aa": len(model["protein_seq"]),
            "tblastn_pident": 100.0,
            "tblastn_qcov": 100.0,
            "cds_seq": model["cds_seq"],
            "protein_seq": model["protein_seq"],
            "touches_gap_or_edge": False,
            "notes": "Primary reference ND90Pr",
        }

    q_prot_fa = work_dir / f"{tname}_ref_prot.faa"
    if not q_prot_fa.exists():
        q_prot_fa.write_text(f">{tname}\n{model['protein_seq']}\n", encoding="utf-8")

    q_gen_fa = work_dir / f"{tname}_ref_genomic.fna"
    if not q_gen_fa.exists():
        q_gen_fa.write_text(f">{tname}\n{model['genomic_seq']}\n", encoding="utf-8")

    tb_out = work_dir / f"{iso_acc}_{tname}_tblastn.tsv"
    if tb_out.exists() and tb_out.stat().st_size > 0:
        from brd.integrations.blast import BLAST_COLS
        tb_df = pd.read_csv(tb_out, sep="\t", names=BLAST_COLS)
    else:
        tb_df = run_blast(tblastn_bin, q_prot_fa, blast_db_prefix, tb_out, evalue=1e-5, max_target_seqs=10, num_threads=4, log_path=log_path)

    bn_out = work_dir / f"{iso_acc}_{tname}_genomic_blastn.tsv"
    if bn_out.exists() and bn_out.stat().st_size > 0:
        from brd.integrations.blast import BLAST_COLS
        bn_df = pd.read_csv(bn_out, sep="\t", names=BLAST_COLS)
    else:
        bn_df = run_blast(blastn_bin, q_gen_fa, blast_db_prefix, bn_out, evalue=1e-10, max_target_seqs=10, num_threads=4, log_path=log_path)

    if tb_df.empty and bn_df.empty:
        return {
            "target": tname,
            "isolate_id": iso_id,
            "assembly_accession": iso_acc,
            "status": "TRUE_ABSENCE_NO_BLAST_HIT",
            "classification": "True gene absence (no tblastn/blastn hit)",
            "scaffold": "",
            "start": 0,
            "end": 0,
            "strand": "",
            "cds_length_nt": 0,
            "protein_length_aa": 0,
            "tblastn_pident": 0.0,
            "tblastn_qcov": 0.0,
            "cds_seq": "",
            "protein_seq": "",
            "touches_gap_or_edge": False,
            "notes": "Zero hits against full assembly",
        }

    # Select best locus from genomic blastn (preferred) or tblastn
    if not bn_df.empty:
        top_bn = bn_df.sort_values("bitscore", ascending=False).iloc[0]
        best_scf = str(top_bn["sseqid"])
        same_scf_bn = bn_df[bn_df["sseqid"] == best_scf]
        loc_min = int(min(same_scf_bn["sstart"].min(), same_scf_bn["send"].min()))
        loc_max = int(max(same_scf_bn["sstart"].max(), same_scf_bn["send"].max()))
        strand = "+" if top_bn["sstrand"] == "plus" else "-"
    else:
        top_tb = tb_df.sort_values("bitscore", ascending=False).iloc[0]
        best_scf = str(top_tb["sseqid"])
        same_scf_tb = tb_df[tb_df["sseqid"] == best_scf]
        loc_min = int(min(same_scf_tb["sstart"].min(), same_scf_tb["send"].min()))
        loc_max = int(max(same_scf_tb["sstart"].max(), same_scf_tb["send"].max()))
        strand = "+" if int(top_tb["sstart"]) <= int(top_tb["send"]) else "-"

    # Extract bounding locus (+/- 300 bp flank) from isolate genome
    scf_rec = None
    for rec in SeqIO.parse(str(iso_fna), "fasta"):
        if rec.id == best_scf:
            scf_rec = rec
            break
    if scf_rec is None:
        raise ValueError(f"Scaffold {best_scf} not found in {iso_fna}")

    pad = 150
    ext_start = max(1, loc_min - pad)
    ext_end = min(len(scf_rec.seq), loc_max + pad)
    locus_seq = scf_rec.seq[ext_start - 1 : ext_end]
    if strand == "-":
        locus_seq = locus_seq.reverse_complement()
    locus_str = str(locus_seq).upper()

    touches_edge = (loc_min <= 50) or (loc_max >= len(scf_rec.seq) - 50) or ("NNN" in locus_str)

    # Full-gene semi-global alignment: anchor ref genomic_seq (exons + introns) inside locus_str
    # then project exact ND90Pr exon intervals onto isolate coordinates
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    # Semi-global (target=locus_str > query=ref_gen): in Biopython PairwiseAligner,
    # gaps in query at ends of longer target are "deletion" (set 0.0 for free flank overhang),
    # while gaps in target at ends of query are "insertion" (keep -5/-1 so full ref_gen is anchored).
    aligner.open_left_deletion_score = 0.0
    aligner.extend_left_deletion_score = 0.0
    aligner.open_right_deletion_score = 0.0
    aligner.extend_right_deletion_score = 0.0
    aligner.open_left_insertion_score = -5.0
    aligner.extend_left_insertion_score = -1.0
    aligner.open_right_insertion_score = -5.0
    aligner.extend_right_insertion_score = -1.0

    ref_gen = model["genomic_seq"]
    aln = aligner.align(locus_str, ref_gen)[0]
    # Walk aligned segments (target_blocks, query_blocks) to build coordinate map from ref_gen -> locus_str
    t_blocks = aln.aligned[0]
    q_blocks = aln.aligned[1]
    ex_intervals = model["exon_rel_intervals"]

    # Map each ref_gen position inside an exon to the corresponding isolate bases (preserving isolate insertions inside exons)
    reconstructed_exons = []
    for ex_s, ex_e in ex_intervals:
        mapped_positions = []
        for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
            ov_s = max(qs, ex_s)
            ov_e = min(qe, ex_e)
            if ov_s < ov_e:
                t_ov_s = ts + (ov_s - qs)
                t_ov_e = ts + (ov_e - qs)
                mapped_positions.append((t_ov_s, t_ov_e))
        if mapped_positions:
            # Drop only stray micro-blocks (<=6 bp) separated by large (>30 bp) jumps (never drop the main exon block)
            if (
                len(mapped_positions) > 1
                and (mapped_positions[1][0] - mapped_positions[0][1]) > 30
                and (mapped_positions[0][1] - mapped_positions[0][0]) <= 6
            ):
                mapped_positions = mapped_positions[1:]
            if (
                len(mapped_positions) > 1
                and (mapped_positions[-1][0] - mapped_positions[-2][1]) > 30
                and (mapped_positions[-1][1] - mapped_positions[-1][0]) <= 6
            ):
                mapped_positions = mapped_positions[:-1]
            # Reconstruct exon preserving small internal indels (<=30 bp)
            ex_seq_parts = []
            prev_end = mapped_positions[0][0]
            for ts_i, te_i in mapped_positions:
                if ts_i > prev_end and (ts_i - prev_end) <= 30:
                    ex_seq_parts.append(locus_str[prev_end:ts_i])
                ex_seq_parts.append(locus_str[ts_i:te_i])
                prev_end = te_i
            reconstructed_exons.append("".join(ex_seq_parts))

    spliced_cds = "".join(reconstructed_exons)
    val = validate_cds_translation(spliced_cds)

    tb_pident = float(tb_df["pident"].max()) if not tb_df.empty else 0.0
    tb_qcov = float(tb_df["qcovs"].max()) if not tb_df.empty else 0.0
    prot_seq = ""
    if len(spliced_cds) >= 3 and not val["has_ambiguous_bases"]:
        trim = len(spliced_cds) - (len(spliced_cds) % 3)
        prot_seq = str(Seq(spliced_cds[:trim]).translate(table=1)).rstrip("*")

    if len(reconstructed_exons) == len(model["exon_seqs"]) and val["mod3_ok"] and val["internal_stop_count"] == 0 and not val["has_ambiguous_bases"]:
        status = "PRESENT_COMPLETE_CDS"
        cls_label = "Present (Complete CDS verified)"
    elif touches_edge:
        status = "ASSEMBLY_GAP_OR_CONTIG_EDGE"
        cls_label = "Incomplete assembly (contig edge or N-gap)"
    else:
        status = "PRESENT_PARTIAL_OR_VARIANT_CDS"
        cls_label = "Present (Partial/frameshift/variant CDS)"

    return {
        "target": tname,
        "isolate_id": iso_id,
        "assembly_accession": iso_acc,
        "status": status,
        "classification": cls_label,
        "scaffold": best_scf,
        "start": loc_min,
        "end": loc_max,
        "strand": strand,
        "cds_length_nt": len(spliced_cds),
        "protein_length_aa": len(prot_seq),
        "tblastn_pident": round(tb_pident, 2),
        "tblastn_qcov": round(tb_qcov, 2),
        "cds_seq": spliced_cds,
        "protein_seq": prot_seq,
        "touches_gap_or_edge": touches_edge,
        "notes": f"exons_matched={len(reconstructed_exons)}/{len(model['exon_seqs'])}; internal_stops={val['internal_stop_count']}; mod3={val['mod3_ok']}",
    }


def run_m04(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    cmd_log = cfg.output_dir / "logs" / "commands.jsonl"
    logger = get_logger("BRD.M04", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    m04_dir = cfg.module_dir("04_comparative_genomics")
    of_dir = m04_dir / "orthofinder"
    dmnd_dir = m04_dir / "diamond"
    mafft_dir = m04_dir / "mafft"
    core_dir = m04_dir / "core_accessory"
    ortho_dir = m04_dir / "target_orthologs"
    cons_dir = m04_dir / "nucleotide_conservation"
    fig_dir = m04_dir / "figures"
    for d in (of_dir, dmnd_dir, mafft_dir, core_dir, ortho_dir, cons_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    df_iso = pd.read_csv(cfg.output_dir / "01_data_acquisition" / "isolate_metadata.tsv", sep="\t")

    # 1. Orthogroup clustering across annotated proteomes + Botrytis cinerea (LUCID comparator)
    proteomes: Dict[str, Path] = {}
    for _, r in df_iso.iterrows():
        pp = str(r.get("protein_path", ""))
        if pp and pp != "nan" and Path(pp).exists():
            proteomes[f"Bsor_{r['isolate_id']}"] = Path(pp)
    bc_pep = Path(cfg.tool_path("lucid_repo")) / "data" / "references" / "Botrytis_cinerea.ASM83294v1.pep.all.fa"
    if bc_pep.exists():
        proteomes["Bcin_ASM83294v1"] = bc_pep

    logger.info("Running comparative orthogroup clustering across %d annotated proteomes...", len(proteomes))
    og_df, pa_df, og_method = run_comparative_orthogroups(
        proteomes,
        diamond_bin=cfg.tool_path("diamond"),
        out_dir=of_dir,
        log_path=cmd_log,
    )
    # Classify core vs accessory across B. sorokiniana annotated proteomes
    bsor_cols = [c for c in pa_df.columns if c.startswith("Bsor_")]
    core_rows = []
    for _, r in pa_df.iterrows():
        bs_present = sum(1 for c in bsor_cols if int(r[c]) > 0)
        tot_bs = sum(int(r[c]) for c in bsor_cols)
        if tot_bs == 0:
            continue
        cls_og = "Core (all annotated B. sorokiniana)" if bs_present == len(bsor_cols) else "Accessory (annotated B. sorokiniana)"
        core_rows.append({"Orthogroup": r["Orthogroup"], "bsor_proteomes_present": bs_present, "total_bsor_genes": tot_bs, "category": cls_og})
    core_df = pd.DataFrame(core_rows)
    core_p = core_dir / "core_accessory_genes.tsv"
    core_df.to_csv(core_p, sep="\t", index=False)

    # 2. Extract ND90Pr reference target models
    ref_row = df_iso[df_iso["assembly_accession"] == "GCF_000338995.1"].iloc[0]
    ref_pkg = Path(str(ref_row["gff_path"])).parent
    models = _extract_reference_target_models(
        cfg,
        gff_path=Path(str(ref_row["gff_path"])),
        genome_fna=Path(str(ref_row["genome_fasta"])),
        cds_fna=ref_pkg / "cds_from_genomic.fna",
        prot_faa=ref_pkg / "protein.faa",
    )

    # 3. Map all 4 targets across ALL 19 isolate genomes via BLAST+ & splice-aware exon extraction
    db_dir = cfg.output_dir / "intermediate" / "isolate_blast_dbs"
    work_dir = cfg.output_dir / "intermediate" / "target_extraction"
    db_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    makeblastdb_bin = cfg.tool_path("makeblastdb")
    tblastn_bin = cfg.tool_path("tblastn")
    blastn_bin = cfg.tool_path("blastn")
    mafft_bin = cfg.tool_path("mafft")

    homolog_rows = []
    for _, r in df_iso.iterrows():
        acc = str(r["assembly_accession"])
        iso_id = str(r["isolate_id"])
        iso_fna = Path(str(r["genome_fasta"]))
        db_pref = db_dir / acc
        make_blast_db(makeblastdb_bin, iso_fna, db_pref, dbtype="nucl", log_path=cmd_log)
        for tname, model in models.items():
            res = _splice_align_target_in_genome(
                model, iso_id, acc, iso_fna, db_pref, tblastn_bin, blastn_bin, work_dir, cmd_log
            )
            # Lookup orthogroup for ND90Pr protein
            og_match = og_df[og_df["Bsor_ND90Pr"].fillna("").str.contains(model["protein_id"], regex=False)]
            res["orthogroup_id"] = str(og_match["Orthogroup"].iloc[0]) if not og_match.empty else "Singleton"
            homolog_rows.append(res)

    homolog_df = pd.DataFrame(homolog_rows)
    homolog_p = ortho_dir / "target_homologs.tsv"
    homolog_df.drop(columns=["cds_seq", "protein_seq"]).to_csv(homolog_p, sep="\t", index=False)

    # 4. Run MAFFT on CDS and Protein per target & compute nucleotide conservation
    cds_aln_dir = ortho_dir / "target_cds_alignments"
    prot_aln_dir = ortho_dir / "target_protein_alignments"
    cds_aln_dir.mkdir(parents=True, exist_ok=True)
    prot_aln_dir.mkdir(parents=True, exist_ok=True)

    cons_summary_rows = []
    pairwise_all_rows = []
    per_pos_dict: Dict[str, pd.DataFrame] = {}

    for tname, model in models.items():
        sub = homolog_df[(homolog_df["target"] == tname) & (homolog_df["cds_length_nt"] >= 100)].copy()
        cds_in = cds_aln_dir / f"{tname}_isolates_cds.fasta"
        prot_in = prot_aln_dir / f"{tname}_isolates_prot.fasta"
        cds_recs = []
        prot_recs = []
        for _, hr in sub.iterrows():
            sid = f"{hr['isolate_id']}|{hr['assembly_accession']}"
            cds_recs.append(SeqRecord(Seq(hr["cds_seq"]), id=sid, description=hr["status"]))
            if hr["protein_seq"]:
                prot_recs.append(SeqRecord(Seq(hr["protein_seq"]), id=sid, description=hr["status"]))
        SeqIO.write(cds_recs, str(cds_in), "fasta")
        SeqIO.write(prot_recs, str(prot_in), "fasta")

        cds_aln_out = mafft_dir / f"{tname}_cds_aligned.fasta"
        prot_aln_out = mafft_dir / f"{tname}_prot_aligned.fasta"
        cds_msa = run_mafft(mafft_bin, cds_in, cds_aln_out, log_path=cmd_log)
        prot_msa = run_mafft(mafft_bin, prot_in, prot_aln_out, log_path=cmd_log)
        # Copy aligned files into target_cds_alignments / target_protein_alignments too
        (cds_aln_dir / f"{tname}_cds_aligned.fasta").write_text(cds_aln_out.read_text(encoding="utf-8"), encoding="utf-8")
        (prot_aln_dir / f"{tname}_prot_aligned.fasta").write_text(prot_aln_out.read_text(encoding="utf-8"), encoding="utf-8")

        ref_sid = f"ND90Pr|GCF_000338995.1"
        cons = analyze_msa_conservation(cds_msa, ref_id=ref_sid, window_size=21)
        per_pos_df = cons["per_position"]
        per_pos_df.to_csv(cons_dir / f"{tname}_per_position_conservation.tsv", sep="\t", index=False)
        per_pos_dict[tname] = per_pos_df

        for pw in cons["pairwise_stats"]:
            pairwise_all_rows.append({"target": tname, **pw})

        # Also check protein conservation
        prot_cons = analyze_msa_conservation(prot_msa, ref_id=ref_sid, window_size=10)

        blocks = cons["conserved_blocks_ge21nt"]
        longest_block = max((b[2] for b in blocks), default=0)
        blocks_ge150 = [b for b in blocks if b[2] >= 150]

        cons_summary_rows.append(
            {
                "target": tname,
                "orthogroup_id": sub["orthogroup_id"].iloc[0] if not sub.empty else "NA",
                "isolates_evaluated": len(df_iso),
                "isolates_with_gene": len(sub),
                "isolates_complete_cds": int((sub["status"].str.startswith("PRESENT_COMPLETE")).sum()),
                "ref_cds_length_nt": cons["ref_length_nt"],
                "mean_pairwise_nt_identity_pct": cons["mean_pairwise_identity_pct"],
                "min_pairwise_nt_identity_pct": cons["min_pairwise_identity_pct"],
                "mean_pairwise_aa_identity_pct": prot_cons["mean_pairwise_identity_pct"],
                "snp_sites": cons["snp_sites"],
                "indel_columns": cons["indel_columns"],
                "conserved_blocks_ge21nt_count": len(blocks),
                "conserved_blocks_ge150nt_count": len(blocks_ge150),
                "longest_100pct_conserved_block_nt": longest_block,
                "mean_21mer_isolate_conservation_pct": cons["mean_21mer_isolate_conservation_pct"],
                "pct_100pct_conserved_21mers": cons["pct_100pct_conserved_21mers"],
                "orthology_method": og_method,
            }
        )

    # 5. Also analyze CsPKS1 supplied partial CDS records (KF769134-KF769141 + HQ830033)
    pks1_sup_fa = cfg.output_dir / "01_data_acquisition" / "external_references" / "CsPKS1_supplied_records.fasta"
    if pks1_sup_fa.exists():
        pks1_sup_aln = mafft_dir / "CsPKS1_supplied_records_aligned.fasta"
        run_mafft(mafft_bin, pks1_sup_fa, pks1_sup_aln, log_path=cmd_log)

    cons_df = pd.DataFrame(cons_summary_rows)
    cons_p = cons_dir / "target_nucleotide_conservation.tsv"
    cons_df.to_csv(cons_p, sep="\t", index=False)
    pd.DataFrame(pairwise_all_rows).to_csv(cons_dir / "target_pairwise_nucleotide_identity.tsv", sep="\t", index=False)

    # Copy mandatory final report tables
    fin_dir = cfg.module_dir("10_final_reports")
    og_df.to_csv(fin_dir / "BRD_Orthogroups.tsv", sep="\t", index=False)
    homolog_df.drop(columns=["cds_seq", "protein_seq"]).to_csv(fin_dir / "BRD_Target_Orthologs.tsv", sep="\t", index=False)
    cons_df.to_csv(fin_dir / "BRD_Target_Nucleotide_Conservation.tsv", sep="\t", index=False)

    # 6. Plot Multi-Isolate Nucleotide Conservation Figures (PNG + SVG)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    axes_flat = axes.flatten()
    colors = {"CsPKS1": "#7c3aed", "CsFUS3": "#2563eb", "CsSTE11": "#059669", "CsSp1": "#dc2626"}
    for idx, tname in enumerate(["CsPKS1", "CsFUS3", "CsSTE11", "CsSp1"]):
        ax = axes_flat[idx]
        dfp = per_pos_dict.get(tname)
        if dfp is not None and not dfp.empty:
            ax.plot(dfp["ref_pos_1based"], dfp["conservation_frac"] * 100.0, color=colors[tname], linewidth=1.2)
            snps = dfp[dfp["is_snp"]]
            if not snps.empty:
                ax.scatter(snps["ref_pos_1based"], snps["conservation_frac"] * 100.0, color="#ef4444", s=18, zorder=5, label=f"SNP sites (n={len(snps)})")
            row_c = cons_df[cons_df["target"] == tname].iloc[0]
            ax.set_title(
                f"{tname} ({row_c['isolates_with_gene']}/{row_c['isolates_evaluated']} isolates; "
                f"mean id={row_c['mean_pairwise_nt_identity_pct']}%, longest 100% block={row_c['longest_100pct_conserved_block_nt']} nt)",
                fontsize=10,
                fontweight="bold",
            )
            ax.set_xlabel("Reference CDS Position (nt)")
            ax.set_ylabel("Nucleotide Conservation (%)")
            ax.set_ylim(50, 102)
            ax.grid(True, linestyle="--", alpha=0.35)
            ax.legend(loc="lower right", fontsize=8)

    fig.suptitle("Multi-Isolate CDS Nucleotide Conservation Across B. sorokiniana Isolates", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(fig_dir / "multi_isolate_nucleotide_conservation.png", dpi=300)
    fig.savefig(fig_dir / "multi_isolate_nucleotide_conservation.svg")
    plt.close(fig)

    prov.record("M04", "Comparative Genomics", str(homolog_p), cons_p, "DIAMOND+BLAST+MAFFT", "v2.2.8/2.17.0+/7.526", f"{len(df_iso)}_isolates")
    logger.info("M04 complete: nucleotide conservation computed across %d isolates for all 4 targets.", len(df_iso))
    return {
        "orthogroups": of_dir / "Orthogroups.tsv",
        "core_accessory": core_p,
        "target_homologs": homolog_p,
        "target_nucleotide_conservation": cons_p,
    }
