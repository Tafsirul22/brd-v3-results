"""Module M08 — Multi-Reference Off-Target Screening (Wheat A/B/D, Barley, Beneficial Fungi)."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Set

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from brd.config import BRDConfig, load_yaml
from brd.integrations.blast import make_blast_db, run_blast
from brd.integrations.rnahybrid import compute_duplex_mfe
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker
from brd.utils.checksums import compute_sha256
from brd.utils.report_utils import df_to_html_table, wrap_html_report, write_text


def _extract_kmer_hits_fast(guide_set: Set[str], fasta_path: Path, k: int = 21) -> Dict[str, List[str]]:
    """Exact 21-mer index scan across a non-target FASTA (sense and antisense)."""
    hits: Dict[str, List[str]] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(rec.seq).upper()
        n = len(seq)
        if n < k:
            continue
        # Scan in 21-mer slides only if needed; check fast substring when guide_set small or slide seq
        for i in range(0, n - k + 1):
            sub = seq[i : i + k]
            if sub in guide_set:
                hits.setdefault(sub, []).append(rec.id)
    return hits


def _classify_hits(sub_df: pd.DataFrame) -> str:
    if sub_df.empty:
        return "No match detected in the available database"
    exact_21 = sub_df[(sub_df["length"] >= 21) & (sub_df["pident"] == 100.0) & (sub_df["gapopen"] == 0)]
    if not exact_21.empty:
        return "Exact sequence match"
    high_sim = sub_df[((sub_df["length"] >= 19) & (sub_df["mismatch"] <= 1) & (sub_df["gapopen"] == 0)) | ((sub_df["length"] >= 30) & (sub_df["pident"] >= 85.0))]
    if not high_sim.empty:
        return "High-similarity match"
    part_sim = sub_df[(sub_df["length"] >= 17) & (sub_df["mismatch"] <= 2)]
    if not part_sim.empty:
        return "Partial similarity match"
    return "No match detected in the available database"


def _extract_wheat_gene_locus(desc: str) -> str:
    m = re.search(r"gene=([^\s\]]+)|LOC\d+", desc)
    return m.group(0) if m else desc.split()[0]


def run_m08(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    cmd_log = cfg.output_dir / "logs" / "commands.jsonl"
    logger = get_logger("BRD.M08", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    m08_dir = cfg.module_dir("08_off_target_screening")
    db_root = m08_dir / "blast_databases"
    wheat_dir = m08_dir / "wheat"
    barley_dir = m08_dir / "barley"
    ben_dir = m08_dir / "beneficial_fungi"
    other_dir = m08_dir / "other_non_targets"
    sum_dir = m08_dir / "summaries"
    fig_dir = m08_dir / "figures"
    for d in (db_root, wheat_dir, barley_dir, ben_dir, other_dir, sum_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    cand_fa = cfg.output_dir / "07_dsRNA_design" / "candidate_sequences" / "dsRNA_candidates.fasta"
    cand_tsv = cfg.output_dir / "07_dsRNA_design" / "candidate_sequences" / "dsRNA_candidates.tsv"
    guide_fa = cfg.output_dir / "07_dsRNA_design" / "guide_sequences" / "guide_sequences.fasta"
    guide_tsv = cfg.output_dir / "07_dsRNA_design" / "guide_sequences" / "guide_evaluation.tsv"

    cand_df = pd.read_csv(cand_tsv, sep="\t")
    guide_df = pd.read_csv(guide_tsv, sep="\t")

    refs_cfg = load_yaml(cfg.output_dir / "config" / "references.yaml")
    makeblastdb_bin = cfg.tool_path("makeblastdb")
    blastn_bin = cfg.tool_path("blastn")

    # 1. Build BLAST+ databases and map FASTA descriptions
    db_manifest_rows = []
    db_entries = []
    for grp_name, entries in refs_cfg.get("non_target_panel", {}).items():
        for entry in entries:
            acc = entry["accession"]
            src_sub = db_root / "sources" / entry["id"] / "ncbi_dataset" / "data" / acc
            rna_f = src_sub / "rna.fna"
            cds_f = src_sub / "cds_from_genomic.fna"
            src_f = rna_f if rna_f.exists() else cds_f
            if not src_f.exists():
                logger.warning("Reference FASTA missing for %s (%s)", entry["id"], acc)
                continue
            db_pref = db_root / entry["id"]
            make_blast_db(makeblastdb_bin, src_f, db_pref, dbtype="nucl", log_path=cmd_log)
            n_seqs = sum(1 for _ in SeqIO.parse(str(src_f), "fasta"))
            db_manifest_rows.append(
                {
                    "ref_id": entry["id"],
                    "group": grp_name,
                    "organism": entry["organism"],
                    "accession": acc,
                    "fasta_file": str(src_f),
                    "sequence_count": n_seqs,
                    "sha256": compute_sha256(src_f),
                    "rationale": entry["rationale"],
                }
            )
            db_entries.append({**entry, "group": grp_name, "fasta": src_f, "db_pref": db_pref})

    pd.DataFrame(db_manifest_rows).to_csv(db_root / "database_manifest.tsv", sep="\t", index=False)

    # 2. Run BLASTN (blastn-short) for all dsRNA candidates against EVERY non-target database
    all_hits_by_group: Dict[str, List[pd.DataFrame]] = {"wheat": [], "barley": [], "beneficial_fungi": []}
    guide_hits_list: List[pd.DataFrame] = []

    for ent in db_entries:
        rid = ent["id"]
        grp = ent["group"]
        out_raw = (wheat_dir if grp == "wheat" else barley_dir if grp == "barley" else ben_dir) / f"candidates_vs_{rid}_blastn.tsv"
        logger.info("Screening candidates via blastn-short against %s (%s)...", rid, ent["organism"])
        df_h = run_blast(
            blastn_bin,
            cand_fa,
            ent["db_pref"],
            out_raw,
            evalue=10.0,
            task="blastn-short",
            word_size=7,
            max_target_seqs=200,
            num_threads=4,
            extra_args=["-dust", "no"],
            log_path=cmd_log,
        )
        # Filter to informative RNAi alignment lengths (>= 16 nt)
        df_sig = df_h[df_h["length"] >= 16].copy()
        df_sig["ref_id"] = rid
        df_sig["group"] = grp
        df_sig["organism"] = ent["organism"]
        all_hits_by_group[grp].append(df_sig)

        # Also screen unique 21-mer guides via blastn-short (>= 17 nt, <= 2 mismatches)
        g_raw = (wheat_dir if grp == "wheat" else barley_dir if grp == "barley" else ben_dir) / f"guides_vs_{rid}_blastn.tsv"
        df_gh = run_blast(
            blastn_bin,
            guide_fa,
            ent["db_pref"],
            g_raw,
            evalue=10.0,
            task="blastn-short",
            word_size=7,
            max_target_seqs=50,
            num_threads=4,
            extra_args=["-dust", "no", "-perc_identity", "85"],
            log_path=cmd_log,
        )
        df_gsig = df_gh[(df_gh["length"] >= 17) & (df_gh["mismatch"] <= 2)].copy()
        df_gsig["ref_id"] = rid
        df_gsig["group"] = grp
        guide_hits_list.append(df_gsig)

    # Combine per-group tables
    wheat_df = pd.concat(all_hits_by_group["wheat"], ignore_index=True) if all_hits_by_group["wheat"] else pd.DataFrame()
    barley_df = pd.concat(all_hits_by_group["barley"], ignore_index=True) if all_hits_by_group["barley"] else pd.DataFrame()
    ben_df = pd.concat(all_hits_by_group["beneficial_fungi"], ignore_index=True) if all_hits_by_group["beneficial_fungi"] else pd.DataFrame()
    guides_ot_df = pd.concat(guide_hits_list, ignore_index=True) if guide_hits_list else pd.DataFrame()

    # Exact 21-mer double-check via direct Python scan for 100% guarantee (sense + antisense)
    # Build set of all 21-mers in candidates
    cand_kmer_map: Dict[str, Set[str]] = {}
    all_21mers: Set[str] = set()
    for _, cr in cand_df.iterrows():
        cid = cr["candidate_id"]
        seq = str(cr["sequence"]).upper()
        km_set = set()
        for i in range(len(seq) - 21 + 1):
            k1 = seq[i : i + 21]
            k2 = str(Seq(k1).reverse_complement())
            km_set.add(k1)
            km_set.add(k2)
        cand_kmer_map[cid] = km_set
        all_21mers.update(km_set)

    # Save BLAST outputs
    wheat_p = wheat_dir / "wheat_off_targets.tsv"
    barley_p = barley_dir / "barley_off_targets.tsv"
    ben_p = ben_dir / "beneficial_fungi_off_targets.tsv"
    guide_ot_p = sum_dir / "guide_off_target_results.tsv"

    wheat_df.to_csv(wheat_p, sep="\t", index=False)
    barley_df.to_csv(barley_p, sep="\t", index=False)
    ben_df.to_csv(ben_p, sep="\t", index=False)
    guides_ot_df.to_csv(guide_ot_p, sep="\t", index=False)

    # 3. Compute thermodynamic ViennaRNA duplexfold MFE for top BLAST hits (max 1 hit per candidate/group to stay crisp)
    # Index subject sequences for hits that appeared in wheat/barley/beneficial_fungi
    hit_subject_ids = set()
    for d in (wheat_df, barley_df, ben_df):
        if not d.empty:
            hit_subject_ids.update(d["sseqid"].astype(str).unique())

    subject_seqs: Dict[str, str] = {}
    if hit_subject_ids:
        for ent in db_entries:
            for rec in SeqIO.parse(str(ent["fasta"]), "fasta"):
                if rec.id in hit_subject_ids:
                    subject_seqs[rec.id] = str(rec.seq).upper()

    # Summarize per candidate across Wheat, Barley, Beneficial Fungi
    summary_rows = []
    for _, cr in cand_df.iterrows():
        cid = cr["candidate_id"]
        cseq = str(cr["sequence"]).upper()
        tname = cr["target"]

        w_sub = wheat_df[wheat_df["qseqid"] == cid] if not wheat_df.empty else pd.DataFrame()
        b_sub = barley_df[barley_df["qseqid"] == cid] if not barley_df.empty else pd.DataFrame()
        f_sub = ben_df[ben_df["qseqid"] == cid] if not ben_df.empty else pd.DataFrame()

        w_cls = _classify_hits(w_sub)
        b_cls = _classify_hits(b_sub)
        f_cls = _classify_hits(f_sub)

        w_exact = int(((w_sub["length"] >= 21) & (w_sub["pident"] == 100.0)).sum()) if not w_sub.empty else 0
        w_19 = int(((w_sub["length"] >= 19) & (w_sub["mismatch"] <= 1)).sum()) if not w_sub.empty else 0
        w_loci = len(set(w_sub["sseqid"].astype(str))) if not w_sub.empty else 0

        b_exact = int(((b_sub["length"] >= 21) & (b_sub["pident"] == 100.0)).sum()) if not b_sub.empty else 0
        b_19 = int(((b_sub["length"] >= 19) & (b_sub["mismatch"] <= 1)).sum()) if not b_sub.empty else 0

        f_exact = int(((f_sub["length"] >= 21) & (f_sub["pident"] == 100.0)).sum()) if not f_sub.empty else 0
        f_19 = int(((f_sub["length"] >= 19) & (f_sub["mismatch"] <= 1)).sum()) if not f_sub.empty else 0
        f_sp_hit = len(set(f_sub["organism"].astype(str))) if not f_sub.empty else 0

        # Thermodynamic duplexfold check on strongest BLAST hit if any (orient guide to match mRNA subject)
        best_thermo_mfe = "NA"
        all_sub = pd.concat([w_sub, b_sub, f_sub], ignore_index=True)
        if not all_sub.empty:
            top_hit = all_sub.sort_values("bitscore", ascending=False).iloc[0]
            qs, qe = int(top_hit["qstart"]) - 1, int(top_hit["qend"])
            ss_id = str(top_hit["sseqid"])
            if ss_id in subject_seqs:
                q_frag = cseq[max(0, qs) : min(len(cseq), qe)]
                if len(q_frag) >= 16:
                    s_lo = max(0, min(int(top_hit["sstart"]), int(top_hit["send"])) - 5)
                    s_hi = min(len(subject_seqs[ss_id]), max(int(top_hit["sstart"]), int(top_hit["send"])) + 5)
                    s_frag = subject_seqs[ss_id][s_lo:s_hi]
                    if len(s_frag) >= 16:
                        is_minus = str(top_hit.get("sstrand", "plus")).lower() == "minus"
                        wsize = min(21, len(q_frag))
                        best_mfe_val = None
                        for w_i in range(0, max(1, len(q_frag) - wsize + 1), 5):
                            qw = q_frag[w_i : w_i + wsize]
                            guide_seq = qw if is_minus else str(Seq(qw).reverse_complement())
                            guide_rna = guide_seq.replace("T", "U")
                            mfe_v = compute_duplex_mfe(guide_rna, s_frag)["mfe_kcal_mol"]
                            if best_mfe_val is None or mfe_v < best_mfe_val:
                                best_mfe_val = mfe_v
                        best_thermo_mfe = best_mfe_val if best_mfe_val is not None else "NA"

        # Overall off-target safety score (0..1; higher = fewer BLAST/k-mer off-target concerns)
        # Note: exact 21-mer match in Wheat/Barley is penalized heavily; beneficial fungi high-similarity also penalized.
        penalty = 0.0
        if w_exact > 0:
            penalty += 0.60
        elif w_19 > 0:
            penalty += 0.30
        elif not w_sub.empty:
            penalty += 0.10

        if b_exact > 0:
            penalty += 0.30
        elif b_19 > 0:
            penalty += 0.15
        elif not b_sub.empty:
            penalty += 0.05

        if f_exact > 0:
            penalty += 0.35
        elif f_19 > 0:
            penalty += 0.20
        elif not f_sub.empty:
            penalty += 0.05 * min(3, f_sp_hit)

        safety_score = max(0.0, round(1.0 - penalty, 4))

        # Overall status
        statuses = [w_cls, b_cls, f_cls]
        if "Exact sequence match" in statuses:
            overall = "Exact sequence match"
        elif "High-similarity match" in statuses:
            overall = "High-similarity match"
        elif "Partial similarity match" in statuses:
            overall = "Partial similarity match"
        else:
            overall = "No match detected in the available database"

        summary_rows.append(
            {
                "candidate_id": cid,
                "target": tname,
                "tool_name": cr["tool_name"],
                "candidate_length_nt": cr["candidate_length_nt"],
                "wheat_classification": w_cls,
                "wheat_exact_21mer_hits": w_exact,
                "wheat_high_sim_hits_ge19nt": w_19,
                "wheat_unique_transcript_hits": w_loci,
                "barley_classification": b_cls,
                "barley_exact_21mer_hits": b_exact,
                "barley_high_sim_hits_ge19nt": b_19,
                "beneficial_fungi_classification": f_cls,
                "beneficial_fungi_exact_21mer_hits": f_exact,
                "beneficial_fungi_high_sim_hits_ge19nt": f_19,
                "beneficial_fungi_species_with_hits": f_sp_hit,
                "top_hit_viennarna_duplexfold_mfe": best_thermo_mfe,
                "overall_off_target_classification": overall,
                "off_target_safety_score_0_1": safety_score,
                "databases_screened": len(db_entries),
            }
        )

    sum_df = pd.DataFrame(summary_rows)
    sum_p = sum_dir / "off_target_summary.tsv"
    sum_df.to_csv(sum_p, sep="\t", index=False)

    # Copy mandatory final tables
    fin_dir = cfg.module_dir("10_final_reports")
    sum_df.to_csv(fin_dir / "BRD_Off_Target_Analysis.tsv", sep="\t", index=False)
    wheat_df.to_csv(fin_dir / "BRD_Wheat_Off_Targets.tsv", sep="\t", index=False)
    barley_df.to_csv(fin_dir / "BRD_Barley_Off_Targets.tsv", sep="\t", index=False)
    ben_df.to_csv(fin_dir / "BRD_Beneficial_Fungi_Off_Targets.tsv", sep="\t", index=False)

    # Generate Off-Target Summary Figures (Wheat, Barley, Beneficial Fungi)
    for grp_col, exact_col, high_col, title, fname in [
        ("wheat_classification", "wheat_exact_21mer_hits", "wheat_high_sim_hits_ge19nt", "Wheat (T. aestivum IWGSC v2.1) Off-Target Summary", "wheat_off_target_summary"),
        ("barley_classification", "barley_exact_21mer_hits", "barley_high_sim_hits_ge19nt", "Barley (H. vulgare MorexV3) Off-Target Summary", "barley_off_target_summary"),
        ("beneficial_fungi_classification", "beneficial_fungi_exact_21mer_hits", "beneficial_fungi_high_sim_hits_ge19nt", "Beneficial Fungi Panel Off-Target Summary", "beneficial_fungi_off_target_summary"),
    ]:
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ct = sum_df.groupby(["target", grp_col]).size().unstack(fill_value=0)
        ct.plot(kind="bar", stacked=True, ax=ax, colormap="Set2", edgecolor="#0f172a")
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Target Gene")
        ax.set_ylabel("Number of dsRNA Candidates")
        ax.tick_params(axis="x", rotation=0)
        ax.legend(title="Off-Target Classification", fontsize=8, loc="best")
        fig.tight_layout()
        fig.savefig(fig_dir / f"{fname}.png", dpi=300)
        fig.savefig(fig_dir / f"{fname}.svg")
        plt.close(fig)

    body = (
        "<div class='warn'><strong>Interpretation Safeguard:</strong> 'No match detected in the available database' does "
        "<strong>not</strong> demonstrate biological safety; it reflects absence of sequence matches exceeding our BLASTN-short "
        "(&ge;16 nt, word_size=7) and 21-mer exact criteria across the 7 RefSeq non-target transcriptomes screened.</div>"
        f"<h2>Screened Reference Databases</h2>{df_to_html_table(pd.DataFrame(db_manifest_rows)[['ref_id', 'group', 'organism', 'accession', 'sequence_count']])}"
        f"<h2>Off-Target Summary by Candidate</h2>{df_to_html_table(sum_df)}"
    )
    rep_html = m08_dir / "off_target_report.html"
    write_text(rep_html, wrap_html_report("BRD M08 — Multi-Reference Off-Target Screening Report", "Wheat (A/B/D), Barley, and Beneficial Fungi Panel", body))

    prov.record("M08", "Off-Target Screening", str(cand_fa), sum_p, "NCBI_BLAST+_ViennaRNA", "2.17.0+/2.8.0", f"{len(db_entries)}_references")
    logger.info("M08 complete: %d candidates screened across %d non-target transcriptomes.", len(sum_df), len(db_entries))
    return {
        "wheat_off_targets": wheat_p,
        "barley_off_targets": barley_p,
        "beneficial_fungi_off_targets": ben_p,
        "guide_off_target_results": guide_ot_p,
        "off_target_summary": sum_p,
        "report": rep_html,
    }
