"""Module M07 — dsRNA Candidate Design (LUCID/sifi21_plus, SIREN, dsRNAmax, and BRD Isolate-Aware Design)."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from brd.config import BRDConfig
from brd.integrations.blast import make_blast_db, run_blast
from brd.integrations.dsRNAmax import run_dsrnamax
from brd.integrations.lucid import design_lucid_windows, score_sirnas_lucid_sifi21
from brd.integrations.siren import evaluate_sirna_duplexfold_offtargets, run_siren_design_cli, run_siren_prefilter_cli
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker
from brd.utils.file_validation import validate_nucleotide_seq
from brd.utils.report_utils import df_to_html_table, wrap_html_report, write_text


def _max_homopolymer(seq: str) -> int:
    best = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return best


def _extract_isolate_cds_map(cfg: BRDConfig, target_name: str) -> Dict[str, str]:
    aln_fa = cfg.output_dir / "04_comparative_genomics" / "target_orthologs" / "target_cds_alignments" / f"{target_name}_isolates_cds.fasta"
    out = {}
    for rec in SeqIO.parse(str(aln_fa), "fasta"):
        out[rec.id] = str(rec.seq).upper().replace("-", "")
    return out


def _compute_candidate_isolate_conservation(cand_seq: str, iso_cds_map: Dict[str, str], k: int = 21) -> Dict[str, object]:
    """Compute exact 21-mer isolate coverage and per-isolate conservation vector for a dsRNA candidate."""
    seq = cand_seq.upper()
    kmers = [seq[i : i + k] for i in range(0, len(seq) - k + 1) if "N" not in seq[i : i + k]]
    total_k = max(len(kmers), 1)
    per_iso = {}
    iso_covered_ge90 = 0
    iso_covered_100 = 0
    for iso_id, iso_seq in iso_cds_map.items():
        hits = sum(1 for km in kmers if km in iso_seq)
        frac = hits / total_k
        per_iso[iso_id] = round(frac, 4)
        if frac >= 0.90:
            iso_covered_ge90 += 1
        if frac == 1.0:
            iso_covered_100 += 1
    n_iso = max(len(iso_cds_map), 1)
    mean_kmer_cov = sum(per_iso.values()) / n_iso * 100.0
    min_kmer_cov = min(per_iso.values()) * 100.0 if per_iso else 0.0
    return {
        "n_isolates_evaluated": len(iso_cds_map),
        "isolates_ge90pct_21mer": iso_covered_ge90,
        "isolates_100pct_21mer": iso_covered_100,
        "mean_21mer_conservation_pct": round(mean_kmer_cov, 2),
        "min_21mer_conservation_pct": round(min_kmer_cov, 2),
        "per_isolate_fractions": per_iso,
    }


def _build_target_self_tx_ids(cfg: BRDConfig) -> Dict[str, set]:
    """Map each target name to its own ND90Pr locus_tag and RefSeq mRNA accession (e.g. XM_007698023.1)."""
    out: Dict[str, set] = {t: {tcfg["nd90pr_locus_tag"]} for t, tcfg in cfg.raw["targets"].items()}
    iso_p = cfg.output_dir / "01_data_acquisition" / "isolate_metadata.tsv"
    if not iso_p.exists():
        return out
    meta = pd.read_csv(iso_p, sep="\t")
    ref_rows = meta[meta["assembly_accession"] == "GCF_000338995.1"]
    if ref_rows.empty:
        return out
    gff = Path(str(ref_rows.iloc[0]["gff_path"]))
    if not gff.exists():
        return out
    pid_to_target = {tcfg["nd90pr_protein_id"]: t for t, tcfg in cfg.raw["targets"].items()}
    with open(gff, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            m_p = re.search(r"Name=([^;]+)", parts[8])
            if not m_p or m_p.group(1) not in pid_to_target:
                continue
            tname = pid_to_target[m_p.group(1)]
            m_par = re.search(r"Parent=rna-([^;]+)", parts[8])
            if m_par:
                out[tname].add(m_par.group(1))
    return out


def _index_kmer_hits_in_transcriptome(query_kmers: set, bsor_rnas: Dict[str, str], k: int = 21) -> Dict[str, set]:
    """Single-pass O(N_transcriptome) exact 21-mer index for all candidate k-mers (both strands)."""
    hits: Dict[str, set] = {}
    if not query_kmers:
        return hits
    for rid, rseq in bsor_rnas.items():
        n = len(rseq)
        if n < k:
            continue
        for i in range(0, n - k + 1):
            sub = rseq[i : i + k]
            if sub in query_kmers:
                hits.setdefault(sub, set()).add(rid)
    return hits


def _check_self_paralog_21mers(
    cand_seq: str,
    self_ids: set,
    kmer_index: Dict[str, set],
    k: int = 21,
) -> int:
    """Count OTHER B. sorokiniana transcripts sharing at least one 21-mer (sense or antisense)."""
    matched_rids = set()
    for i in range(len(cand_seq) - k + 1):
        km = cand_seq[i : i + k]
        rc = str(Seq(km).reverse_complement())
        for q in (km, rc):
            rids = kmer_index.get(q)
            if rids:
                for rid in rids:
                    if not any(sid in rid for sid in self_ids):
                        matched_rids.add(rid)
    return len(matched_rids)


def run_m07(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    cmd_log = cfg.output_dir / "logs" / "commands.jsonl"
    logger = get_logger("BRD.M07", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    m07_dir = cfg.module_dir("07_dsRNA_design")
    lucid_dir = m07_dir / "lucid"
    siren_dir = m07_dir / "siren"
    dsmax_dir = m07_dir / "dsRNAmax"
    cand_dir = m07_dir / "candidate_sequences"
    guide_dir = m07_dir / "guide_sequences"
    comp_dir = m07_dir / "candidate_comparison"
    fig_dir = m07_dir / "figures"
    for d in (lucid_dir, siren_dir, dsmax_dir, cand_dir, guide_dir, comp_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    python_bin = cfg.tool_path("python")
    siren_repo = Path(cfg.tool_path("siren_repo"))
    dsmax_bin = cfg.tool_path("dsrnamax")

    # Load ND90Pr reference transcriptome and exact GFF transcript IDs for self/paralog screening
    ref_rna_fa = cfg.output_dir / "03_annotation" / "transcripts" / "genome_transcripts.fasta"
    bsor_rnas = {r.description: str(r.seq).upper() for r in SeqIO.parse(str(ref_rna_fa), "fasta")}
    target_self_tx_ids = _build_target_self_tx_ids(cfg)

    target_cds_map: Dict[str, str] = {}
    target_meta_map: Dict[str, dict] = {}
    for tname, tcfg in cfg.raw["targets"].items():
        cds_fa = cfg.output_dir / "06_target_prioritisation" / "target_sequences" / f"{tname}_cds.fasta"
        rec = list(SeqIO.parse(str(cds_fa), "fasta"))[0]
        target_cds_map[tname] = str(rec.seq).upper()
        target_meta_map[tname] = tcfg

    # Build combined target + B. sorokiniana transcript FASTA for SIREN prefilter
    siren_targets_fa = siren_dir / "bsor_targets_plus_transcriptome.fasta"
    with open(siren_targets_fa, "w", encoding="utf-8") as fh:
        for tname, seq in target_cds_map.items():
            fh.write(f">{tname}_CDS\n{seq}\n")
        for r in SeqIO.parse(str(ref_rna_fa), "fasta"):
            fh.write(f">{r.id} {r.description}\n{str(r.seq).upper()}\n")

    raw_candidates: List[dict] = []
    all_guides_rows: List[dict] = []
    all_guide_recs: List[SeqRecord] = []
    iso_cons_matrix_rows: List[dict] = []
    tool_comp_rows: List[dict] = []

    for tname, cds_seq in target_cds_map.items():
        logger.info("Designing dsRNA candidates for %s (CDS len = %d nt)...", tname, len(cds_seq))
        iso_cds_map = _extract_isolate_cds_map(cfg, tname)
        locus_tag = target_meta_map[tname]["nd90pr_locus_tag"]

        # --- TOOL 1: LUCID (sifi21_plus + ViennaRNA accessibility + Primer3, 200-400 nt baseline) ---
        lucid_df, sirna_df = design_lucid_windows(tname, cds_seq, min_len=200, max_len=400, step_len=50, step_pos=25)
        sirna_df.to_csv(lucid_dir / f"{tname}_sifi21_sirnas.tsv", sep="\t", index=False)
        lucid_df.to_csv(lucid_dir / f"{tname}_lucid_windows.tsv", sep="\t", index=False)

        # Select top 4 non-overlapping/distinct LUCID windows by efficient_sirna_fraction + accessibility
        lucid_ranked = lucid_df.sort_values(["efficient_sirna_fraction", "mean_accessibility"], ascending=False)
        for idx_l, (_, lr) in enumerate(lucid_ranked.head(4).iterrows(), 1):
            raw_candidates.append(
                {
                    "target": tname,
                    "tool": "LUCID_sifi21_plus",
                    "tool_version": "0.1.0+ViennaRNA_2.8.0",
                    "tool_score": float(lr["efficient_sirna_fraction"]),
                    "tool_score_metric": "efficient_sirna_fraction",
                    "start_1based": int(lr["start_1based"]),
                    "end_1based": int(lr["end_1based"]),
                    "sequence": lr["sequence"],
                    "mean_accessibility": float(lr["mean_accessibility"]),
                    "efficient_sirna_frac": float(lr["efficient_sirna_fraction"]),
                    "primer_left": lr["primer_left"],
                    "primer_right": lr["primer_right"],
                    "amplicon_bp": lr["expected_amplicon_bp"],
                }
            )

        # --- TOOL 2: SIREN (siren_prefilter + ViennaRNA duplexfold + siren_designVIII CLI) ---
        t_siren_dir = siren_dir / tname
        t_siren_dir.mkdir(parents=True, exist_ok=True)
        pref_fa = run_siren_prefilter_cli(python_bin, siren_repo, siren_targets_fa, f"{tname}_CDS", t_siren_dir, log_path=cmd_log)
        ot_sum_df = evaluate_sirna_duplexfold_offtargets(f"{tname}_CDS", cds_seq, pref_fa, sirna_size=21, step=2, mfe_threshold=-25.0)
        ot_sum_p = t_siren_dir / "off_targets_summary.tsv"
        ot_sum_df.to_csv(ot_sum_p, sep="\t", index=False)

        single_fa = t_siren_dir / f"{tname}_CDS.fasta"
        single_fa.write_text(f">{tname}_CDS\n{cds_seq}\n", encoding="utf-8")
        siren_df = run_siren_design_cli(
            python_bin,
            siren_repo,
            tname,
            single_fa,
            ot_sum_p,
            t_siren_dir / f"{tname}_siren_design.tsv",
            base_rnai_length=250,
            threads=2,
            log_path=cmd_log,
        )
        if not siren_df.empty:
            for _, sr in siren_df.head(4).iterrows():
                m_pos = re.search(r"_(\d+)-(\d+)_len", str(sr["Name"]))
                s_pos = int(m_pos.group(1)) if m_pos else 1
                e_pos = int(m_pos.group(2)) if m_pos else len(sr["Sequence"])
                raw_candidates.append(
                    {
                        "target": tname,
                        "tool": "SIREN",
                        "tool_version": "0.1.6+ViennaRNA_duplexfold",
                        "tool_score": float(sr["Score"]),
                        "tool_score_metric": "SIREN_cumulative_offtarget_penalty (0=clean)",
                        "start_1based": s_pos,
                        "end_1based": e_pos,
                        "sequence": str(sr["Sequence"]).upper(),
                        "mean_accessibility": float("nan"),
                        "efficient_sirna_frac": float("nan"),
                        "primer_left": str(sr.get("Primer_Left", "NA")),
                        "primer_right": str(sr.get("Primer_Right", "NA")),
                        "amplicon_bp": str(sr.get("Expected_Amplicon_Size", "")),
                    }
                )

        # --- TOOL 3: dsRNAmax v1.1.14 (multi-isolate pan-CDS k-mer maximisation) ---
        t_ds_dir = dsmax_dir / tname
        t_ds_dir.mkdir(parents=True, exist_ok=True)
        pan_fa = t_ds_dir / f"{tname}_pan_isolates.fasta"
        with open(pan_fa, "w", encoding="utf-8") as fh:
            for iso_id, iseq in iso_cds_map.items():
                if len(iseq) >= 150:
                    fh.write(f">{iso_id}\n{iseq}\n")
        # Construct length 300 nt (or 200 nt for CsSp1 which is 270 nt)
        ds_len = 300 if len(cds_seq) >= 350 else 200
        # Use prefiltered self-paralog off-target file (excluding target_name and its own RefSeq mRNA) for dsRNAmax
        self_ids_t = target_self_tx_ids.get(tname, {locus_tag})
        ot_ds_fa = t_ds_dir / f"{tname}_paralog_offtargets.fasta"
        with open(ot_ds_fa, "w", encoding="utf-8") as fh:
            for rec in SeqIO.parse(str(pref_fa), "fasta"):
                hdr = f"{rec.id} {rec.description}"
                if tname not in rec.id and not any(sid in hdr for sid in self_ids_t):
                    fh.write(f">{rec.id}\n{str(rec.seq).upper()}\n")
        ot_arg = [ot_ds_fa] if ot_ds_fa.stat().st_size > 0 else None
        ds_res = run_dsrnamax(
            dsmax_bin,
            pan_fa,
            t_ds_dir / f"{tname}_dsrnamax.csv",
            t_ds_dir / f"{tname}_dsrnamax.log",
            off_targets=ot_arg,
            construct_len=ds_len,
            kmer_len=21,
            ot_kmer_len=21,
            iterations=100,
            log_path=cmd_log,
        )
        if ds_res["sense_arm_sequence"]:
            # Locate best alignment span on ND90Pr CDS if contiguous
            s_arm = ds_res["sense_arm_sequence"].upper()
            pos_idx = cds_seq.find(s_arm[:30])
            raw_candidates.append(
                {
                    "target": tname,
                    "tool": "dsRNAmax",
                    "tool_version": "1.1.14",
                    "tool_score": float(ds_len - 20),
                    "tool_score_metric": "21nt_kmer_coverage_max",
                    "start_1based": (pos_idx + 1) if pos_idx != -1 else 1,
                    "end_1based": (pos_idx + len(s_arm)) if pos_idx != -1 else len(s_arm),
                    "sequence": s_arm,
                    "mean_accessibility": float("nan"),
                    "efficient_sirna_frac": float("nan"),
                    "primer_left": "NA (synthetic k-mer construct)",
                    "primer_right": "NA (synthetic k-mer construct)",
                    "amplicon_bp": "",
                }
            )

        # --- TOOL 4: BRD Isolate-Aware Multi-Length Sliding Window Scan (150-550 nt) ---
        cand_lengths = [l for l in cfg.raw["dsrna_design"]["candidate_lengths_nt"] if l <= len(cds_seq)]
        if not cand_lengths and len(cds_seq) >= 150:
            cand_lengths = [len(cds_seq)]
        brd_wins = []
        for wlen in cand_lengths:
            for start in range(0, len(cds_seq) - wlen + 1, 25):
                subseq = cds_seq[start : start + wlen]
                icons = _compute_candidate_isolate_conservation(subseq, iso_cds_map, k=21)
                sub_sifi = sirna_df[(sirna_df["sirna_pos_1based"] >= start + 1) & (sirna_df["sirna_pos_1based"] + 20 <= start + wlen)]
                eff_frac = float(sub_sifi["is_efficient_sifi21"].mean()) if not sub_sifi.empty else 0.0
                mean_acc = float(sub_sifi["accessibility_prob"].mean()) if not sub_sifi.empty else 0.0
                brd_wins.append(
                    {
                        "start": start + 1,
                        "end": start + wlen,
                        "wlen": wlen,
                        "seq": subseq,
                        "cons_pct": icons["mean_21mer_conservation_pct"],
                        "eff_frac": eff_frac,
                        "mean_acc": mean_acc,
                        "score": icons["mean_21mer_conservation_pct"] / 100.0 * 0.6 + eff_frac * 0.25 + mean_acc * 0.15,
                    }
                )
        # Keep top 5 across distinct lengths/positions
        brd_df = pd.DataFrame(brd_wins).sort_values("score", ascending=False)
        chosen_starts = []
        for _, br in brd_df.iterrows():
            if len(chosen_starts) >= 5:
                break
            if any(abs(br["start"] - cs) < 60 and br["wlen"] == cw for cs, cw in chosen_starts):
                continue
            chosen_starts.append((br["start"], br["wlen"]))
            raw_candidates.append(
                {
                    "target": tname,
                    "tool": "BRD_IsolateAware",
                    "tool_version": "1.0.0",
                    "tool_score": round(float(br["score"]), 4),
                    "tool_score_metric": "BRD_conservation_accessibility_sifi_composite",
                    "start_1based": int(br["start"]),
                    "end_1based": int(br["end"]),
                    "sequence": br["seq"],
                    "mean_accessibility": round(float(br["mean_acc"]), 4),
                    "efficient_sirna_frac": round(float(br["eff_frac"]), 4),
                    "primer_left": "NA",
                    "primer_right": "NA",
                    "amplicon_bp": "",
                }
            )

        # Tool comparison metrics per target
        tool_comp_rows.append(
            {
                "target": tname,
                "cds_length_nt": len(cds_seq),
                "lucid_windows_evaluated": len(lucid_df),
                "siren_candidates_evaluated": len(siren_df),
                "dsrnamax_construct_len": ds_res["construct_len"],
                "brd_windows_evaluated": len(brd_df),
                "lucid_note": "sifi21_plus + ViennaRNA accessibility (200-400 nt)" if len(cds_seq) >= 400 else f"CDS {len(cds_seq)} nt < 400 nt max; bounded to <= {len(cds_seq)} nt",
                "siren_note": "siren_prefilter + ViennaRNA duplexfold + siren_designVIII",
                "dsrnamax_note": f"dsRNAmax v1.1.14 pan-isolate k-mer assembly ({ds_len} nt)",
            }
        )

    # Also run optional chimeric multi-target dsRNAmax (CsFUS3 + CsSTE11)
    chim_dir = dsmax_dir / "chimeric_CsFUS3_CsSTE11"
    chim_dir.mkdir(parents=True, exist_ok=True)
    chim_in = chim_dir / "CsFUS3_CsSTE11_targets.fasta"
    chim_in.write_text(f">CsFUS3_CDS\n{target_cds_map['CsFUS3']}\n>CsSTE11_CDS\n{target_cds_map['CsSTE11']}\n", encoding="utf-8")
    chim_res = run_dsrnamax(
        dsmax_bin,
        chim_in,
        chim_dir / "chimeric_dsrnamax.csv",
        chim_dir / "chimeric_dsrnamax.log",
        construct_len=350,
        kmer_len=21,
        ot_kmer_len=21,
        iterations=150,
        log_path=cmd_log,
    )

    # Deduplicate identical sequences while preserving tool attribution & assign stable BRD candidate IDs
    seen_seqs: Dict[str, str] = {}
    final_cand_rows: List[dict] = []
    final_cand_recs: List[SeqRecord] = []
    guide_seen = set()

    # Load all target sirna_dfs for guide property lookup
    sifi_lookup: Dict[str, dict] = {}
    for tname in target_cds_map:
        sdf = pd.read_csv(lucid_dir / f"{tname}_sifi21_sirnas.tsv", sep="\t")
        for _, sr in sdf.iterrows():
            sifi_lookup[f"{tname}:{sr['sirna_seq']}"] = sr.to_dict()

    # Single-pass O(N_transcriptome) index across all candidate 21-mers for fast paralog counting
    iso_cds_cache = {tname: _extract_isolate_cds_map(cfg, tname) for tname in target_cds_map}
    all_cand_kmers = set()
    for rc in raw_candidates:
        s = rc["sequence"].upper()
        for i in range(max(0, len(s) - 20)):
            km = s[i : i + 21]
            all_cand_kmers.add(km)
            all_cand_kmers.add(str(Seq(km).reverse_complement()))
    paralog_kmer_index = _index_kmer_hits_in_transcriptome(all_cand_kmers, bsor_rnas, k=21)

    cand_counter: Dict[str, int] = {}
    for rc in raw_candidates:
        tname = rc["target"]
        seq = rc["sequence"].upper()
        ok_seq, seq_stats = validate_nucleotide_seq(seq)
        if not ok_seq or len(seq) < 150:
            continue
        if seq in seen_seqs:
            continue
        cand_counter[tname] = cand_counter.get(tname, 0) + 1
        cid = f"BRD_{tname}_CAND{cand_counter[tname]:03d}"
        seen_seqs[seq] = cid

        cds_ref = target_cds_map[tname]
        iso_cds_map = iso_cds_cache[tname]
        icons = _compute_candidate_isolate_conservation(seq, iso_cds_map, k=21)
        locus_tag = target_meta_map[tname]["nd90pr_locus_tag"]
        paralog_hits = _check_self_paralog_21mers(seq, target_self_tx_ids.get(tname, {locus_tag}), paralog_kmer_index, k=21)
        max_hp = _max_homopolymer(seq)

        # Verify against reference CDS
        exact_in_ref = seq in cds_ref
        ref_id_pct = 100.0 if exact_in_ref else icons["per_isolate_fractions"].get("ND90Pr|GCF_000338995.1", 0.0) * 100.0

        # Evaluate 21-mer guides inside this candidate
        n_guides = max(0, len(seq) - 21 + 1)
        eff_guides = 0
        acc_vals = []
        for g_i in range(n_guides):
            g_seq = seq[g_i : g_i + 21]
            info = sifi_lookup.get(f"{tname}:{g_seq}")
            if info:
                if info["is_efficient_sifi21"]:
                    eff_guides += 1
                acc_vals.append(float(info["accessibility_prob"]))
            if g_seq not in guide_seen:
                guide_seen.add(g_seq)
                gid = f"GUIDE_{tname}_{len(guide_seen):05d}"
                all_guide_recs.append(SeqRecord(Seq(g_seq), id=gid, description=f"target={tname} candidate={cid}"))
                all_guides_rows.append(
                    {
                        "guide_id": gid,
                        "target": tname,
                        "first_candidate_id": cid,
                        "guide_sequence_sense_21nt": g_seq,
                        "guide_antisense_rna": str(Seq(g_seq).reverse_complement()).replace("T", "U"),
                        "gc_percent": round((g_seq.count("G") + g_seq.count("C")) / 21.0 * 100.0, 2),
                        "is_efficient_sifi21": info["is_efficient_sifi21"] if info else False,
                        "accessibility_prob": info["accessibility_prob"] if info else "NA",
                        "delta_end_stability": info["delta_end_stability"] if info else "NA",
                    }
                )

        mean_acc = float(rc["mean_accessibility"]) if pd.notna(rc["mean_accessibility"]) else (sum(acc_vals) / len(acc_vals) if acc_vals else 0.0)
        eff_frac = float(rc["efficient_sirna_frac"]) if pd.notna(rc["efficient_sirna_frac"]) else (eff_guides / max(n_guides, 1))
        qc_pass = (
            150 <= len(seq) <= 550
            and seq_stats["n_bases"] == 0
            and 30.0 <= seq_stats["gc_percent"] <= 65.0
            and max_hp <= 5
            and paralog_hits == 0
        )

        # Save per-isolate conservation row for bootstrap in M09
        iso_row = {"candidate_id": cid}
        iso_row.update(icons["per_isolate_fractions"])
        iso_cons_matrix_rows.append(iso_row)

        final_cand_rows.append(
            {
                "candidate_id": cid,
                "target": tname,
                "reference_assembly": "GCF_000338995.1",
                "isolate": "ND90Pr",
                "gene_id": target_meta_map[tname]["nd90pr_gene_id"],
                "locus_tag": locus_tag,
                "cds_length_nt": len(cds_ref),
                "candidate_start_nt": rc["start_1based"],
                "candidate_end_nt": rc["end_1based"],
                "candidate_length_nt": len(seq),
                "gc_percent": seq_stats["gc_percent"],
                "n_bases": seq_stats["n_bases"],
                "max_homopolymer": max_hp,
                "exact_cds_substring": exact_in_ref,
                "ref_identity_or_21mer_pct": round(ref_id_pct, 2),
                "isolates_evaluated": icons["n_isolates_evaluated"],
                "isolates_ge90pct_21mer": icons["isolates_ge90pct_21mer"],
                "isolates_100pct_21mer": icons["isolates_100pct_21mer"],
                "conservation_pct": icons["mean_21mer_conservation_pct"],
                "min_isolate_conservation_pct": icons["min_21mer_conservation_pct"],
                "tool_name": rc["tool"],
                "tool_version": rc["tool_version"],
                "tool_score": rc["tool_score"],
                "tool_score_metric": rc["tool_score_metric"],
                "total_21mer_guides": n_guides,
                "efficient_sifi21_guides": eff_guides,
                "efficient_guide_fraction": round(eff_frac, 4),
                "mean_target_accessibility": round(mean_acc, 4),
                "self_paralog_21mer_hits": paralog_hits,
                "primer_left": rc["primer_left"],
                "primer_right": rc["primer_right"],
                "expected_amplicon_bp": rc["amplicon_bp"],
                "seq_qc_pass": qc_pass,
                "off_target_status": "PENDING_M08",
                "confidence_classification": "High sequence conservation + QC pass" if (qc_pass and icons["mean_21mer_conservation_pct"] >= 95.0) else "Moderate/Conditional",
                "sequence": seq,
            }
        )
        final_cand_recs.append(SeqRecord(Seq(seq), id=cid, description=f"target={tname} tool={rc['tool']} len={len(seq)} cons={icons['mean_21mer_conservation_pct']}%"))

    cand_df = pd.DataFrame(final_cand_rows)
    guides_df = pd.DataFrame(all_guides_rows)
    comp_df = pd.DataFrame(tool_comp_rows)
    iso_mat_df = pd.DataFrame(iso_cons_matrix_rows)

    cand_fa_p = cand_dir / "dsRNA_candidates.fasta"
    cand_tsv_p = cand_dir / "dsRNA_candidates.tsv"
    guide_fa_p = guide_dir / "guide_sequences.fasta"
    guide_tsv_p = guide_dir / "guide_evaluation.tsv"
    comp_tsv_p = comp_dir / "candidate_tool_comparison.tsv"
    iso_mat_p = cand_dir / "candidate_isolate_conservation_matrix.tsv"

    SeqIO.write(final_cand_recs, str(cand_fa_p), "fasta")
    cand_df.to_csv(cand_tsv_p, sep="\t", index=False)
    SeqIO.write(all_guide_recs, str(guide_fa_p), "fasta")
    guides_df.to_csv(guide_tsv_p, sep="\t", index=False)
    comp_df.to_csv(comp_tsv_p, sep="\t", index=False)
    iso_mat_df.to_csv(iso_mat_p, sep="\t", index=False)

    # Copy to 10_final_reports
    fin_dir = cfg.module_dir("10_final_reports")
    SeqIO.write(final_cand_recs, str(fin_dir / "BRD_dsRNA_Candidates.fasta"), "fasta")
    cand_df.to_csv(fin_dir / "BRD_dsRNA_Candidates.tsv", sep="\t", index=False)
    SeqIO.write(all_guide_recs, str(fin_dir / "BRD_Guide_Sequences.fasta"), "fasta")
    guides_df.to_csv(fin_dir / "BRD_Guide_Evaluation.tsv", sep="\t", index=False)
    comp_df.to_csv(fin_dir / "BRD_Tool_Comparison.tsv", sep="\t", index=False)

    # Figures (PNG + SVG): length distribution, GC distribution, conservation comparison
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    axes[0].hist(cand_df["candidate_length_nt"], bins=10, color="#2563eb", edgecolor="#0f172a", alpha=0.85)
    axes[0].set_xlabel("Candidate Length (nt)")
    axes[0].set_ylabel("Candidate Count")
    axes[0].set_title("Candidate dsRNA Length Distribution")

    axes[1].hist(cand_df["gc_percent"], bins=10, color="#059669", edgecolor="#0f172a", alpha=0.85)
    axes[1].axvline(30, color="#dc2626", linestyle="--", linewidth=1)
    axes[1].axvline(65, color="#dc2626", linestyle="--", linewidth=1)
    axes[1].set_xlabel("GC Content (%)")
    axes[1].set_ylabel("Candidate Count")
    axes[1].set_title("Candidate GC% Distribution")

    for tname, col in zip(["CsPKS1", "CsFUS3", "CsSTE11", "CsSp1"], ["#7c3aed", "#2563eb", "#059669", "#dc2626"]):
        sub = cand_df[cand_df["target"] == tname]
        if not sub.empty:
            axes[2].scatter(sub["conservation_pct"], sub["efficient_guide_fraction"] * 100.0, label=tname, color=col, s=55, alpha=0.85)
    axes[2].set_xlabel("Mean Multi-Isolate 21-mer Conservation (%)")
    axes[2].set_ylabel("sifi21_plus Efficient Guide Fraction (%)")
    axes[2].set_title("Conservation vs Guide Efficiency")
    axes[2].legend(fontsize=8)
    axes[2].grid(True, linestyle="--", alpha=0.35)

    fig.suptitle(f"BRD M07 dsRNA Candidate Characteristics (n = {len(cand_df)} candidates across 4 targets)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(fig_dir / "candidate_characteristics.png", dpi=300)
    fig.savefig(fig_dir / "candidate_characteristics.svg")
    plt.close(fig)

    # Individual named figures requested in Section 49
    for col_name, fname, title, xlabel in [
        ("candidate_length_nt", "candidate_length_distribution", "Candidate dsRNA Length Distribution", "Length (nt)"),
        ("gc_percent", "candidate_gc_distribution", "Candidate dsRNA GC Content Distribution", "GC Content (%)"),
        ("conservation_pct", "candidate_conservation_comparison", "Candidate 21-mer Conservation Across 19 Isolates", "Mean 21-mer Conservation (%)"),
    ]:
        f_i, a_i = plt.subplots(figsize=(7, 4.5))
        for tname in ["CsPKS1", "CsFUS3", "CsSTE11", "CsSp1"]:
            vals = cand_df[cand_df["target"] == tname][col_name]
            if not vals.empty:
                a_i.hist(vals, bins=8, alpha=0.55, label=tname)
        a_i.set_xlabel(xlabel)
        a_i.set_ylabel("Number of Candidates")
        a_i.set_title(title)
        a_i.legend()
        f_i.tight_layout()
        f_i.savefig(fig_dir / f"{fname}.png", dpi=300)
        f_i.savefig(fig_dir / f"{fname}.svg")
        plt.close(f_i)

    body = (
        f"<div class='warn'><strong>Scope Boundary:</strong> Computational dsRNA candidate generation evaluates sequence conservation, "
        f"thermodynamic asymmetry, accessibility, and primer feasibility. It does <strong>not</strong> establish exogenous dsRNA uptake, "
        f"Dicer processing, gene silencing, or disease suppression in <em>B. sorokiniana</em>.</div>"
        f"<h2>Tool Comparison Summary</h2>{df_to_html_table(comp_df)}"
        f"<h2>Top Designed Candidates</h2>{df_to_html_table(cand_df.drop(columns=['sequence']).head(40))}"
        f"<figure><img src='figures/candidate_characteristics.png'/><figcaption>Length, GC%, and conservation vs efficient siRNA fraction.</figcaption></figure>"
    )
    rep_html = m07_dir / "candidate_design_report.html"
    write_text(rep_html, wrap_html_report("BRD M07 — dsRNA Candidate Design Report", "LUCID/sifi21_plus, SIREN, dsRNAmax, and BRD Isolate-Aware Design", body))

    prov.record("M07", "dsRNA Candidate Design", str(ref_rna_fa), cand_tsv_p, "LUCID+SIREN+dsRNAmax+ViennaRNA", "1.0.0", f"{len(cand_df)}_candidates")
    logger.info("M07 complete: %d unique dsRNA candidates and %d unique 21-mer guides generated.", len(cand_df), len(guides_df))
    return {
        "candidates_fasta": cand_fa_p,
        "candidates_tsv": cand_tsv_p,
        "guides_fasta": guide_fa_p,
        "guides_tsv": guide_tsv_p,
        "tool_comparison": comp_tsv_p,
        "isolate_conservation_matrix": iso_mat_p,
        "report": rep_html,
    }
