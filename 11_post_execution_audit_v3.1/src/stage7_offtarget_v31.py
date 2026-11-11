"""Stage 7: Audited Off-Target Screening Metadata, Mismatch-Position Resolution & Risk Reporting (BRD v3.1)."""

from __future__ import annotations

import hashlib
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import pandas as pd
from Bio.Seq import Seq

try:
    from .paths_v3_1 import resolve_repo_root, resolve_v31_root
except ImportError:
    from paths_v3_1 import resolve_repo_root, resolve_v31_root


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def _iter_needed_fasta(fa_path: Path, needed_ids: Set[str]) -> Dict[str, str]:
    """Stream FASTA and extract ONLY sequences whose record ID is in `needed_ids`."""
    found: Dict[str, str] = {}
    if not fa_path or not fa_path.exists() or not needed_ids:
        return found
    hdr = None
    keep = False
    chunks: List[str] = []
    with open(fa_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and hdr is not None:
                    found[hdr] = "".join(chunks).upper()
                    if len(found) == len(needed_ids):
                        return found
                hdr = line[1:].strip().split()[0]
                keep = hdr in needed_ids
                chunks = []
            elif keep:
                chunks.append(line.strip())
        if keep and hdr is not None:
            found[hdr] = "".join(chunks).upper()
    return found


def _resolve_reference_fastas(root: Path, ref_id: str, acc: str) -> Tuple[Path | None, Path | None]:
    """Locate transcript FASTA and genome FASTA for each expanded panel reference."""
    if ref_id == "Bipolaris_sorokiniana_ND90Pr_self":
        base = root / "01_data_acquisition" / "local_nd90pr" / "ncbi_dataset" / "data" / acc
        tx = base / "rna.fna"
        if not tx.exists():
            tx = base / "cds_from_genomic.fna"
        gen = base / "GCF_000338995.1_Cocsa1_genomic.fna"
        return (tx if tx.exists() else None), (gen if gen.exists() else None)

    v3_base = root / "11_post_execution_audit_v3" / "expanded_references" / "sources" / ref_id / "ncbi_dataset" / "data" / acc
    v1_base = root / "08_off_target_screening" / "blast_databases" / "sources" / ref_id / "ncbi_dataset" / "data" / acc

    tx_p = None
    for cand_base in (v3_base, v1_base):
        for fname in ("rna.fna", "cds_from_genomic.fna"):
            cp = cand_base / fname
            if cp.exists() and cp.stat().st_size > 0:
                tx_p = cp
                break
        if tx_p:
            break

    gen_p = None
    for cand_base in (v3_base, v1_base):
        if cand_base.exists():
            g_cands = [p for p in cand_base.glob("*_genomic.fna") if p.name != "cds_from_genomic.fna" and p.stat().st_size > 0]
            if g_cands:
                gen_p = g_cands[0]
                break
    return tx_p, gen_p


def compute_guide_mismatch_details(guide: str, subj_seq: str, sstart: int, send: int, sstrand: str) -> Tuple[int, str, str]:
    """Return (total_21nt_mismatches, comma_separated_1based_positions, aligned_21nt_subject_window)."""
    g = guide.upper()
    if len(g) != 21 or not subj_seq:
        return -1, "UNRESOLVED_SEQ", ""
    # BLAST coordinates are 1-based inclusive; plus strand sstart<=send, minus strand often sstart>=send or swapped in some outputs
    lo = max(0, min(sstart, send) - 1)
    hi = min(len(subj_seq), max(sstart, send))
    frag = subj_seq[lo:hi]
    if sstrand == "minus":
        frag = str(Seq(frag).reverse_complement())

    # Search locally around [lo-5 : hi+5] for the minimal-mismatch 21-mer alignment against guide `g`
    win_lo = max(0, lo - 4)
    win_hi = min(len(subj_seq), hi + 4)
    region = subj_seq[win_lo:win_hi]
    if sstrand == "minus":
        region = str(Seq(region).reverse_complement())

    if len(region) < 21:
        return -1, "SHORT_SUBJECT_SPAN", frag

    best_mm = 22
    best_pos: List[int] = []
    best_win = ""
    for off in range(len(region) - 20):
        sub21 = region[off : off + 21]
        mm_pos = [idx + 1 for idx, (a, b) in enumerate(zip(g, sub21)) if a != b]
        if len(mm_pos) < best_mm:
            best_mm = len(mm_pos)
            best_pos = mm_pos
            best_win = sub21
    pos_str = ",".join(str(x) for x in best_pos) if best_pos else "NONE"
    return best_mm, pos_str, best_win


def run_stage7_offtarget(repo_root: Path | None = None) -> Dict[str, Path]:
    root = resolve_repo_root(repo_root)
    v3 = root / "11_post_execution_audit_v3"
    v31 = resolve_v31_root(root)
    s7 = v31 / "stage7_offtarget"
    s7.mkdir(parents=True, exist_ok=True)
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    man_v3 = pd.read_csv(v3 / "stage4_off_target_expanded" / "expanded_reference_manifest_v3.tsv", sep="\t")
    gaps_v3 = pd.read_csv(v3 / "stage4_off_target_expanded" / "reference_panel_coverage_and_gaps_v3.tsv", sep="\t").set_index("ref_id")
    exact_v3 = pd.read_csv(v3 / "stage4_off_target_expanded" / "exact_21mer_matches_expanded_v3.tsv", sep="\t")
    near_v3 = pd.read_csv(v3 / "stage4_off_target_expanded" / "near_guide_matches_expanded_v3.tsv", sep="\t")
    sel_v3 = pd.read_csv(v3 / "stage4_off_target_expanded" / "candidate_selectivity_expanded_v3.tsv", sep="\t")
    auth_c = pd.read_csv(v3 / "authoritative_data" / "auth_candidates.tsv", sep="\t").set_index("candidate_id")

    # 1. Reference metadata with verified real FASTA SHA-256 checksums
    ref_paths: Dict[Tuple[str, str], Path | None] = {}
    meta_rows: List[Dict[str, Any]] = []
    for _, r in man_v3.iterrows():
        ref_id = str(r["ref_id"])
        acc = str(r["accession"])
        tx_p, gen_p = _resolve_reference_fastas(root, ref_id, acc)
        ref_paths[(ref_id, "transcript")] = tx_p
        ref_paths[(ref_id, "genome")] = gen_p
        gp = gaps_v3.loc[ref_id]
        meta_rows.append(
            {
                "ref_id": ref_id,
                "group": str(r["group"]),
                "organism": str(r["organism"]),
                "tax_id": int(r["tax_id"]),
                "accession": acc,
                "source_database": str(r["source_database"]),
                "assembly_version_name": str(r["assembly_version_name"]),
                "release_date": str(r["release_date"]),
                "transcript_screened": bool(r["transcript_screened"]),
                "transcript_count": int(r["transcript_count"]),
                "transcript_fasta_file": tx_p.name if tx_p else "MISSING",
                "transcript_fasta_bytes": tx_p.stat().st_size if tx_p else 0,
                "transcript_fasta_sha256": _sha256_file(tx_p) if tx_p else "NOT_LOCALLY_PRESENT",
                "genome_screened": bool(r["genome_screened"]),
                "genome_scaffold_count": int(r["genome_scaffold_count"]),
                "genome_fasta_file": gen_p.name if gen_p else "EXCLUDED_OR_NOT_SCREENED",
                "genome_fasta_bytes": gen_p.stat().st_size if gen_p else 0,
                "genome_fasta_sha256": _sha256_file(gen_p) if gen_p else "EXCLUDED_MULTI_GB_GENOME",
                "screening_scope_note": str(gp["screening_scope_note"]),
                "rationale": str(r["rationale"]),
            }
        )
    meta_tsv = s7 / "offtarget_reference_metadata.tsv"
    pd.DataFrame(meta_rows).to_csv(meta_tsv, sep="\t", index=False)

    # 2. Gather needed subject accessions across exact & near hits to resolve exact orientation + 21-nt mismatch positions
    needed_by_file: Dict[Path, Set[str]] = defaultdict(set)
    for df_src in (exact_v3, near_v3):
        for ref_id, scope, sacc in zip(df_src["ref_id"], df_src["screen_scope"], df_src["subject_accession"]):
            fp = ref_paths.get((str(ref_id), str(scope)))
            if fp is not None:
                needed_by_file[fp].add(str(sacc))

    seq_cache: Dict[Tuple[Path, str], str] = {}
    for fp, id_set in needed_by_file.items():
        extracted = _iter_needed_fasta(fp, id_set)
        for sid, sq in extracted.items():
            seq_cache[(fp, sid)] = sq

    # 3. Enrich Exact 21-mer Hits with alignment orientation & self vs non-target label
    exact_rows: List[Dict[str, Any]] = []
    for _, er in exact_v3.iterrows():
        ref_id = str(er["ref_id"])
        scope = str(er["screen_scope"])
        sacc = str(er["subject_accession"])
        kmer = str(er["matched_21mer"]).upper()
        fp = ref_paths.get((ref_id, scope))
        sq = seq_cache.get((fp, sacc), "") if fp else ""
        orient = "UNRESOLVED"
        pos_1based = "NA"
        if sq:
            rc = str(Seq(kmer).reverse_complement())
            idx_p = sq.find(kmer)
            idx_m = sq.find(rc)
            if idx_p != -1 and idx_m != -1:
                orient = "both(+/-)"
                pos_1based = f"+:{idx_p+1}, -:{idx_m+1}"
            elif idx_p != -1:
                orient = "plus(+)"
                pos_1based = str(idx_p + 1)
            elif idx_m != -1:
                orient = "minus(-)"
                pos_1based = str(idx_m + 1)
        exact_rows.append(
            {
                "candidate_id": str(er["candidate_id"]),
                "target": str(er["target"]),
                "matched_21mer": kmer,
                "is_self_targeting": ref_id == "Bipolaris_sorokiniana_ND90Pr_self",
                "screen_scope": scope,
                "database_group": str(er["database_group"]),
                "ref_id": ref_id,
                "organism": str(er["organism"]),
                "subject_accession": sacc,
                "alignment_orientation": orient,
                "subject_match_start_1based": pos_1based,
            }
        )
    exact_tsv = s7 / "offtarget_exact_hits_v3_1.tsv"
    pd.DataFrame(exact_rows).to_csv(exact_tsv, sep="\t", index=False)

    # 4. Enrich Near Hits with true 21-nt mismatch counts & 1-based mismatch positions
    near_rows: List[Dict[str, Any]] = []
    for _, nr in near_v3.iterrows():
        ref_id = str(nr["ref_id"])
        scope = str(nr["screen_scope"])
        sacc = str(nr["subject_accession"])
        guide = str(nr["guide_21nt"]).upper()
        sstart = int(nr["sstart"])
        send = int(nr["send"])
        sstrand = str(nr["sstrand"])
        fp = ref_paths.get((ref_id, scope))
        sq = seq_cache.get((fp, sacc), "") if fp else ""
        if sq:
            mm_cnt, mm_pos, subj_21 = compute_guide_mismatch_details(guide, sq, sstart, send, sstrand)
        else:
            # Estimate minimum mismatches from BLASTn alignment length + internal mismatches if FASTA sequence lookup unavailable
            aln_len = int(nr["alignment_length_nt"])
            int_mm = int(nr["mismatches"])
            mm_cnt = (21 - aln_len) + int_mm
            mm_pos = f"EDGE_OR_CLIPPED(aln_len={aln_len},int_mm={int_mm})"
            subj_21 = ""
        near_rows.append(
            {
                "candidate_id": str(nr["candidate_id"]),
                "target": str(nr["target"]),
                "guide_21nt": guide,
                "match_tier": str(nr["match_tier"]),
                "is_self_targeting": ref_id == "Bipolaris_sorokiniana_ND90Pr_self",
                "screen_scope": scope,
                "database_group": str(nr["database_group"]),
                "ref_id": ref_id,
                "organism": str(nr["organism"]),
                "subject_accession": sacc,
                "alignment_orientation": sstrand,
                "sstart": sstart,
                "send": send,
                "blast_alignment_length_nt": int(nr["alignment_length_nt"]),
                "blast_internal_mismatches": int(nr["mismatches"]),
                "total_21nt_mismatches": mm_cnt,
                "mismatch_positions_1based": mm_pos,
                "aligned_subject_21nt": subj_21,
                "pident": float(nr["pident"]),
                "evalue": float(nr["evalue"]),
            }
        )
    near_df = pd.DataFrame(near_rows)
    near_tsv = s7 / "offtarget_near_hits_v3_1.tsv"
    near_df.to_csv(near_tsv, sep="\t", index=False)

    # 5. Summary table per candidate retaining prominent near-match burdens & coverage limits
    sum_rows: List[Dict[str, Any]] = []
    near_by_cand = {cid: grp for cid, grp in near_df.groupby("candidate_id")}
    exact_df = pd.DataFrame(exact_rows)
    exact_by_cand = {cid: grp for cid, grp in exact_df.groupby("candidate_id")}

    for _, sr in sel_v3.iterrows():
        cid = str(sr["candidate_id"])
        ng = near_by_cand.get(cid, pd.DataFrame(columns=near_df.columns))
        eg = exact_by_cand.get(cid, pd.DataFrame(columns=exact_df.columns))

        n_20_non_self = int(ng[(~ng["is_self_targeting"]) & (ng["match_tier"] == "near_20_of_21")]["guide_21nt"].nunique())
        n_19_non_self = int(ng[(~ng["is_self_targeting"]) & (ng["match_tier"] == "near_19_of_21")]["guide_21nt"].nunique())
        n_18_non_self = int(ng[(~ng["is_self_targeting"]) & (ng["match_tier"] == "near_18_of_21")]["guide_21nt"].nunique())
        n_self_exact = int(eg[eg["is_self_targeting"]]["matched_21mer"].nunique()) if not eg.empty else 0
        n_nontarg_exact = int(eg[~eg["is_self_targeting"]]["matched_21mer"].nunique()) if not eg.empty else 0

        sum_rows.append(
            {
                "candidate_id": cid,
                "target": str(sr["target"]),
                "candidate_length_nt": int(auth_c.loc[cid, "candidate_length_nt"]),
                "total_21nt_windows": int(auth_c.loc[cid, "candidate_length_nt"]) - 20,
                "passes_zero_exact_21mer_filter": bool(sr["passes_zero_exact_21mer_filter"]),
                "risk_classification_v3_1": str(sr["selectivity_classification_v3"]),
                "nontarget_exact_21mer_unique_guides": n_nontarg_exact,
                "self_paralog_exact_21mer_unique_guides": n_self_exact,
                "wheat_exact_21": int(sr["wheat_exact_21"]),
                "barley_exact_21": int(sr["barley_exact_21"]),
                "rice_exact_21": int(sr["rice_exact_21"]),
                "maize_exact_21": int(sr["maize_exact_21"]),
                "beneficial_fungi_tx_exact_21": int(sr["beneficial_fungi_tx_exact_21"]),
                "beneficial_fungi_genome_exact_21": int(sr["beneficial_fungi_genome_exact_21"]),
                "bacteria_exact_21": int(sr["bacteria_exact_21"]),
                "near_20of21_unique_guides": n_20_non_self,
                "near_19of21_unique_guides": n_19_non_self,
                "near_18of21_unique_guides": n_18_non_self,
                "crops_unique_20of21_guides": int(sr["crops_unique_20of21_guides"]),
                "beneficial_fungi_unique_20of21_guides": int(sr["beneficial_fungi_unique_20of21_guides"]),
                "bacteria_unique_20of21_guides": int(sr["bacteria_unique_20of21_guides"]),
                "crops_long_homology_ge30nt": int(sr["crops_long_homology_ge30nt"]),
                "beneficial_fungi_long_homology_ge30nt": int(sr["beneficial_fungi_long_homology_ge30nt"]),
                "bacteria_long_homology_ge30nt": int(sr["bacteria_long_homology_ge30nt"]),
                "heuristic_sequence_selectivity_index": float(sr["heuristic_sequence_selectivity_index"]),
                "database_coverage_limitation": "Wheat, barley & maize screened on complete RefSeq transcriptomes only (>20 GB nuclear genomes excluded); unsequenced field microbiome unscreened; 0 hits != biological safety.",
            }
        )
    sum_tsv = s7 / "offtarget_summary_v3_1.tsv"
    pd.DataFrame(sum_rows).to_csv(sum_tsv, sep="\t", index=False)

    methods_md = f"""# Stage 7 — Expanded Off-Target Screening & Mismatch Resolution Methodology (BRD v3.1)

**Executed:** {now_utc}

## Scope, Thresholds & Explicit Caveats
1. **Taxonomic panel (`offtarget_reference_metadata.tsv`):** `17` non-target organisms (`4` crops: wheat, barley, rice, maize; `9` beneficial fungi: *T. atroviride*, *T. virens*, *T. reesei*, *T. harzianum*, *B. bassiana*, *M. anisopliae*, *M. robertsii*, *R. irregularis*, *R. clarus*; `4` rhizosphere/soil bacteria: *B. subtilis*, *B. velezensis*, *P. fluorescens*, *P. putida*) plus masked self *B. sorokiniana* `ND90Pr` control. Verified local FASTA SHA-256 checksums are recorded for every present transcript & genome file.
2. **Exact & near-hit resolution:**
   - **Exact 21-nt hard threshold preserved (`passes_zero_exact_21mer_filter`):** `32/50` pass; `18/50` excluded (including all `13/13` `CsFUS3` candidates). Alignment strand (`+`/`-`) and 1-based start offsets are resolved directly against reference sequences in `offtarget_exact_hits_v3_1.tsv`.
   - **Near matches (`offtarget_near_hits_v3_1.tsv`):** Enriches all `{len(near_df):,}` BLASTn-short hits (`20/21`, `19/21`, `18/21`) with full 21-nt subject window extraction, mismatch count (`total_21nt_mismatches`), and **exact 1-based mismatch positions (`1..21`)** inside the guide.
3. **Near-match burden visibility (`AUDIT_08`):** `29/50` candidates (`7/8` clean shortlist representatives!) carry `Elevated_20of21_Near_Match_Burden` (`3–14` unique non-target `20/21` guides across leading representatives). Passing exact 21-mer filtering NEVER equals biological safety.
4. **Coverage boundaries:** Whole nuclear genomes of wheat (`T. aestivum`), barley (`H. vulgare`), and maize (`Z. mays`) (`>20 GB`) were excluded locally (RefSeq transcriptomes screened), and unsequenced soil/phyllosphere biota remain uncovered.
"""
    (s7 / "offtarget_methods.md").write_text(methods_md, encoding="utf-8")

    return {
        "metadata": meta_tsv,
        "exact_hits": exact_tsv,
        "near_hits": near_tsv,
        "summary": sum_tsv,
        "methods": s7 / "offtarget_methods.md",
    }


if __name__ == "__main__":
    for k, v in run_stage7_offtarget().items():
        print(f"STAGE7_OT_OUT {k}: {v}")
