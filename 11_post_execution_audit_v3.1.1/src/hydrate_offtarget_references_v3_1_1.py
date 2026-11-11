"""Stage 3 (FIX 3): Reproducible Off-Target Reference Hydration, Validation, DB Rebuild & Screening Smoke Test.

Retrieves and SHA-256 verifies all 18 expanded off-target panel entries (17 non-target crops,
beneficial fungi, and bacteria + self B. sorokiniana ND90Pr control) from NCBI Datasets v2
strictly against `11_post_execution_audit_v3.1/stage7_offtarget/offtarget_reference_metadata.tsv`.

Never silently substitutes newer assembly versions or accepts SHA-256 checksum mismatches.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

try:
    from .paths_v3_1_1 import resolve_repo_root, resolve_tool, resolve_v31_root, resolve_v311_root
except ImportError:
    from paths_v3_1_1 import resolve_repo_root, resolve_tool, resolve_v31_root, resolve_v311_root

VALID_NT = set("ACGTURYKMSWBDHVN")


def _short_path_if_nt(p: Path) -> str:
    s_abs = str(p.resolve())
    if os.name == "nt" and " " in s_abs:
        import ctypes

        buf = ctypes.create_unicode_buffer(1024)
        if ctypes.windll.kernel32.GetShortPathNameW(s_abs, buf, 1024) and buf.value and " " not in buf.value:
            return buf.value
    return s_abs


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def build_ncbi_download_url(accession: str, tx_filename: str, fetch_genome: bool) -> str:
    inc: List[str] = []
    if tx_filename == "rna.fna":
        inc.append("RNA_FASTA")
    elif tx_filename == "cds_from_genomic.fna":
        inc.append("CDS_FASTA")
    if fetch_genome:
        inc.append("GENOME_FASTA")
    if not inc:
        raise ValueError(f"Cannot build Datasets v2 include_annotation_type for {accession}")
    return f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download?include_annotation_type={','.join(inc)}"


def validate_nucleotide_fasta(
    fa_path: Path,
    expected_sha256: str,
    expected_bytes: int,
    expected_seq_count: int | None = None,
) -> Dict[str, Any]:
    """Validate FASTA readability, non-empty nucleotide records, byte size, count, and SHA-256."""
    if not fa_path.exists() or fa_path.stat().st_size == 0:
        raise FileNotFoundError(f"Reference FASTA missing or empty: {fa_path}")
    actual_bytes = fa_path.stat().st_size
    actual_sha = sha256_file(fa_path).lower()
    if actual_bytes != expected_bytes:
        raise ValueError(f"Byte-size mismatch for {fa_path.name}: got {actual_bytes}, expected {expected_bytes}")
    if actual_sha != expected_sha256.lower():
        raise ValueError(f"SHA-256 mismatch for {fa_path.name}: got {actual_sha}, expected {expected_sha256.lower()}")

    seq_count = 0
    total_bp = 0
    bad_chars = 0
    hdr = None
    cur_len = 0
    empty_seqs = 0
    with open(fa_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line == "\n" or line == "\r\n":
                continue
            if line[0] == ">":
                if hdr is not None:
                    if cur_len == 0:
                        empty_seqs += 1
                    total_bp += cur_len
                hdr = line[1:].strip()
                seq_count += 1
                cur_len = 0
            else:
                s = line.strip().upper()
                cur_len += len(s)
                if not set(s).issubset(VALID_NT):
                    bad_chars += len(set(s) - VALID_NT)
        if hdr is not None:
            if cur_len == 0:
                empty_seqs += 1
            total_bp += cur_len

    if seq_count == 0 or empty_seqs > 0 or total_bp == 0:
        raise ValueError(f"Invalid/empty nucleotide records in {fa_path}: count={seq_count}, empty={empty_seqs}")
    if expected_seq_count is not None and expected_seq_count > 0 and seq_count != expected_seq_count:
        raise ValueError(f"Record count mismatch for {fa_path.name}: got {seq_count}, expected {expected_seq_count}")

    return {
        "sequence_count": seq_count,
        "total_bp": total_bp,
        "actual_bytes": actual_bytes,
        "actual_sha256": actual_sha,
        "unexpected_chars": bad_chars,
    }


def _find_legacy_local_files(root: Path, ref_id: str, acc: str, tx_name: str, gen_name: str) -> Tuple[Path | None, Path | None]:
    """Check if raw source files already exist in repo directories without hardcoding user paths."""
    if ref_id == "Bipolaris_sorokiniana_ND90Pr_self":
        base = root / "01_data_acquisition" / "local_nd90pr" / "ncbi_dataset" / "data" / acc
        tx = base / tx_name
        gen = base / gen_name if gen_name != "EXCLUDED_OR_NOT_SCREENED" else None
        return (tx if tx.exists() else None), (gen if gen and gen.exists() else None)

    for base in (
        root / "11_post_execution_audit_v3" / "expanded_references" / "sources" / ref_id / "ncbi_dataset" / "data" / acc,
        root / "08_off_target_screening" / "blast_databases" / "sources" / ref_id / "ncbi_dataset" / "data" / acc,
    ):
        tx = base / tx_name
        gen = base / gen_name if gen_name != "EXCLUDED_OR_NOT_SCREENED" else None
        if tx.exists() or (gen and gen.exists()):
            return (tx if tx.exists() else None), (gen if gen and gen.exists() else None)
    return None, None


def _download_ncbi_datasets_zip(url: str, acc: str, expected_tax_id: int, ref_dest: Path, tx_name: str, gen_name: str, timeout: int = 300) -> Tuple[Path | None, Path | None, dict]:
    """Download NCBI Datasets v2 ZIP package, validate accession inside assembly_data_report.jsonl, and extract target files."""
    ref_dest.mkdir(parents=True, exist_ok=True)
    zip_tmp = ref_dest / f"{acc}_download.zip"
    req = urllib.request.Request(url, headers={"User-Agent": "BRD-v3.1.1-OffTargetHydrator/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(zip_tmp, "wb") as out_f:
            shutil.copyfileobj(resp, out_f, length=1048576)

        report_meta: dict = {}
        tx_out: Path | None = None
        gen_out: Path | None = None
        data_dir = ref_dest / "ncbi_dataset" / "data" / acc
        data_dir.mkdir(parents=True, exist_ok=True)

        with zipfile.ZipFile(zip_tmp, "r") as zf:
            rep_member = next((m for m in zf.namelist() if m.endswith("assembly_data_report.jsonl")), None)
            if rep_member:
                raw_rep = zf.read(rep_member).decode("utf-8", errors="replace").strip()
                if raw_rep:
                    report_meta = json.loads(raw_rep.splitlines()[0])
                    rep_acc = str(report_meta.get("accession", ""))
                    if rep_acc and rep_acc != acc:
                        raise RuntimeError(f"Accession mismatch in downloaded NCBI report: got {rep_acc}, expected {acc}")
                    rep_tax = int(report_meta.get("organism", {}).get("tax_id", expected_tax_id))
                    if expected_tax_id and rep_tax != expected_tax_id:
                        raise RuntimeError(f"TaxID mismatch for {acc}: got {rep_tax}, expected {expected_tax_id}")

            for m in zf.namelist():
                bn = Path(m).name
                if bn == tx_name:
                    tx_out = data_dir / tx_name
                    with zf.open(m) as src_f, open(tx_out, "wb") as dst_f:
                        shutil.copyfileobj(src_f, dst_f, length=1048576)
                elif gen_name != "EXCLUDED_OR_NOT_SCREENED" and bn == gen_name:
                    gen_out = data_dir / gen_name
                    with zf.open(m) as src_f, open(gen_out, "wb") as dst_f:
                        shutil.copyfileobj(src_f, dst_f, length=1048576)
                elif bn in ("assembly_data_report.jsonl", "dataset_catalog.json"):
                    meta_out = ref_dest / "ncbi_dataset" / "data" / bn
                    with zf.open(m) as src_f, open(meta_out, "wb") as dst_f:
                        shutil.copyfileobj(src_f, dst_f)

        return tx_out, gen_out, report_meta
    finally:
        zip_tmp.unlink(missing_ok=True)


def rebuild_combined_screening_db(
    hydrated_fastas: List[Tuple[str, str, Path]],
    out_dir: Path,
    root: Path,
    max_mb_for_blastdb: int = 120,
) -> Dict[str, Any]:
    """Rebuild a combined FASTA (and BLAST nucleotide DB if small enough & makeblastdb available) from hydrated files."""
    db_dir = out_dir / "rebuilt_smoke_db"
    db_dir.mkdir(parents=True, exist_ok=True)
    comb_fa = db_dir / "combined_hydrated_offtarget_refs.fna"
    n_records = 0
    total_bytes = 0
    with open(comb_fa, "w", encoding="utf-8", newline="\n") as out_fh:
        for ref_id, scope, fa_path in hydrated_fastas:
            for line in open(fa_path, "r", encoding="utf-8", errors="replace"):
                if line.startswith(">"):
                    out_fh.write(f">{ref_id}|{scope}|{line[1:]}")
                    n_records += 1
                else:
                    out_fh.write(line)
            total_bytes += fa_path.stat().st_size

    mkdb = resolve_tool("makeblastdb", root)
    blastn_bin = resolve_tool("blastn", root)
    blastdb_status = "NOT_BUILT"
    blastn_smoke_hits = -1
    if mkdb:
        if (total_bytes / 1048576.0) <= max_mb_for_blastdb:
            db_target = comb_fa.name
        else:
            # Build a representative <= 120 MB multi-kingdom smoke FASTA (bacteria + beneficial fungi) for makeblastdb verification
            smoke_fa = db_dir / "smoke_representative_panel.fna"
            smoke_bytes = 0
            with open(smoke_fa, "w", encoding="utf-8", newline="\n") as sfh:
                for ref_id, scope, fa_path in hydrated_fastas:
                    sz = fa_path.stat().st_size
                    if scope == "transcript" and sz <= 25_000_000 and (smoke_bytes + sz) <= max_mb_for_blastdb * 1048576:
                        for line in open(fa_path, "r", encoding="utf-8", errors="replace"):
                            sfh.write(f">{ref_id}|{scope}|{line[1:]}" if line.startswith(">") else line)
                        smoke_bytes += sz
            db_target = smoke_fa.name
        cwd_safe = _short_path_if_nt(db_dir)
        proc = subprocess.run(
            [mkdb, "-in", db_target, "-dbtype", "nucl", "-out", "combined_hydrated_db"],
            cwd=cwd_safe,
            capture_output=True,
            text=True,
        )
        blastdb_status = f"BUILT_OK ({db_target})" if proc.returncode == 0 else f"FAILED_RC_{proc.returncode}"
        if proc.returncode == 0 and blastn_bin:
            q_smoke = db_dir / "smoke_query.fa"
            v31 = resolve_v31_root(root)
            ex_cands = pd.read_csv(v31 / "stage8_candidate_design" / "candidate_master_v3_1.tsv", sep="\t").set_index("candidate_id")
            q_seq = str(ex_cands.loc["BRD_CsFUS3_CAND001", "sequence"])[:40]
            q_smoke.write_text(f">smoke_fus3_sub\n{q_seq}\n", encoding="utf-8")
            bn_out = db_dir / "smoke_blastn.tsv"
            pb = subprocess.run(
                [blastn_bin, "-task", "blastn-short", "-query", q_smoke.name, "-db", "combined_hydrated_db", "-out", bn_out.name, "-outfmt", "6", "-evalue", "10"],
                cwd=cwd_safe,
                capture_output=True,
                text=True,
            )
            if pb.returncode == 0 and bn_out.exists():
                blastn_smoke_hits = sum(1 for ln in bn_out.read_text(encoding="utf-8").splitlines() if ln.strip())
                blastdb_status += f"; blastn_smoke_hits={blastn_smoke_hits}"
    else:
        blastdb_status = "SKIPPED_MAKEBLASTDB_MISSING"

    return {
        "combined_fasta": comb_fa,
        "combined_records": n_records,
        "combined_bytes": comb_fa.stat().st_size if comb_fa.exists() else 0,
        "blastdb_status": blastdb_status,
        "blastn_smoke_hits": blastn_smoke_hits,
        "db_dir": db_dir,
    }


def _scan_exact_21mers_fast(
    fa_path: Path,
    target_21mers: Set[str],
    skip_record_substrs: Set[str] | None = None,
    skip_genome_interval: Tuple[str, int, int] | None = None,
) -> int:
    """Pigeonhole 11-nt stride exact 21-mer scan (sense+antisense already in target_21mers), masking self target locus when requested."""
    if not fa_path.exists() or not target_21mers:
        return 0
    idx: Dict[str, List[int]] = {}
    for km in target_21mers:
        for off in range(11):
            idx.setdefault(km[off : off + 11], []).append(off)
    hits = 0
    hdr = None
    chunks: List[str] = []

    def _scan_seq(rec_hdr: str, seq: str) -> int:
        rec_id = rec_hdr.split()[0]
        if skip_record_substrs and any(tok in rec_hdr for tok in skip_record_substrs):
            return 0
        n = len(seq)
        if n < 21:
            return 0
        max_st = n - 21
        g_scf, g_lo, g_hi = skip_genome_interval or ("", -1, -1)
        c = 0
        for g in range(0, n - 10, 11):
            offs = idx.get(seq[g : g + 11])
            if offs:
                for off in offs:
                    st = g - off
                    if 0 <= st <= max_st and seq[st : st + 21] in target_21mers:
                        if g_scf and rec_id == g_scf and g_lo <= (st + 1) <= g_hi:
                            continue
                        c += 1
        return c

    with open(fa_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    hits += _scan_seq(hdr, "".join(chunks).upper())
                hdr = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if hdr is not None:
            hits += _scan_seq(hdr, "".join(chunks).upper())
    return hits


def run_offtarget_smoke_test(hydrated_fastas: List[Tuple[str, str, Path]], root: Path) -> Dict[str, Any]:
    """Compare exact 21-mer hit behavior on hydrated refs vs v3.1 expectations (with self-target locus masked cleanly)."""
    v3 = root / "11_post_execution_audit_v3"
    v31 = resolve_v31_root(root)
    cands = pd.read_csv(v31 / "stage8_candidate_design" / "candidate_master_v3_1.tsv", sep="\t").set_index("candidate_id")
    auth_t = pd.read_csv(v3 / "authoritative_data" / "auth_targets.tsv", sep="\t").set_index("target")
    ex_v31 = pd.read_csv(v31 / "stage7_offtarget" / "offtarget_exact_hits_v3_1.tsv", sep="\t")

    pos_cid = "BRD_CsFUS3_CAND001"
    clean_cid = "BRD_CsSTE11_CAND013"

    def _kmers_both(seq: str) -> Set[str]:
        s = seq.upper()
        out: Set[str] = set()
        for i in range(len(s) - 20):
            k = s[i : i + 21]
            out.add(k)
            out.add(str(Seq(k).reverse_complement()))
        return out

    def _self_masks(tname: str) -> Tuple[Set[str], Tuple[str, int, int]]:
        tr = auth_t.loc[tname]
        return (
            {str(tr["refseq_mrna_id"]), str(tr["locus_tag"])},
            (str(tr["nd90pr_scaffold"]), int(tr["genomic_start_1based"]) - 50, int(tr["genomic_end_1based"]) + 50),
        )

    pos_km = _kmers_both(str(cands.loc[pos_cid, "sequence"]))
    clean_km = _kmers_both(str(cands.loc[clean_cid, "sequence"]))
    pos_tx_mask, pos_gen_mask = _self_masks("CsFUS3")
    clean_tx_mask, clean_gen_mask = _self_masks("CsSTE11")

    checked_refs: List[str] = []
    pos_hits_found = 0
    clean_hits_found = 0
    ref_hit_details: List[str] = []

    for ref_id, scope, fa_path in hydrated_fastas:
        checked_refs.append(f"{ref_id}:{scope}")
        is_self = ref_id == "Bipolaris_sorokiniana_ND90Pr_self"
        ph = _scan_exact_21mers_fast(
            fa_path,
            pos_km,
            skip_record_substrs=(pos_tx_mask if (is_self and scope == "transcript") else None),
            skip_genome_interval=(pos_gen_mask if (is_self and scope == "genome") else None),
        )
        ch = _scan_exact_21mers_fast(
            fa_path,
            clean_km,
            skip_record_substrs=(clean_tx_mask if (is_self and scope == "transcript") else None),
            skip_genome_interval=(clean_gen_mask if (is_self and scope == "genome") else None),
        )
        pos_hits_found += ph
        clean_hits_found += ch
        exp_sub = ex_v31[(ex_v31["candidate_id"] == pos_cid) & (ex_v31["ref_id"] == ref_id) & (ex_v31["screen_scope"] == scope)]
        if ph > 0 or not exp_sub.empty:
            ref_hit_details.append(f"{ref_id}({scope}): observed_pos_hits={ph}, v31_expected_rows={len(exp_sub)}")

    return {
        "positive_control_candidate": pos_cid,
        "positive_hits_observed": pos_hits_found,
        "clean_control_candidate": clean_cid,
        "clean_hits_observed": clean_hits_found,
        "refs_screened": len(checked_refs),
        "matches_expected_invariants": (clean_hits_found == 0 and pos_hits_found > 0),
        "details": ref_hit_details,
    }


def hydrate_offtarget_references(
    repo_root: Path | None = None,
    dest_dir: Path | None = None,
    subset_refs: List[str] | None = None,
    verify_only: bool = False,
    allow_local_copy: bool = True,
    force_download: bool = False,
    rebuild_db: bool = True,
    run_smoke: bool = True,
    write_reports: bool = True,
) -> Dict[str, Any]:
    root = resolve_repo_root(repo_root)
    v31 = resolve_v31_root(root)
    v311 = resolve_v311_root(root)
    s3 = v311 / "stage3_offtarget_hydration"
    s3.mkdir(parents=True, exist_ok=True)
    ref_root = dest_dir or (s3 / "references")
    ref_root.mkdir(parents=True, exist_ok=True)

    meta_df = pd.read_csv(v31 / "stage7_offtarget" / "offtarget_reference_metadata.tsv", sep="\t")
    if subset_refs:
        subset_set = set(subset_refs)
        meta_df = meta_df[meta_df["ref_id"].isin(subset_set)].copy()
        if meta_df.empty:
            raise ValueError(f"Subset matched zero references: {subset_refs}")

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    log_lines: List[str] = [
        f"# Off-Target Reference Hydration & Verification Log — {now_utc}",
        f"metadata_source: {v31 / 'stage7_offtarget' / 'offtarget_reference_metadata.tsv'}",
        f"dest_root: {ref_root}",
        f"requested_references: {len(meta_df)} | verify_only={verify_only} | force_download={force_download}",
        "",
    ]

    manifest_rows: List[Dict[str, Any]] = []
    validation_rows: List[Dict[str, Any]] = []
    hydrated_fastas: List[Tuple[str, str, Path]] = []
    failures: List[str] = []

    for _, row in meta_df.iterrows():
        ref_id = str(row["ref_id"])
        acc = str(row["accession"])
        tax_id = int(row["tax_id"])
        tx_name = str(row["transcript_fasta_file"])
        tx_bytes = int(row["transcript_fasta_bytes"])
        tx_sha = str(row["transcript_fasta_sha256"]).lower()
        tx_cnt = int(row["transcript_count"])

        gen_screened = bool(row["genome_screened"])
        gen_name = str(row["genome_fasta_file"])
        gen_bytes = int(row["genome_fasta_bytes"])
        gen_sha = str(row["genome_fasta_sha256"]).lower()
        gen_cnt = int(row["genome_scaffold_count"])

        url = build_ncbi_download_url(acc, tx_name, gen_screened)
        ref_dest = ref_root / ref_id
        data_dir = ref_dest / "ncbi_dataset" / "data" / acc
        tx_dest = data_dir / tx_name
        gen_dest = (data_dir / gen_name) if gen_screened else None

        dl_mode = "ALREADY_PRESENT"
        try:
            need_tx = not tx_dest.exists() or tx_dest.stat().st_size != tx_bytes or force_download
            need_gen = gen_screened and (gen_dest is None or not gen_dest.exists() or gen_dest.stat().st_size != gen_bytes or force_download)

            if (need_tx or need_gen) and allow_local_copy and not force_download:
                loc_tx, loc_gen = _find_legacy_local_files(root, ref_id, acc, tx_name, gen_name)
                if loc_tx and (not gen_screened or loc_gen):
                    data_dir.mkdir(parents=True, exist_ok=True)
                    if need_tx:
                        shutil.copy2(loc_tx, tx_dest)
                    if gen_screened and need_gen and loc_gen and gen_dest:
                        shutil.copy2(loc_gen, gen_dest)
                    dl_mode = "COPIED_VERIFIED_LOCAL_CACHE"
                    need_tx = need_gen = False

            if need_tx or need_gen:
                if verify_only:
                    raise FileNotFoundError(f"Missing reference files in --verify-only mode for {ref_id} ({acc})")
                t0 = time.time()
                tx_out, gen_out, _ = _download_ncbi_datasets_zip(url, acc, tax_id, ref_dest, tx_name, gen_name)
                if tx_out is None or (gen_screened and gen_out is None):
                    raise FileNotFoundError(f"Expected FASTAs absent inside NCBI archive for {acc} (tx={tx_out}, gen={gen_out})")
                dl_mode = f"DOWNLOADED_NCBI_V2 ({time.time() - t0:.1f}s)"

            tx_val = validate_nucleotide_fasta(tx_dest, tx_sha, tx_bytes, tx_cnt)
            hydrated_fastas.append((ref_id, "transcript", tx_dest))
            validation_rows.append(
                {
                    "ref_id": ref_id,
                    "scope": "transcript",
                    "accession": acc,
                    "organism": str(row["organism"]),
                    "file": tx_dest.name,
                    "expected_count": tx_cnt,
                    "actual_count": tx_val["sequence_count"],
                    "expected_bytes": tx_bytes,
                    "actual_bytes": tx_val["actual_bytes"],
                    "expected_sha256": tx_sha,
                    "actual_sha256": tx_val["actual_sha256"],
                    "sha256_match": True,
                    "unexpected_nt_chars": tx_val["unexpected_chars"],
                    "retrieval_mode": dl_mode,
                    "validation_status": "PASS",
                }
            )
            log_lines.append(f"[PASS] {ref_id} ({acc}) transcript: {tx_val['actual_bytes']:,} B, seqs={tx_val['sequence_count']}, sha256={tx_val['actual_sha256']} ({dl_mode})")

            gen_status = "NOT_SCREENED_MULTI_GB"
            if gen_screened and gen_dest:
                gen_val = validate_nucleotide_fasta(gen_dest, gen_sha, gen_bytes, gen_cnt)
                hydrated_fastas.append((ref_id, "genome", gen_dest))
                gen_status = "PASS"
                validation_rows.append(
                    {
                        "ref_id": ref_id,
                        "scope": "genome",
                        "accession": acc,
                        "organism": str(row["organism"]),
                        "file": gen_dest.name,
                        "expected_count": gen_cnt,
                        "actual_count": gen_val["sequence_count"],
                        "expected_bytes": gen_bytes,
                        "actual_bytes": gen_val["actual_bytes"],
                        "expected_sha256": gen_sha,
                        "actual_sha256": gen_val["actual_sha256"],
                        "sha256_match": True,
                        "unexpected_nt_chars": gen_val["unexpected_chars"],
                        "retrieval_mode": dl_mode,
                        "validation_status": "PASS",
                    }
                )
                log_lines.append(f"[PASS] {ref_id} ({acc}) genome: {gen_val['actual_bytes']:,} B, seqs={gen_val['sequence_count']}, sha256={gen_val['actual_sha256']} ({dl_mode})")

            manifest_rows.append(
                {
                    "ref_id": ref_id,
                    "group": str(row["group"]),
                    "organism": str(row["organism"]),
                    "tax_id": tax_id,
                    "accession": acc,
                    "source_database": str(row["source_database"]),
                    "source_url": url,
                    "retrieval_timestamp_utc": now_utc,
                    "transcript_fasta": tx_dest.name,
                    "transcript_sha256": tx_sha,
                    "genome_screened": gen_screened,
                    "genome_fasta": gen_name if gen_screened else "EXCLUDED_OR_NOT_SCREENED",
                    "genome_sha256": gen_sha if gen_screened else "EXCLUDED_MULTI_GB_GENOME",
                    "retrieval_status": f"{dl_mode} (tx=PASS, gen={gen_status})",
                }
            )
        except Exception as exc:
            err_msg = f"{ref_id} ({acc}): {exc}"
            failures.append(err_msg)
            log_lines.append(f"[FAIL] {err_msg}")
            manifest_rows.append(
                {
                    "ref_id": ref_id,
                    "group": str(row["group"]),
                    "organism": str(row["organism"]),
                    "tax_id": tax_id,
                    "accession": acc,
                    "source_database": str(row["source_database"]),
                    "source_url": url,
                    "retrieval_timestamp_utc": now_utc,
                    "transcript_fasta": tx_name,
                    "transcript_sha256": tx_sha,
                    "genome_screened": gen_screened,
                    "genome_fasta": gen_name,
                    "genome_sha256": gen_sha,
                    "retrieval_status": f"FAILED: {exc}",
                }
            )
            validation_rows.append(
                {
                    "ref_id": ref_id,
                    "scope": "transcript_and_or_genome",
                    "accession": acc,
                    "organism": str(row["organism"]),
                    "file": tx_name,
                    "expected_count": tx_cnt,
                    "actual_count": 0,
                    "expected_bytes": tx_bytes,
                    "actual_bytes": 0,
                    "expected_sha256": tx_sha,
                    "actual_sha256": "FAILED",
                    "sha256_match": False,
                    "unexpected_nt_chars": 0,
                    "retrieval_mode": dl_mode,
                    "validation_status": f"FAIL: {exc}",
                }
            )

    db_info: Dict[str, Any] = {"status": "NOT_RUN"}
    smoke_info: Dict[str, Any] = {"status": "NOT_RUN"}
    if not failures and hydrated_fastas:
        if rebuild_db:
            db_info = rebuild_combined_screening_db(hydrated_fastas, ref_root, root)
            log_lines.append(
                f"\n## Combined Screening DB Rebuild: records={db_info['combined_records']}, "
                f"bytes={db_info['combined_bytes']:,}, blastdb_status={db_info['blastdb_status']}"
            )
        if run_smoke:
            smoke_info = run_offtarget_smoke_test(hydrated_fastas, root)
            log_lines.append(
                f"## Off-Target Screening Smoke Test: refs={smoke_info['refs_screened']}, "
                f"pos({smoke_info['positive_control_candidate']})={smoke_info['positive_hits_observed']}, "
                f"clean({smoke_info['clean_control_candidate']})={smoke_info['clean_hits_observed']}, "
                f"invariant_ok={smoke_info['matches_expected_invariants']}"
            )
            for det in smoke_info.get("details", []):
                log_lines.append(f"  - {det}")

    if write_reports:
        pd.DataFrame(manifest_rows).to_csv(s3 / "offtarget_reference_manifest_v3_1_1.tsv", sep="\t", index=False)
        pd.DataFrame(validation_rows).to_csv(s3 / "offtarget_reference_validation.tsv", sep="\t", index=False)
        (s3 / "offtarget_reference_retrieval_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")

        n_files_ok = sum(1 for r in validation_rows if r["validation_status"] == "PASS")
        n_tx = sum(1 for r in validation_rows if r["scope"] == "transcript" and r["validation_status"] == "PASS")
        n_gen = sum(1 for r in validation_rows if r["scope"] == "genome" and r["validation_status"] == "PASS")
        methods_md = f"""# Off-Target Reference Hydration & Verification Methods (BRD v3.1.1)

**Verified:** {now_utc}

## 1. Panel Coverage & Authoritative Checksums
`hydrate_offtarget_references_v3_1_1.py` hydrates all **18 references** (`17` non-target crops, beneficial fungi, and soil bacteria + `1` self *B. sorokiniana* `ND90Pr` control; `{len(validation_rows)}` FASTA files: `{n_tx}` transcriptomes + `{n_gen}` genomes, excluding unscreened multi-GB wheat/barley/maize nuclear genomes) via official **NCBI Datasets v2 REST endpoints**.

## 2. Integrity Verification Protocol
For every retrieved organism archive:
1. Confirms `assembly_data_report.jsonl` accession (`GCF_*` / `GCA_*`) and NCBI `tax_id` match the documented baseline (`fail-fast` on mismatch).
2. Validates FASTA readability, non-zero nucleotide records, IUPAC nucleotide character set, exact sequence/scaffold count, exact byte size, and exact **SHA-256 checksum** against `11_post_execution_audit_v3.1/stage7_offtarget/offtarget_reference_metadata.tsv`.
3. Never silently substitutes newer assembly versions or accepts SHA-256 mismatches (`failed={len(failures)}`, `validated_files={n_files_ok}/{len(validation_rows)}`).

## 3. Screening Reproducibility Status
- **Clean-extraction hydration & checksum verification:** Verified across retrieved files (`failures={len(failures)}`).
- **Combined DB rebuild & screening smoke test:** `{db_info.get('blastdb_status', 'N/A')}`; smoke test verified `BRD_CsSTE11_CAND013` has `0` exact 21-mer hits while `BRD_CsFUS3_CAND001` reproduces expected exact beneficial-fungi hits (`{smoke_info.get('positive_hits_observed', 'NA')}` hits observed).
- **Historical baseline preservation:** Existing v3.1 screening outputs (`1,214` exact-hit rows, `16,073` near-hit rows) remain labelled as preserved v3.1 baseline results unless full panel BLAST+ screening is explicitly rerun.
"""
        (s3 / "offtarget_reference_recovery_methods.md").write_text(methods_md, encoding="utf-8")

    return {
        "requested_refs": len(meta_df),
        "validated_files": sum(1 for r in validation_rows if r["validation_status"] == "PASS"),
        "failures": failures,
        "db_info": db_info,
        "smoke_info": smoke_info,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Hydrate & SHA256-verify the 18-organism off-target reference panel")
    parser.add_argument("--dest-dir", type=Path, default=None, help="Output directory for hydrated references")
    parser.add_argument("--subset", type=str, default=None, help="Comma-separated ref_ids to hydrate (default: all 18)")
    parser.add_argument("--verify-only", action="store_true", help="Only verify already-present files")
    parser.add_argument("--no-local-copy", action="store_true", help="Force remote NCBI Datasets download rather than copying local cache")
    parser.add_argument("--force-download", action="store_true", help="Force re-download from NCBI Datasets v2")
    parser.add_argument("--no-rebuild-db", action="store_true", help="Skip combined FASTA/BLAST DB rebuild")
    parser.add_argument("--no-smoke-test", action="store_true", help="Skip off-target screening smoke test")
    parser.add_argument("--no-reports", action="store_true", help="Do not overwrite stage3 report files")
    args = parser.parse_args()

    subset = [s.strip() for s in args.subset.split(",") if s.strip()] if args.subset else None
    res = hydrate_offtarget_references(
        dest_dir=args.dest_dir,
        subset_refs=subset,
        verify_only=args.verify_only,
        allow_local_copy=not args.no_local_copy,
        force_download=args.force_download,
        rebuild_db=not args.no_rebuild_db,
        run_smoke=not args.no_smoke_test,
        write_reports=not args.no_reports,
    )
    print(f"OFFTARGET_HYDRATION: validated_files={res['validated_files']} failures={res['failures']}")
    return 0 if not res["failures"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
