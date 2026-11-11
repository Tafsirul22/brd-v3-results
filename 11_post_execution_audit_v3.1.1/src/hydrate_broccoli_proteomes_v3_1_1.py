"""Stage 2 (FIX 2): Deterministic Recovery, Reconstruction & Validation of Broccoli v2.0.1 Proteome Inputs.

Reconstructs the exact historical proteome inputs omitted from 11_post_execution_audit_v3.1.zip:
  - Bsor_ND90Pr.faa          (from NCBI Datasets v2 GCF_000338995.1 PROT_FASTA -> 12,210 proteins)
  - Bsor_Sacc._Shoemaker.faa (from NCBI Datasets v2 GCA_013416765.1 PROT_FASTA -> 10,755 proteins)
  - Bcin_ASM83294v1.faa      (from pinned Brunxi/LUCID@37910bd6744b2ce6ce4c9b2582819b130e7ea509 Ensembl Bcin PEP -> 13,749 proteins)
  - all_proteomes_combined_v3.faa (exact combined FASTA with `{sp}||{id}` headers -> 36,714 proteins)

Validates format, non-empty sequences, unique identifiers, exact sequence counts, exact byte sizes,
and exact SHA-256 checksums against historical v3/v3.1 hashes, and optionally runs a 200-protein/species
Broccoli smoke test when `broccoli` + `diamond` executables are present.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import os
import shutil
import subprocess
import time
import urllib.request
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd
from Bio import SeqIO

try:
    from .paths_v3_1_1 import resolve_repo_root, resolve_tool, resolve_v311_root
except ImportError:
    from paths_v3_1_1 import resolve_repo_root, resolve_tool, resolve_v311_root

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUO*")

PROTEOME_SPECS: List[Dict[str, Any]] = [
    {
        "species_id": "Bcin_ASM83294v1",
        "organism": "Botrytis cinerea B05.10 (outgroup)",
        "accession": "ASM83294v1 (Ensembl Fungi / LUCID commit 37910bd6744b2ce6ce4c9b2582819b130e7ea509)",
        "annotation_release": "Ensembl Fungi ASM83294v1 (mirrored in Brunxi/LUCID@37910bd)",
        "source_type": "github_raw_pinned",
        "source_url": "https://raw.githubusercontent.com/Brunxi/LUCID/37910bd6744b2ce6ce4c9b2582819b130e7ea509/data/references/Botrytis_cinerea.ASM83294v1.pep.all.fa",
        "raw_filename": "Botrytis_cinerea.ASM83294v1.pep.all.fa",
        "expected_raw_sha256_alternatives": [
            "d897d558115de86c1122058d3b58c5cc386f53c1f78dda18376efec407571da7",  # LF raw from GitHub (64 hex chars)
            "a24542913ab0310ddfb0f45283aa8a46bf5828286cdf5d6965d6ec228d90a26d",  # CRLF git checkout on Windows
        ],
        "expected_protein_count": 13749,
        "staged_filename": "Bcin_ASM83294v1.faa",
        "expected_staged_bytes": 7034175,
        "expected_staged_sha256_crlf": "6fccfbf667b2c63695f8d2dbc603c24e413db53a47acdd2799c218cdd948ffa2",
        "redistribution_note": "Public Ensembl Fungi / GPL-free sequence annotation asset mirrored at pinned commit; hydrated dynamically via script.",
    },
    {
        "species_id": "Bsor_ND90Pr",
        "organism": "Bipolaris sorokiniana ND90Pr",
        "accession": "GCF_000338995.1",
        "annotation_release": "NCBI RefSeq Cocsa1 (2013-02-01)",
        "source_type": "ncbi_datasets_prot",
        "source_url": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000338995.1/download?include_annotation_type=PROT_FASTA",
        "raw_filename": "ncbi_dataset/data/GCF_000338995.1/protein.faa",
        "expected_raw_sha256_alternatives": [
            "cacc4be7fa1511474f949a7b0eec2ba863faf625ef3a8d49a3feacbabb79f708",
        ],
        "expected_protein_count": 12210,
        "staged_filename": "Bsor_ND90Pr.faa",
        "expected_staged_bytes": 5653448,
        "expected_staged_sha256_crlf": "2e8783a00d21843ca48f3b5ffc1fc2cf0b4abe486affb8a3cc887ec2acb54c39",
        "redistribution_note": "NCBI RefSeq Public Domain protein annotation; hydrated dynamically via NCBI Datasets v2 API.",
    },
    {
        "species_id": "Bsor_Sacc._Shoemaker",
        "organism": "Bipolaris sorokiniana Sacc. Shoemaker",
        "accession": "GCA_013416765.1",
        "annotation_release": "NCBI GenBank ASM1341676v1 (2020-06-29)",
        "source_type": "ncbi_datasets_prot",
        "source_url": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_013416765.1/download?include_annotation_type=PROT_FASTA",
        "raw_filename": "ncbi_dataset/data/GCA_013416765.1/protein.faa",
        "expected_raw_sha256_alternatives": [
            "5e4d2cf135d516d6f2ab627a0613d0b13365f32a6c181bfeea0e67ca988944b8",
        ],
        "expected_protein_count": 10755,
        "staged_filename": "Bsor_Sacc._Shoemaker.faa",
        "expected_staged_bytes": 5107801,
        "expected_staged_sha256_crlf": "583e29570f85cc85797b1b3a6773d6c454dcf9c5ad25d17269e6ccbd1cf740a3",
        "redistribution_note": "NCBI GenBank annotation; hydrated dynamically via NCBI Datasets v2 API.",
    },
]

EXPECTED_COMBINED_FILENAME = "all_proteomes_combined_v3.faa"
EXPECTED_COMBINED_COUNT = 36714
EXPECTED_COMBINED_BYTES = 18424497
EXPECTED_COMBINED_SHA256_CRLF = "b45191ae98961fef4f1e1c67d8486093f21cc392fc00864b88600a7763a39fd1"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def _download_raw_proteome(spec: Dict[str, Any], raw_dest: Path, timeout: int = 120) -> Tuple[Path, str]:
    raw_dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(spec["source_url"], headers={"User-Agent": "BRD-v3.1.1-ProteomeHydrator/1.0"})
    if spec["source_type"] == "github_raw_pinned":
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(raw_dest, "wb") as out_f:
            shutil.copyfileobj(resp, out_f)
    elif spec["source_type"] == "ncbi_datasets_prot":
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            payload = resp.read()
        with zipfile.ZipFile(io.BytesIO(payload)) as zf:
            target_member = next((m for m in zf.namelist() if m.endswith("protein.faa")), None)
            if not target_member:
                raise RuntimeError(f"protein.faa not found in NCBI Datasets ZIP for {spec['accession']}")
            with zf.open(target_member) as src_f, open(raw_dest, "wb") as dst_f:
                shutil.copyfileobj(src_f, dst_f)
    else:
        raise ValueError(f"Unsupported source_type: {spec['source_type']}")

    actual_raw_sha = sha256_file(raw_dest).lower()
    if actual_raw_sha not in [s.lower() for s in spec["expected_raw_sha256_alternatives"]]:
        raw_dest.unlink(missing_ok=True)
        raise RuntimeError(
            f"Raw SHA256 mismatch for {spec['species_id']} ({spec['accession']}): "
            f"got {actual_raw_sha}, expected one of {spec['expected_raw_sha256_alternatives']}"
        )
    return raw_dest, actual_raw_sha


def validate_protein_fasta(faa_path: Path, expected_count: int | None = None, require_species_prefix: str | None = None) -> Dict[str, Any]:
    if not faa_path.exists() or faa_path.stat().st_size == 0:
        raise RuntimeError(f"Protein FASTA missing or empty: {faa_path}")
    ids: List[str] = []
    empty_seqs = 0
    invalid_chars = 0
    total_aa = 0
    species_labels: Counter[str] = Counter()

    for rec in SeqIO.parse(str(faa_path), "fasta"):
        pid = rec.id.strip()
        seq = str(rec.seq).strip().upper()
        if not seq:
            empty_seqs += 1
        bad = set(seq) - VALID_AA
        if bad:
            invalid_chars += len(bad)
        total_aa += len(seq)
        ids.append(pid)
        if "||" in pid:
            species_labels[pid.split("||", 1)[0]] += 1

    if not ids:
        raise RuntimeError(f"FASTA contains zero parseable protein records: {faa_path}")
    dup_ids = sum(c - 1 for c in Counter(ids).values() if c > 1)
    if dup_ids > 0:
        raise RuntimeError(f"Duplicate protein identifiers detected ({dup_ids}) in {faa_path}")
    if empty_seqs > 0:
        raise RuntimeError(f"Empty protein sequences detected ({empty_seqs}) in {faa_path}")
    if expected_count is not None and len(ids) != expected_count:
        raise RuntimeError(f"Protein count mismatch for {faa_path}: got {len(ids)}, expected {expected_count}")
    if require_species_prefix is not None and require_species_prefix not in species_labels:
        raise RuntimeError(f"Missing expected species prefix {require_species_prefix!r} in {faa_path}")

    return {
        "protein_count": len(ids),
        "unique_id_count": len(set(ids)),
        "duplicate_id_count": dup_ids,
        "empty_sequence_count": empty_seqs,
        "unexpected_char_kinds": invalid_chars,
        "total_aa": total_aa,
        "species_prefix_counts": dict(species_labels),
    }


def run_broccoli_smoke_test(staged_faa_files: List[Path], out_dir: Path, root: Path) -> Dict[str, Any]:
    broccoli_bin = resolve_tool("broccoli", root)
    diamond_bin = resolve_tool("diamond", root)
    if not broccoli_bin or not diamond_bin:
        return {
            "status": "SKIPPED_BINARY_ABSENT",
            "broccoli_bin": broccoli_bin or "MISSING",
            "diamond_bin": diamond_bin or "MISSING",
            "elapsed_sec": 0.0,
            "og_count": 0,
            "pair_count": 0,
            "log_excerpt": "Broccoli/DIAMOND binary not found in extracted environment.",
        }

    smoke_root = out_dir / "_broccoli_smoke_v3_1_1"
    smoke_in = smoke_root / "in"
    smoke_out = smoke_root / "out"
    shutil.rmtree(smoke_root, ignore_errors=True)
    smoke_in.mkdir(parents=True, exist_ok=True)
    smoke_out.mkdir(parents=True, exist_ok=True)

    for faa in sorted(staged_faa_files):
        with open(smoke_in / faa.name, "w", encoding="utf-8", newline="\n") as fh:
            for idx, rec in enumerate(SeqIO.parse(str(faa), "fasta")):
                if idx >= 150:
                    break
                fh.write(f">{rec.id}\n{str(rec.seq)}\n")

    path_env = os.pathsep.join(
        [str(Path(diamond_bin).parent), str(Path(broccoli_bin).parent), os.environ.get("PATH", "")]
    )
    env = {**os.environ, "PATH": path_env}
    cmd = [broccoli_bin, "-dir", str(smoke_in), "-phylogenies", "nj", "-nb_hits", "6", "-threads", "4"]
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(smoke_out), env=env, capture_output=True, text=True)
    elapsed = round(time.time() - t0, 2)
    og_txt = smoke_out / "dir_step3" / "table_OGs_protein_names.txt"
    pr_txt = smoke_out / "dir_step4" / "orthologous_pairs.txt"
    og_cnt = sum(1 for ln in og_txt.read_text(encoding="utf-8").splitlines() if ln.strip() and not ln.startswith("#")) if og_txt.exists() else 0
    pr_cnt = sum(1 for ln in pr_txt.read_text(encoding="utf-8").splitlines() if ln.strip() and not ln.startswith("#")) if pr_txt.exists() else 0
    return {
        "status": "PASSED" if (proc.returncode == 0 and og_cnt > 0 and pr_cnt > 0) else f"FAILED_RC_{proc.returncode}",
        "broccoli_bin": broccoli_bin,
        "diamond_bin": diamond_bin,
        "elapsed_sec": elapsed,
        "og_count": og_cnt,
        "pair_count": pr_cnt,
        "log_excerpt": (proc.stdout + "\n" + proc.stderr).strip(),
    }


def hydrate_and_validate_proteomes(
    repo_root: Path | None = None,
    dest_dir: Path | None = None,
    force_download: bool = False,
    run_smoke: bool = False,
    write_reports: bool = True,
) -> Dict[str, Any]:
    root = resolve_repo_root(repo_root)
    v311 = resolve_v311_root(root)
    s2 = v311 / "stage2_proteome_hydration"
    s2.mkdir(parents=True, exist_ok=True)
    work_dir = dest_dir or (s2 / "hydrated_proteomes")
    raw_dir = work_dir / "raw_sources"
    staged_dir = work_dir / "input_proteomes"
    raw_dir.mkdir(parents=True, exist_ok=True)
    staged_dir.mkdir(parents=True, exist_ok=True)

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    log_lines: List[str] = [
        f"# Broccoli Proteome Recovery & Validation Log — {now_utc}",
        f"dest_dir: {work_dir}",
        "",
    ]

    manifest_rows: List[Dict[str, Any]] = []
    val_rows: List[Dict[str, Any]] = []
    combined_path = work_dir / EXPECTED_COMBINED_FILENAME
    staged_paths: List[Path] = []

    with open(combined_path, "w", encoding="utf-8", newline="\r\n") as comb_fh:
        for spec in sorted(PROTEOME_SPECS, key=lambda x: x["species_id"]):
            sp = spec["species_id"]
            raw_dest = raw_dir / f"{sp}_{Path(spec['raw_filename']).name}"
            if force_download or not raw_dest.exists() or sha256_file(raw_dest).lower() not in [
                s.lower() for s in spec["expected_raw_sha256_alternatives"]
            ]:
                _, raw_sha = _download_raw_proteome(spec, raw_dest)
                retrieval_mode = "DOWNLOADED_AND_VERIFIED"
            else:
                raw_sha = sha256_file(raw_dest).lower()
                retrieval_mode = "VERIFIED_EXISTING_RAW"

            staged_fp = staged_dir / spec["staged_filename"]
            staged_paths.append(staged_fp)
            with open(staged_fp, "w", encoding="utf-8", newline="\r\n") as sp_fh:
                for rec in SeqIO.parse(str(raw_dest), "fasta"):
                    seq_comb = str(rec.seq)
                    seq_staged = seq_comb.strip().upper()
                    comb_fh.write(f">{sp}||{rec.id}\n{seq_comb}\n")
                    sp_fh.write(f">{rec.id}\n{seq_staged}\n")

            v_info = validate_protein_fasta(staged_fp, expected_count=spec["expected_protein_count"])
            st_sha = sha256_file(staged_fp).lower()
            st_bytes = staged_fp.stat().st_size
            expected_st_sha = spec["expected_staged_sha256_crlf"].lower()
            sha_ok = st_sha == expected_st_sha
            if not sha_ok:
                raise RuntimeError(
                    f"Reconstructed {spec['staged_filename']} SHA-256 mismatch: {st_sha} != {expected_st_sha}"
                )

            log_lines.append(
                f"[{sp}] {retrieval_mode} raw_sha256={raw_sha} -> {spec['staged_filename']} "
                f"({st_bytes:,} B, seqs={v_info['protein_count']}, staged_sha256={st_sha}, exact_match={sha_ok})"
            )
            manifest_rows.append(
                {
                    "species_id": sp,
                    "organism": spec["organism"],
                    "accession": spec["accession"],
                    "annotation_release": spec["annotation_release"],
                    "source_type": spec["source_type"],
                    "source_url": spec["source_url"],
                    "retrieval_timestamp_utc": now_utc,
                    "raw_expected_sha256": spec["expected_raw_sha256_alternatives"][0],
                    "raw_actual_sha256": raw_sha,
                    "reconstructed_file": spec["staged_filename"],
                    "expected_bytes": spec["expected_staged_bytes"],
                    "actual_bytes": st_bytes,
                    "expected_sha256": expected_st_sha,
                    "actual_sha256": st_sha,
                    "redistribution_note": spec["redistribution_note"],
                }
            )
            val_rows.append(
                {
                    "file": spec["staged_filename"],
                    "species_or_scope": sp,
                    "protein_count": v_info["protein_count"],
                    "unique_id_count": v_info["unique_id_count"],
                    "duplicate_id_count": v_info["duplicate_id_count"],
                    "empty_sequence_count": v_info["empty_sequence_count"],
                    "unexpected_aa_chars": v_info["unexpected_char_kinds"],
                    "size_bytes": st_bytes,
                    "sha256": st_sha,
                    "matches_historical_sha256": sha_ok,
                    "validation_status": "PASS",
                }
            )

    # Validate reconstructed combined proteome FASTA
    comb_info = validate_protein_fasta(combined_path, expected_count=EXPECTED_COMBINED_COUNT)
    for sp_req, exp_cnt in ((s["species_id"], s["expected_protein_count"]) for s in PROTEOME_SPECS):
        got_cnt = comb_info["species_prefix_counts"].get(sp_req, 0)
        if got_cnt != exp_cnt:
            raise RuntimeError(f"Species {sp_req} count mismatch inside {EXPECTED_COMBINED_FILENAME}: {got_cnt} != {exp_cnt}")

    comb_sha = sha256_file(combined_path).lower()
    comb_bytes = combined_path.stat().st_size
    if comb_sha != EXPECTED_COMBINED_SHA256_CRLF:
        raise RuntimeError(
            f"{EXPECTED_COMBINED_FILENAME} SHA-256 mismatch: {comb_sha} != {EXPECTED_COMBINED_SHA256_CRLF}"
        )

    log_lines.append(
        f"[COMBINED] {EXPECTED_COMBINED_FILENAME} ({comb_bytes:,} B, seqs={comb_info['protein_count']}, "
        f"species={comb_info['species_prefix_counts']}, sha256={comb_sha}, exact_match=True)"
    )

    val_rows.append(
        {
            "file": EXPECTED_COMBINED_FILENAME,
            "species_or_scope": ";".join(f"{k}:{v}" for k, v in sorted(comb_info["species_prefix_counts"].items())),
            "protein_count": comb_info["protein_count"],
            "unique_id_count": comb_info["unique_id_count"],
            "duplicate_id_count": comb_info["duplicate_id_count"],
            "empty_sequence_count": comb_info["empty_sequence_count"],
            "unexpected_aa_chars": comb_info["unexpected_char_kinds"],
            "size_bytes": comb_bytes,
            "sha256": comb_sha,
            "matches_historical_sha256": True,
            "validation_status": "PASS",
        }
    )

    smoke_res: Dict[str, Any] = {"status": "NOT_REQUESTED"}
    if run_smoke:
        smoke_res = run_broccoli_smoke_test(staged_paths, work_dir, root)
        log_lines.extend(
            [
                "",
                f"## Broccoli Smoke Test (150 proteins/species) — status={smoke_res['status']} "
                f"elapsed={smoke_res['elapsed_sec']}s ogs={smoke_res['og_count']} pairs={smoke_res['pair_count']}",
                smoke_res["log_excerpt"],
            ]
        )

    if write_reports:
        pd.DataFrame(manifest_rows).to_csv(s2 / "broccoli_proteome_manifest_v3_1_1.tsv", sep="\t", index=False)
        pd.DataFrame(val_rows).to_csv(s2 / "broccoli_proteome_validation.tsv", sep="\t", index=False)
        (s2 / "broccoli_proteome_retrieval_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")

        methods_md = f"""# Broccoli Proteome Recovery & Validation Methods (BRD v3.1.1)

**Verified:** {now_utc}

## 1. Investigation Findings (Exact Historical Provenance)
Inspecting `11_post_execution_audit_v3/src/stage2_genome_qc_comparative.py` (lines 155–180) and `11_post_execution_audit_v3.1/src/stage5_broccoli_orthology.py` (lines 52–125) identified all three source proteomes behind the 36,714-protein Broccoli `v2.0.1` run:
1. **`Bsor_ND90Pr.faa` (`12,210` proteins):** Extracted from NCBI RefSeq `GCF_000338995.1` (`Cocsa1`, `protein.faa`, raw SHA-256 `cacc4be7fa1511474f949a7b0eec2ba863faf625ef3a8d49a3feacbabb79f708`).
2. **`Bsor_Sacc._Shoemaker.faa` (`10,755` proteins):** Extracted from NCBI GenBank `GCA_013416765.1` (`ASM1341676v1`, `protein.faa`, raw SHA-256 `5e4d2cf135d516d6f2ab627a0613d0b13365f32a6c181bfeea0e67ca988944b8`).
3. **`Bcin_ASM83294v1.faa` (`13,749` proteins):** Extracted from pinned upstream `Brunxi/LUCID` commit `37910bd6744b2ce6ce4c9b2582819b130e7ea509` (`data/references/Botrytis_cinerea.ASM83294v1.pep.all.fa`, Ensembl Fungi `ASM83294v1`, raw LF SHA-256 `d897d558115de86c1122058d3b58c5cc386f53c1f78dda18376efec407571da7`).

## 2. Deterministic Reconstruction
Writing `all_proteomes_combined_v3.faa` in sorted species order (`>{{species}}||{{protein_id}}\\n{{seq}}\\n`, CRLF line endings) and splitting per species with description-free `>{{protein_id}}\\n{{SEQ}}\\n` reproduces all four historical SHA-256 checksums **byte-for-byte (`100%`)**:
- `Bcin_ASM83294v1.faa`: `6fccfbf667b2c63695f8d2dbc603c24e413db53a47acdd2799c218cdd948ffa2` (`13,749` proteins, `7,034,175` B)
- `Bsor_ND90Pr.faa`: `2e8783a00d21843ca48f3b5ffc1fc2cf0b4abe486affb8a3cc887ec2acb54c39` (`12,210` proteins, `5,653,448` B)
- `Bsor_Sacc._Shoemaker.faa`: `583e29570f85cc85797b1b3a6773d6c454dcf9c5ad25d17269e6ccbd1cf740a3` (`10,755` proteins, `5,107,801` B)
- `all_proteomes_combined_v3.faa`: `b45191ae98961fef4f1e1c67d8486093f21cc392fc00864b88600a7763a39fd1` (`36,714` proteins, `18,424,497` B)

## 3. Broccoli Clean-Extraction Reproducibility Boundary
- **Proteome inputs:** `100%` reproducibly recoverable from public NCBI + pinned GitHub endpoints (`hydrate_broccoli_proteomes_v3_1_1.py`).
- **Executable requirement:** Re-running Broccoli requires `broccoli` (`v2.0.1`) + `diamond` (`v2.2.8`). Where installed (e.g., workspace `tools/bin/` or Bioconda), a clean-extraction smoke test runs (`status={smoke_res['status']}`); without `broccoli`/`diamond` binaries on a fresh system, orthology execution remains dependent on installing those external executables.
"""
        (s2 / "broccoli_proteome_recovery_methods.md").write_text(methods_md, encoding="utf-8")

    return {
        "combined_faa": combined_path,
        "staged_files": staged_paths,
        "validation": val_rows,
        "smoke_test": smoke_res,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Hydrate & SHA256-verify Broccoli v2.0.1 proteome inputs")
    parser.add_argument("--dest-dir", type=Path, default=None, help="Output directory for hydrated proteomes")
    parser.add_argument("--force-download", action="store_true", help="Re-download raw proteomes even if cached")
    parser.add_argument("--smoke-test", action="store_true", help="Run small Broccoli smoke test if binaries exist")
    parser.add_argument("--no-reports", action="store_true", help="Skip updating v3.1.1 stage2 report files")
    args = parser.parse_args()

    res = hydrate_and_validate_proteomes(
        dest_dir=args.dest_dir,
        force_download=args.force_download,
        run_smoke=args.smoke_test,
        write_reports=not args.no_reports,
    )
    print(f"PROTEOME_HYDRATION_OK: combined={res['combined_faa']} smoke={res['smoke_test'].get('status')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
