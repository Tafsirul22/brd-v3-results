"""Stage 5 (FIX 5): Complete Release SHA-256 Manifest Builder & Independent Verifier.

Builds `11_post_execution_audit_v3.1.1/hashes_v3_1_1.tsv` covering every packaged release file
(excluding `hashes_v3_1_1.tsv` itself), documents intentionally excluded files in
`excluded_large_files_v3_1_1.tsv`, and provides strict independent verification returning exit code != 0
on any missing file, size mismatch, SHA-256 mismatch, duplicate path, or unlisted file.
"""

from __future__ import annotations

import argparse
import hashlib
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

try:
    from .paths_v3_1_1 import resolve_repo_root, resolve_v311_root
    from .release_layout_v3_1_1 import (
        V31_EXCL_SUBS,
        V311_EXCL_PREFIXES,
        V3_EXCL_SUBS,
        classify_file_category,
        is_included_release_rel,
        iter_release_files,
    )
except ImportError:
    from paths_v3_1_1 import resolve_repo_root, resolve_v311_root
    from release_layout_v3_1_1 import (
        V31_EXCL_SUBS,
        V311_EXCL_PREFIXES,
        V3_EXCL_SUBS,
        classify_file_category,
        is_included_release_rel,
        iter_release_files,
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def build_excluded_manifest(repo_root: Path | None = None) -> Path:
    root = resolve_repo_root(repo_root)
    v311 = resolve_v311_root(root)
    rows: List[Dict[str, Any]] = []

    def _add(rel: str, size: int, reason: str) -> None:
        rows.append({"relative_path": rel, "size_bytes": size, "exclusion_reason": reason})

    # Document excluded genome assemblies, proteomes, off-target refs, and large intermediate caches
    g_man = root / "11_post_execution_audit_v3.1" / "authoritative_data" / "genome_checksums_manifest_v3_1.tsv"
    if g_man.exists():
        gdf = pd.read_csv(g_man, sep="\t")
        for _, gr in gdf.iterrows():
            _add(
                str(gr["relative_genome_path"]),
                int(gr["size_bytes"]),
                f"Large NCBI genome FASTA omitted from compact ZIP; recoverable via hydrate_inputs_v3_1.py (sha256={gr['sha256']})",
            )

    p_man = v311 / "stage2_proteome_hydration" / "broccoli_proteome_manifest_v3_1_1.tsv"
    if p_man.exists():
        for _, pr in pd.read_csv(p_man, sep="\t").iterrows():
            _add(
                f"11_post_execution_audit_v3.1.1/stage2_proteome_hydration/hydrated_proteomes/input_proteomes/{pr['reconstructed_file']}",
                int(pr["expected_bytes"]),
                f"Broccoli input proteome omitted from compact ZIP; recoverable identically via hydrate_broccoli_proteomes_v3_1_1.py (sha256={pr['expected_sha256']})",
            )
    _add(
        "11_post_execution_audit_v3.1.1/stage2_proteome_hydration/hydrated_proteomes/all_proteomes_combined_v3.faa",
        18424497,
        "Combined 3-species proteome FASTA (36,714 seqs) omitted from compact ZIP; recoverable identically via hydrate_broccoli_proteomes_v3_1_1.py (sha256=b45191ae98961fef4f1e1c67d8486093f21cc392fc00864b88600a7763a39fd1)",
    )

    ot_meta = root / "11_post_execution_audit_v3.1" / "stage7_offtarget" / "offtarget_reference_metadata.tsv"
    if ot_meta.exists():
        for _, orow in pd.read_csv(ot_meta, sep="\t").iterrows():
            _add(
                f"11_post_execution_audit_v3.1.1/stage3_offtarget_hydration/references/{orow['ref_id']}/ncbi_dataset/data/{orow['accession']}/{orow['transcript_fasta_file']}",
                int(orow["transcript_fasta_bytes"]),
                f"Multi-100MB off-target transcriptome FASTA omitted from compact ZIP; recoverable via hydrate_offtarget_references_v3_1_1.py (sha256={orow['transcript_fasta_sha256']})",
            )
            if bool(orow["genome_screened"]):
                _add(
                    f"11_post_execution_audit_v3.1.1/stage3_offtarget_hydration/references/{orow['ref_id']}/ncbi_dataset/data/{orow['accession']}/{orow['genome_fasta_file']}",
                    int(orow["genome_fasta_bytes"]),
                    f"Large off-target genome FASTA omitted from compact ZIP; recoverable via hydrate_offtarget_references_v3_1_1.py (sha256={orow['genome_fasta_sha256']})",
                )

    for rel, fp in [
        ("11_post_execution_audit_v3.1/hashes_v3_1.tsv", None),
    ]:
        pass

    for d_root, subs, reason in [
        (root / "11_post_execution_audit_v3", V3_EXCL_SUBS, "Large v3 intermediate BLAST/DIAMOND/source cache excluded from compact ZIP"),
        (root / "11_post_execution_audit_v3.1", V31_EXCL_SUBS, "Large v3.1 Broccoli intermediate (.faa/.fas/dir_step2) excluded from compact ZIP"),
        (v311, V311_EXCL_PREFIXES, "Stage 0 baseline extraction snapshot / hydrated FASTA workspace excluded from compact ZIP"),
    ]:
        if d_root.exists():
            for fp in sorted(d_root.rglob("*")):
                if fp.is_file() and "__pycache__" not in fp.parts and ".pytest_cache" not in fp.parts:
                    rel_d = fp.relative_to(d_root).as_posix()
                    if any(rel_d.startswith(s) or s in rel_d for s in subs):
                        _add(fp.relative_to(root).as_posix(), fp.stat().st_size, reason)

    out_path = v311 / "excluded_large_files_v3_1_1.tsv"
    df = pd.DataFrame(rows).drop_duplicates(subset=["relative_path"]).sort_values("relative_path")
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


def generate_hashes_v3_1_1(repo_root: Path | None = None) -> Path:
    """Build hashes_v3_1_1.tsv covering EVERY release file except hashes_v3_1_1.tsv itself."""
    root = resolve_repo_root(repo_root)
    v311 = resolve_v311_root(root)
    build_excluded_manifest(root)

    rows: List[Dict[str, Any]] = []
    for rel, fp in iter_release_files(root, include_manifest_itself=False):
        rows.append(
            {
                "relative_path": rel,
                "file_size_bytes": fp.stat().st_size,
                "sha256": sha256_file(fp),
                "file_category": classify_file_category(rel),
                "included_in_release": True,
            }
        )

    out_tsv = v311 / "hashes_v3_1_1.tsv"
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    return out_tsv


def verify_manifest(
    repo_root: Path | None = None,
    manifest_path: Path | None = None,
    zip_path: Path | None = None,
    check_unlisted_on_disk: bool = True,
) -> Dict[str, Any]:
    """Independently verify manifest entries against disk and/or a built ZIP archive."""
    root = resolve_repo_root(repo_root)
    v311 = resolve_v311_root(root)
    mpath = manifest_path or (v311 / "hashes_v3_1_1.tsv")
    if not mpath.exists():
        raise FileNotFoundError(f"Manifest missing: {mpath}")

    df = pd.read_csv(mpath, sep="\t", dtype={"file_size_bytes": "int64"})
    req_cols = {"relative_path", "file_size_bytes", "sha256", "file_category", "included_in_release"}
    if not req_cols.issubset(set(df.columns)):
        raise ValueError(f"Manifest missing required columns: {sorted(req_cols - set(df.columns))}")

    rel_list = [str(x) for x in df["relative_path"]]
    counts = Counter(rel_list)
    duplicates = sorted([k for k, v in counts.items() if v > 1])

    self_included = [r for r in rel_list if r == "11_post_execution_audit_v3.1.1/hashes_v3_1_1.tsv"]

    missing_on_disk: List[str] = []
    size_mismatches_disk: List[str] = []
    sha_mismatches_disk: List[str] = []
    verified_disk = 0

    for _, row in df[df["included_in_release"].astype(bool)].iterrows():
        rel = str(row["relative_path"])
        exp_sz = int(row["file_size_bytes"])
        exp_sha = str(row["sha256"]).lower()
        fp = root / rel
        if not fp.is_file():
            missing_on_disk.append(rel)
            continue
        act_sz = fp.stat().st_size
        act_sha = sha256_file(fp).lower()
        if act_sz != exp_sz:
            size_mismatches_disk.append(f"{rel} (got {act_sz} != {exp_sz})")
        elif act_sha != exp_sha:
            sha_mismatches_disk.append(f"{rel} (got {act_sha} != {exp_sha})")
        else:
            verified_disk += 1

    unlisted_on_disk: List[str] = []
    if check_unlisted_on_disk:
        expected_set = set(df[df["included_in_release"].astype(bool)]["relative_path"])
        for rel, _ in iter_release_files(root, include_manifest_itself=False):
            if rel not in expected_set:
                unlisted_on_disk.append(rel)

    zip_missing: List[str] = []
    zip_unlisted: List[str] = []
    zip_mismatches: List[str] = []
    verified_zip = 0
    if zip_path is not None:
        if not zip_path.exists():
            raise FileNotFoundError(f"ZIP archive missing: {zip_path}")
        with zipfile.ZipFile(zip_path, "r") as zf:
            bad_crc = zf.testzip()
            if bad_crc is not None:
                zip_mismatches.append(f"CRC failure: {bad_crc}")
            zip_entries = {n for n in zf.namelist() if not n.endswith("/")}
            expected_with_manifest = set(df[df["included_in_release"].astype(bool)]["relative_path"]) | {
                "11_post_execution_audit_v3.1.1/hashes_v3_1_1.tsv"
            }
            zip_missing = sorted(expected_with_manifest - zip_entries)
            zip_unlisted = sorted(zip_entries - expected_with_manifest)
            for _, row in df[df["included_in_release"].astype(bool)].iterrows():
                rel = str(row["relative_path"])
                if rel not in zip_entries:
                    continue
                data = zf.read(rel)
                if len(data) != int(row["file_size_bytes"]):
                    zip_mismatches.append(f"ZIP size mismatch: {rel}")
                elif sha256_bytes(data).lower() != str(row["sha256"]).lower():
                    zip_mismatches.append(f"ZIP sha256 mismatch: {rel}")
                else:
                    verified_zip += 1

    ok = (
        not duplicates
        and not self_included
        and not missing_on_disk
        and not size_mismatches_disk
        and not sha_mismatches_disk
        and not unlisted_on_disk
        and not zip_missing
        and not zip_unlisted
        and not zip_mismatches
    )
    return {
        "ok": ok,
        "manifest_entries": len(df),
        "verified_disk": verified_disk,
        "duplicates": duplicates,
        "self_included": self_included,
        "missing_on_disk": missing_on_disk,
        "size_mismatches_disk": size_mismatches_disk,
        "sha_mismatches_disk": sha_mismatches_disk,
        "unlisted_on_disk": unlisted_on_disk,
        "verified_zip": verified_zip,
        "zip_missing": zip_missing,
        "zip_unlisted": zip_unlisted,
        "zip_mismatches": zip_mismatches,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build or independently verify hashes_v3_1_1.tsv")
    parser.add_argument("--build", action="store_true", help="Generate hashes_v3_1_1.tsv and excluded_large_files_v3_1_1.tsv")
    parser.add_argument("--zip-file", type=Path, default=None, help="Also verify a release ZIP archive against hashes_v3_1_1.tsv")
    parser.add_argument("--no-unlisted-disk-check", action="store_true", help="Skip unlisted release file scan on disk")
    args = parser.parse_args()

    if args.build:
        p = generate_hashes_v3_1_1()
        print(f"MANIFEST_GENERATED: {p}")

    res = verify_manifest(zip_path=args.zip_file, check_unlisted_on_disk=not args.no_unlisted_disk_check)
    print(f"MANIFEST_VERIFY: {res}")
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
