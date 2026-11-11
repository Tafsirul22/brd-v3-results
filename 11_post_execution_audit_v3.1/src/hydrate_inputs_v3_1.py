"""Stage 13: Reproducible Input Genome Downloader & SHA-256 Verifier for Clean Extractions.

Fixes v3 clean-archive failure AUDIT_05 where GCF_000338995.1_Cocsa1_genomic.fna and the
other 18 B. sorokiniana assemblies were omitted from the ZIP without a hydration workflow.
Uses standard Python stdlib (urllib + zipfile + hashlib) against the official NCBI Datasets
v2 REST endpoint and verifies exact byte sizes and SHA-256 checksums against
`11_post_execution_audit_v3.1/authoritative_data/genome_checksums_manifest_v3_1.tsv`.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Dict, List

import pandas as pd

try:
    from .paths_v3_1 import resolve_repo_root, resolve_v31_root
except ImportError:
    from paths_v3_1 import resolve_repo_root, resolve_v31_root


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def build_genome_manifest(repo_root: Path | None = None) -> Path:
    """Record exact accession, NCBI Datasets v2 URL, relative path, size, and SHA-256 for all 19 genomes."""
    root = resolve_repo_root(repo_root)
    v31 = resolve_v31_root(root)
    auth_dir = v31 / "authoritative_data"
    auth_dir.mkdir(parents=True, exist_ok=True)
    iso_df = pd.read_csv(root / "11_post_execution_audit_v3" / "authoritative_data" / "auth_isolates.tsv", sep="\t")

    rows: List[Dict[str, object]] = []
    for _, r in iso_df.iterrows():
        acc = str(r["assembly_accession"])
        rel = str(r["relative_genome_path"])
        fpath = root / rel
        if not fpath.exists():
            raise FileNotFoundError(f"Expected local genome missing while building manifest: {fpath}")
        rows.append(
            {
                "isolate_id": str(r["isolate_id"]),
                "assembly_accession": acc,
                "relative_genome_path": rel,
                "filename": fpath.name,
                "size_bytes": fpath.stat().st_size,
                "sha256": sha256_file(fpath),
                "ncbi_datasets_api_v2_url": (
                    f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
                    "?include_annotation_type=GENOME_FASTA"
                ),
                "is_primary_reference": acc == "GCF_000338995.1",
            }
        )
    out_tsv = auth_dir / "genome_checksums_manifest_v3_1.tsv"
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    return out_tsv


def hydrate_genomes(
    repo_root: Path | None = None,
    primary_only: bool = False,
    verify_existing: bool = True,
) -> Dict[str, object]:
    """Download missing B. sorokiniana genome FASTAs from NCBI Datasets v2 and verify SHA-256."""
    root = resolve_repo_root(repo_root)
    v31 = resolve_v31_root(root)
    man_path = v31 / "authoritative_data" / "genome_checksums_manifest_v3_1.tsv"
    man = pd.read_csv(man_path, sep="\t")
    if primary_only:
        man = man[man["is_primary_reference"].astype(bool)].copy()

    verified = 0
    downloaded = 0
    failed: List[str] = []

    for _, row in man.iterrows():
        rel = str(row["relative_genome_path"])
        dest = root / rel
        expected_sha = str(row["sha256"]).lower()
        expected_sz = int(row["size_bytes"])
        url = str(row["ncbi_datasets_api_v2_url"])
        target_name = str(row["filename"])

        if dest.exists() and verify_existing:
            if dest.stat().st_size == expected_sz and sha256_file(dest).lower() == expected_sha:
                verified += 1
                continue

        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            with tempfile.TemporaryDirectory() as tmpd:
                zip_tmp = Path(tmpd) / f"{row['assembly_accession']}.zip"
                req = urllib.request.Request(url, headers={"User-Agent": "BRD-v3.1-Hydrator/1.0"})
                with urllib.request.urlopen(req, timeout=180) as resp, open(zip_tmp, "wb") as out_f:
                    shutil.copyfileobj(resp, out_f)
                extracted_fna = None
                with zipfile.ZipFile(zip_tmp) as zf:
                    for member in zf.namelist():
                        if member.endswith(target_name) or (member.endswith(".fna") and "genomic" in member):
                            with zf.open(member) as src_f, open(dest, "wb") as dst_f:
                                shutil.copyfileobj(src_f, dst_f)
                            extracted_fna = dest
                            break
                if extracted_fna is None or not dest.exists():
                    failed.append(f"{row['assembly_accession']}: no genomic .fna inside downloaded archive")
                    continue
                actual_sha = sha256_file(dest).lower()
                if actual_sha != expected_sha:
                    failed.append(f"{row['assembly_accession']}: SHA256 mismatch ({actual_sha} != {expected_sha})")
                    dest.unlink(missing_ok=True)
                    continue
                downloaded += 1
                verified += 1
        except Exception as exc:
            failed.append(f"{row['assembly_accession']}: {exc}")

    return {"total_requested": len(man), "verified": verified, "downloaded": downloaded, "failed": failed}


def main() -> int:
    parser = argparse.ArgumentParser(description="Hydrate and SHA256-verify B. sorokiniana genome FASTAs")
    parser.add_argument("--build-manifest", action="store_true", help="Rebuild genome_checksums_manifest_v3_1.tsv from local inputs")
    parser.add_argument("--primary-only", action="store_true", help="Only fetch/verify GCF_000338995.1_Cocsa1_genomic.fna")
    args = parser.parse_args()

    if args.build_manifest:
        p = build_genome_manifest()
        print(f"MANIFEST_BUILT: {p}")
        return 0

    res = hydrate_genomes(primary_only=args.primary_only)
    print(f"HYDRATION_RESULT: {res}")
    return 0 if not res["failed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
