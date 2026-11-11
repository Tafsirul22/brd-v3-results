"""Stage 0: Inspect supplied 11_post_execution_audit_v3.1.zip, extract clean baseline, compare vs repo, and snapshot rankings."""

from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

try:
    from .paths_v3_1_1 import resolve_repo_root, resolve_v311_root
except ImportError:
    from paths_v3_1_1 import resolve_repo_root, resolve_v311_root


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def run_stage0_v311(repo_root: Path | None = None) -> Dict[str, Any]:
    root = resolve_repo_root(repo_root)
    v311 = resolve_v311_root(root)
    s0 = v311 / "stage0_audit"
    s0.mkdir(parents=True, exist_ok=True)
    extract_dir = v311 / "_v31_extracted_baseline"
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    v31_zip = root / "11_post_execution_audit_v3.1.zip"
    if not v31_zip.exists():
        raise FileNotFoundError(f"Supplied v3.1 release archive not found: {v31_zip}")

    zip_sha256 = sha256_file(v31_zip)
    zip_bytes = v31_zip.stat().st_size

    with zipfile.ZipFile(v31_zip, "r") as zf:
        bad_member = zf.testzip()
        if bad_member is not None:
            raise RuntimeError(f"v3.1 archive CRC failure on member: {bad_member}")
        infos = sorted(zf.infolist(), key=lambda x: x.filename)

        if extract_dir.exists():
            shutil.rmtree(extract_dir)
        extract_dir.mkdir(parents=True, exist_ok=True)
        zf.extractall(extract_dir)

    inv_rows: List[Dict[str, Any]] = []
    cmp_rows: List[Dict[str, Any]] = []
    for info in infos:
        if info.is_dir():
            continue
        rel = info.filename
        ext_fp = extract_dir / rel
        ext_sha = sha256_file(ext_fp)
        inv_rows.append(
            {
                "relative_path": rel,
                "uncompressed_bytes": info.file_size,
                "compressed_bytes": info.compress_size,
                "crc32": f"{info.CRC:08x}",
                "sha256": ext_sha,
            }
        )
        repo_fp = root / rel
        if repo_fp.exists() and repo_fp.is_file():
            repo_sha = sha256_file(repo_fp)
            cmp_status = "IDENTICAL" if repo_sha == ext_sha else "MODIFIED_IN_REPO"
        else:
            repo_sha = ""
            cmp_status = "MISSING_IN_REPO"
        cmp_rows.append(
            {
                "relative_path": rel,
                "zip_sha256": ext_sha,
                "repo_sha256": repo_sha,
                "comparison_status": cmp_status,
            }
        )

    # Also check repo files under 11_post_execution_audit_v3.1 absent from v3.1 ZIP
    zip_names = {r["relative_path"] for r in inv_rows}
    v31_dir = root / "11_post_execution_audit_v3.1"
    if v31_dir.exists():
        for fp in sorted(v31_dir.rglob("*")):
            if fp.is_file():
                rel = fp.relative_to(root).as_posix()
                if rel not in zip_names and "__pycache__" not in fp.parts and ".pytest_cache" not in fp.parts:
                    cmp_rows.append(
                        {
                            "relative_path": rel,
                            "zip_sha256": "",
                            "repo_sha256": sha256_file(fp),
                            "comparison_status": "OMITTED_FROM_V3_1_ZIP",
                        }
                    )

    pd.DataFrame(inv_rows).to_csv(s0 / "v3_1_archive_inventory.tsv", sep="\t", index=False)
    cmp_df = pd.DataFrame(cmp_rows)
    cmp_df.to_csv(s0 / "v3_1_vs_repo_comparison.tsv", sep="\t", index=False)

    # Baseline snapshot of v3.1 outputs, rankings, and key historical scripts (from extracted ZIP!)
    v31_ext = extract_dir / "11_post_execution_audit_v3.1"
    reg_rk = pd.read_csv(v31_ext / "stage10_ranking" / "independent_region_ranking_v3_1.tsv", sep="\t")
    cand_rk = pd.read_csv(v31_ext / "stage10_ranking" / "candidate_ranking_v3_1.tsv", sep="\t")
    mast = pd.read_csv(v31_ext / "stage8_candidate_design" / "candidate_master_v3_1.tsv", sep="\t")

    key_files_to_hash = [
        "11_post_execution_audit_v3.1.zip",
        "11_post_execution_audit/tests/test_audit_v2.py",
        "11_post_execution_audit_v3/tests/test_brd_v3.py",
        "tests/test_brd.py",
        "11_post_execution_audit_v3.1/tests/test_brd_v3_1.py",
        "11_post_execution_audit_v3.1/hashes_v3_1.tsv",
        "11_post_execution_audit_v3.1/V3_1_Final_Validation_Report.md",
        "11_post_execution_audit_v3.1/V3_1_Change_Log.md",
    ]
    baseline_hashes = {
        rel: (sha256_file(root / rel) if (root / rel).exists() else "MISSING")
        for rel in key_files_to_hash
    }

    snapshot = {
        "inspected_utc": now_utc,
        "v31_archive_path": "11_post_execution_audit_v3.1.zip",
        "v31_archive_bytes": zip_bytes,
        "v31_archive_sha256": zip_sha256,
        "v31_archive_crc_test": "PASS (None)",
        "v31_archive_file_count": len(inv_rows),
        "baseline_hashes": baseline_hashes,
        "clean_shortlist_8": reg_rk[reg_rk["passes_zero_exact_21mer_filter"]][
            [
                "independent_region_rank_v3_1",
                "composite_index_rank_overall_12regions",
                "pareto_front_3d_regions",
                "region_id",
                "candidate_id",
                "target",
                "heuristic_composite_prioritization_index",
            ]
        ].to_dict(orient="records"),
        "candidate_master_coord_counts": {
            "total_candidates": int(len(mast)),
            "is_exact_cds_substring_true": int(mast["is_exact_cds_substring"].astype(bool).sum()),
            "is_exact_cds_substring_false": int((~mast["is_exact_cds_substring"].astype(bool)).sum()),
            "synthetic_ids": mast.loc[~mast["is_exact_cds_substring"].astype(bool), "candidate_id"].tolist(),
        },
        "top5_all_candidates": cand_rk.head(5)[
            ["composite_index_rank_all50", "candidate_id", "target", "heuristic_composite_prioritization_index"]
        ].to_dict(orient="records"),
    }
    (s0 / "v3_1_baseline_snapshot.json").write_text(json.dumps(snapshot, indent=2), encoding="utf-8")

    # 6 remediation defect register
    defects = [
        {"fix_id": "FIX_1", "area": "Optional Dependency", "defect_in_v31": "Unconditional `import primer3` in `brd/integrations/lucid.py` (and `m07_dsRNA_design.py`, `m10_reporting.py`) crashes in clean env when `primer3-py` absent", "target_resolution": "Guard optional imports, expose PRIMER3_AVAILABLE, raise actionable RuntimeError when Primer3 called without package, test both installed & absent states"},
        {"fix_id": "FIX_2", "area": "Broccoli Proteomes", "defect_in_v31": "Broccoli inputs (`Bsor_ND90Pr.faa`, `Bsor_Sacc._Shoemaker.faa`, `Bcin_ASM83294v1.faa`, `all_proteomes_combined_v3.faa`) excluded without hydration script", "target_resolution": "Implement `hydrate_broccoli_proteomes_v3_1_1.py` fetching exact NCBI/GitHub sources, reconstructing identical files, verifying SHA-256 & Broccoli smoke test"},
        {"fix_id": "FIX_3", "area": "Off-Target Panel", "defect_in_v31": "Metadata & SHA-256 recorded for 18 panel refs (`17` non-targets + self), but no reference hydration script existed", "target_resolution": "Implement `hydrate_offtarget_references_v3_1_1.py` fetching exact versioned NCBI accessions, verifying length & SHA-256, rebuilding DB & running smoke screen"},
        {"fix_id": "FIX_4", "area": "Regression Testing", "defect_in_v31": "v3.1 ZIP lacked `test_audit_v2.py` (`test_audit_v2.py` also hardcoded `OUT_DIR`), v2/v1 files, and complete clean-extraction test accounting", "target_resolution": "Restore historical tests with relative paths, include required non-huge v2/v1 fixtures, document core vs optional deps & passed vs blocked counts"},
        {"fix_id": "FIX_5", "area": "Manifest Completeness", "defect_in_v31": "`hashes_v3_1.tsv` omitted `V3_1_Execution_Log.txt`, `V3_1_Test_Results.txt`, and `V3_1_Test_Summary.md`", "target_resolution": "Freeze all files & reports before hashing, include ALL release files in `hashes_v3_1_1.tsv`, and verify via standalone script detecting unlisted/missing files"},
        {"fix_id": "FIX_6", "area": "Coordinate Consistency", "defect_in_v31": "`V3_1_Change_Log.md` claimed uniform verification of all 50 candidate coordinates on target FASTA, contradicting report (`49/50` contiguous CDS + `1` synthetic construct `BRD_CsSTE11_CAND009`)", "target_resolution": "Distinguish recorded coordinate assignments (`50`), contiguous CDS (`49`), stitched constructs (`1`: `CAND009`), and unresolved (`0`) consistently across all artifacts"},
    ]
    pd.DataFrame(defects).to_csv(s0 / "v3_1_remediation_register.tsv", sep="\t", index=False)

    n_ident = int((cmp_df["comparison_status"] == "IDENTICAL").sum())
    n_mod = int((cmp_df["comparison_status"] == "MODIFIED_IN_REPO").sum())
    n_omit = int((cmp_df["comparison_status"] == "OMITTED_FROM_V3_1_ZIP").sum())
    report_md = f"""# Stage 0 — Supplied v3.1 Archive Inspection & Baseline Snapshot (BRD v3.1.1)

**Inspected:** {now_utc}  
**Authoritative Archive:** `11_post_execution_audit_v3.1.zip` (`{zip_bytes:,}` bytes, SHA-256 `{zip_sha256}`, `testzip() = None` across `{len(inv_rows)}` entries)

## Comparison vs Workspace Repository
- `IDENTICAL` entries in extracted archive vs repo: `{n_ident}`
- `MODIFIED_IN_REPO` (since archive creation — expected only for Stage 1 fixes applied on `v3.1.1`): `{n_mod}`
- `OMITTED_FROM_V3_1_ZIP` (`input_proteomes/`, `dir_step1/*.fas`, `dir_step2/`, caches): `{n_omit}`

## Preserved v3.1 Baseline Shortlist (8 Clean Independent Regions)
| Rank | Candidate | Target | Region | Score |
|---:|---|---|---|---:|
""" + "\n".join(
        f"| {r['independent_region_rank_v3_1']} | {r['candidate_id']} | {r['target']} | {r['region_id']} | {r['heuristic_composite_prioritization_index']} |"
        for r in snapshot["clean_shortlist_8"]
    ) + "\n"
    (s0 / "Stage0_V3_1_Inspection_Report.md").write_text(report_md, encoding="utf-8")
    return {"snapshot": s0 / "v3_1_baseline_snapshot.json", "inventory": s0 / "v3_1_archive_inventory.tsv", "zip_sha256": zip_sha256}


if __name__ == "__main__":
    out = run_stage0_v311()
    print("STAGE0_V311_DONE:", out)
