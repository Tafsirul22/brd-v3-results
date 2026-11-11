"""Canonical release file enumeration and exclusion definitions for BRD v3.1.1.

Shared identically by `hashes_v3_1_1.tsv`, `verify_manifest_v3_1_1.py`, and
`build_release_zip_v311()` so packaged ZIP entries and SHA-256 manifest entries stay 1:1 synchronized.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

SECRET_NAME_TOKENS = ("secret", "oauth", "token", "credential", ".env")

V3_EXCL_SUBS = (
    "expanded_references/sources/",
    "expanded_references/blast_dbs/",
    "expanded_references/raw_blast/",
    "orthology_diamond_rbh/all_proteomes_combined",
    "orthology_diamond_rbh/all_proteomes_db",
    "orthology_diamond_rbh/diamond_all_vs_all",
    "orthology_diamond_rbh/GCA_",
)

V31_EXCL_SUBS = (
    "stage5_broccoli_orthology/input_proteomes/",
    "stage5_broccoli_orthology/broccoli_run/dir_step1/0.fas",
    "stage5_broccoli_orthology/broccoli_run/dir_step1/1.fas",
    "stage5_broccoli_orthology/broccoli_run/dir_step1/2.fas",
    "stage5_broccoli_orthology/broccoli_run/dir_step2/",
)

V311_EXCL_PREFIXES = (
    "_v31_extracted_baseline/",
    "stage2_proteome_hydration/hydrated_proteomes/",
    "stage3_offtarget_hydration/references/",
    "_tmp",
)


def _has_runtime_cache_part(fp: Path) -> bool:
    return any(p in ("__pycache__", ".pytest_cache", ".git") for p in fp.parts)


def _has_secret_token(fp: Path) -> bool:
    return any(tok in fp.name.lower() for tok in SECRET_NAME_TOKENS)


def classify_file_category(rel: str) -> str:
    if rel.endswith(".py"):
        if "/tests/" in rel or rel.startswith("tests/"):
            return "test_code"
        return "source_code"
    if rel.endswith((".tsv", ".csv", ".json")):
        return "data_or_manifest_table"
    if rel.endswith((".md", ".html", ".pdf", ".txt")):
        if "Log" in rel or "log" in rel:
            return "execution_or_retrieval_log"
        if "Test" in rel or "test" in rel:
            return "test_report"
        return "report_or_documentation"
    if rel.endswith((".fasta", ".fa", ".fna", ".faa")):
        return "sequence_data"
    if rel.endswith((".png", ".svg")):
        return "figure"
    if rel.endswith((".yml", ".yaml")) or rel in (".gitignore", "requirements.txt", "requirements-optional.txt") or rel.endswith(("requirements.txt", "requirements-optional.txt", "environment.yml", "Dockerfile")):
        return "config_or_dependency_spec"
    if ".licenses/" in rel or rel.startswith(".licenses/"):
        return "license_notice"
    return "supporting_artifact"


def is_included_release_rel(rel: str) -> bool:
    """Return True if a normalized POSIX relative path belongs inside the v3.1.1 release archive."""
    norm = rel.lstrip("./")
    if "/__pycache__/" in norm or "/.pytest_cache/" in norm or norm.startswith(("__pycache__/", ".pytest_cache/", ".git/")):
        return False
    if any(tok in Path(norm).name.lower() for tok in SECRET_NAME_TOKENS):
        return False

    exact_root_files = {
        ".gitignore",
        "README.md",
        "requirements.txt",
        "requirements-optional.txt",
        "environment.yml",
        "02_genome_qc/assembly_quality_summary.tsv",
        "08_off_target_screening/blast_databases/database_manifest.tsv",
        "10_final_reports/BRD_Target_Gene_Sequences.fasta",
        "10_final_reports/BRD_Target_Gene_Summary.tsv",
        "10_final_reports/BRD_dsRNA_Candidates.tsv",
    }
    if norm in exact_root_files:
        return True
    if norm.startswith(("brd/", "config/", "tests/", ".licenses/", "11_post_execution_audit/")):
        return True
    if norm.startswith("11_post_execution_audit_v3/"):
        sub = norm[len("11_post_execution_audit_v3/") :]
        return not any(ex in sub for ex in V3_EXCL_SUBS)
    if norm.startswith("11_post_execution_audit_v3.1/"):
        sub = norm[len("11_post_execution_audit_v3.1/") :]
        return not any(sub.startswith(ex) or f"/{ex}" in sub for ex in V31_EXCL_SUBS)
    if norm.startswith("11_post_execution_audit_v3.1.1/"):
        sub = norm[len("11_post_execution_audit_v3.1.1/") :]
        return not any(sub.startswith(ex) or "/_tmp" in sub for ex in V311_EXCL_PREFIXES)
    return False


def iter_release_files(repo_root: Path, include_manifest_itself: bool = True) -> List[Tuple[str, Path]]:
    """Enumerate sorted (relative_posix, Path) pairs for all release files present on disk."""
    candidates: List[Path] = []
    fixed_files = (
        repo_root / ".gitignore",
        repo_root / "README.md",
        repo_root / "requirements.txt",
        repo_root / "requirements-optional.txt",
        repo_root / "environment.yml",
        repo_root / "02_genome_qc" / "assembly_quality_summary.tsv",
        repo_root / "08_off_target_screening" / "blast_databases" / "database_manifest.tsv",
        repo_root / "10_final_reports" / "BRD_Target_Gene_Sequences.fasta",
        repo_root / "10_final_reports" / "BRD_Target_Gene_Summary.tsv",
        repo_root / "10_final_reports" / "BRD_dsRNA_Candidates.tsv",
    )
    for fp in fixed_files:
        if fp.is_file():
            candidates.append(fp)

    for dname in (
        "brd",
        "config",
        "tests",
        ".licenses",
        "11_post_execution_audit",
        "11_post_execution_audit_v3",
        "11_post_execution_audit_v3.1",
        "11_post_execution_audit_v3.1.1",
    ):
        d = repo_root / dname
        if not d.exists():
            continue
        for fp in d.rglob("*"):
            if fp.is_file() and not _has_runtime_cache_part(fp) and not _has_secret_token(fp):
                candidates.append(fp)

    out: List[Tuple[str, Path]] = []
    seen = set()
    for fp in candidates:
        rel = fp.relative_to(repo_root).as_posix()
        if rel in seen:
            continue
        if not include_manifest_itself and rel == "11_post_execution_audit_v3.1.1/hashes_v3_1_1.tsv":
            continue
        if is_included_release_rel(rel):
            seen.add(rel)
            out.append((rel, fp))
    out.sort(key=lambda x: x[0])
    return out
