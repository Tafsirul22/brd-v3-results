"""Portable path and tool resolver for BRD v3.1 (relative/env-resolved, zero hardcoded machine paths)."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Dict


def resolve_repo_root(start: Path | None = None) -> Path:
    env_root = os.environ.get("BRD_OUTPUT_DIR")
    if env_root:
        p = Path(env_root).expanduser().resolve()
        if p.exists():
            return p
    probe = (start or Path(__file__)).resolve()
    for parent in [probe] + list(probe.parents):
        if (parent / "config" / "project_configuration.yaml").exists() and (parent / "brd").is_dir():
            return parent
    raise FileNotFoundError("Could not locate BRD workspace root containing brd/ and config/project_configuration.yaml")


def resolve_v31_root(repo_root: Path | None = None) -> Path:
    root = repo_root or resolve_repo_root()
    v31 = root / "11_post_execution_audit_v3.1"
    v31.mkdir(parents=True, exist_ok=True)
    return v31


def resolve_tool(tool_name: str, repo_root: Path | None = None) -> str | None:
    root = repo_root or resolve_repo_root()
    ext = ".exe" if os.name == "nt" else ""
    local_map: Dict[str, Path] = {
        "broccoli": root / "tools" / "bin" / f"broccoli{ext}",
        "diamond": root / "tools" / "bin" / f"diamond{ext}",
        "datasets": root / "tools" / "bin" / f"datasets{ext}",
        "dsrnamax": root / "tools" / "bin" / f"dsRNAmax{ext}",
        "miniprot": root / "tools" / "bin" / f"miniprot{ext}",
        "hmmsearch": root / "tools" / "bin" / f"hmmsearch{ext}",
        "fasttree": root / "tools" / "bin" / f"FastTree{ext}",
    }
    if tool_name.lower() in local_map and local_map[tool_name.lower()].exists():
        return str(local_map[tool_name.lower()])
    return shutil.which(tool_name)
