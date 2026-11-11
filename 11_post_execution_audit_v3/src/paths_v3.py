"""Portable path and configuration resolver for BRD v3 (no hardcoded machine paths)."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Dict

import yaml


def resolve_output_root(start: Path | None = None) -> Path:
    """Resolve the BRD Output workspace root dynamically or via BRD_OUTPUT_DIR."""
    env_root = os.environ.get("BRD_OUTPUT_DIR")
    if env_root:
        p = Path(env_root).expanduser().resolve()
        if p.exists():
            return p
    probe = (start or Path(__file__)).resolve()
    for parent in [probe] + list(probe.parents):
        if (parent / "config" / "project_configuration.yaml").exists() and (parent / "brd").is_dir():
            return parent
    raise FileNotFoundError("Could not locate BRD Output root containing config/project_configuration.yaml")


def resolve_v3_root(output_root: Path | None = None) -> Path:
    out = output_root or resolve_output_root()
    v3 = out / "11_post_execution_audit_v3"
    v3.mkdir(parents=True, exist_ok=True)
    return v3


def load_yaml_portable(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def resolve_tool_binary(tool_name: str, cfg_tools: Dict[str, str], output_root: Path) -> str | None:
    """Resolve external tool binary across relative workspace tool paths or system PATH."""
    rel_candidates = {
        "diamond": output_root / "tools" / "bin" / ("diamond.exe" if os.name == "nt" else "diamond"),
        "datasets": output_root / "tools" / "bin" / ("datasets.exe" if os.name == "nt" else "datasets"),
        "dsrnamax": output_root / "tools" / "bin" / ("dsRNAmax.exe" if os.name == "nt" else "dsRNAmax"),
        "mafft": output_root / "tools" / "mafft" / "mafft-win" / "mafft.bat",
    }
    if tool_name in rel_candidates and rel_candidates[tool_name].exists():
        return str(rel_candidates[tool_name])
    cfg_val = cfg_tools.get(tool_name, "")
    if cfg_val and Path(cfg_val).exists():
        return str(Path(cfg_val))
    found = shutil.which(tool_name)
    return found
