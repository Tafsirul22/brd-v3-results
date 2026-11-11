"""Stage 12: Environment Validation & Tool Availability Checker for BRD v3.1."""

from __future__ import annotations

import importlib
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

try:
    from .paths_v3_1 import resolve_repo_root, resolve_tool, resolve_v31_root
except ImportError:
    from paths_v3_1 import resolve_repo_root, resolve_tool, resolve_v31_root


def check_python_package(pkg: str, import_name: str, mandatory: bool) -> Dict[str, Any]:
    try:
        mod = importlib.import_module(import_name)
        ver = getattr(mod, "__version__", "INSTALLED")
        return {"component": pkg, "kind": "python_package", "mandatory": mandatory, "status": "AVAILABLE", "version": str(ver)}
    except Exception as exc:
        return {"component": pkg, "kind": "python_package", "mandatory": mandatory, "status": "MISSING" if mandatory else "OPTIONAL_MISSING", "version": str(exc)}


def check_external_tool(tool_key: str, cmd_args: List[str], mandatory: bool, root: Path, note: str) -> Dict[str, Any]:
    bin_path = resolve_tool(tool_key, root)
    if not bin_path:
        return {"component": tool_key, "kind": "external_binary", "mandatory": mandatory, "status": "MISSING" if mandatory else "UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY", "version": note}
    try:
        p = subprocess.run([bin_path, *cmd_args], capture_output=True, text=True, timeout=10)
        out_lines = (p.stdout or p.stderr or "").strip().splitlines()
        ver = out_lines[0][:120] if out_lines else f"rc={p.returncode}"
        return {"component": tool_key, "kind": "external_binary", "mandatory": mandatory, "status": "AVAILABLE", "version": ver}
    except Exception as exc:
        return {"component": tool_key, "kind": "external_binary", "mandatory": mandatory, "status": "ERROR", "version": str(exc)}


def run_environment_validation(repo_root: Path | None = None) -> Dict[str, Any]:
    root = resolve_repo_root(repo_root)
    v31 = resolve_v31_root(root)
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # Verify temporary directory write & cleanup safety
    tmp_ok = False
    with tempfile.TemporaryDirectory() as td:
        tp = Path(td) / "utf8_crlf_check.txt"
        tp.write_text("line1\nline2\n", encoding="utf-8", newline="\n")
        tmp_ok = tp.read_bytes() == b"line1\nline2\n"

    rows = [
        check_python_package("biopython", "Bio", True),
        check_python_package("pandas", "pandas", True),
        check_python_package("numpy", "numpy", True),
        check_python_package("scipy", "scipy", True),
        check_python_package("matplotlib", "matplotlib", True),
        check_python_package("networkx", "networkx", True),
        check_python_package("PyYAML", "yaml", True),
        check_python_package("tabulate", "tabulate", True),
        check_python_package("reportlab", "reportlab", True),
        check_python_package("pytest", "pytest", True),
        check_python_package("ViennaRNA", "RNA", False),
        check_python_package("primer3-py", "primer3", False),
        check_external_tool("broccoli", ["-help"], False, root, "Broccoli v2.0.1 compiled binary"),
        check_external_tool("diamond", ["version"], False, root, "DIAMOND aligner"),
        check_external_tool("blastn", ["-version"], False, root, "NCBI BLAST+"),
        check_external_tool("makeblastdb", ["-version"], False, root, "NCBI makeblastdb"),
        check_external_tool("datasets", ["--version"], False, root, "NCBI Datasets CLI"),
        check_external_tool("miniprot", ["--version"], False, root, "Linux x64 binary required by compleasm"),
        check_external_tool("hmmsearch", ["-h"], False, root, "HMMER 3.4 POSIX binary required by compleasm"),
        check_external_tool("FastTree", ["-help"], False, root, "Optional (Broccoli uses built-in -phylogenies nj)"),
        check_external_tool("docker", ["--version"], False, root, "Docker Desktop unavailable on host"),
        check_external_tool("quast", ["--version"], False, root, "Linux/WSL2 only"),
        check_external_tool("busco", ["--version"], False, root, "Linux/WSL2 only"),
        check_external_tool("orthofinder", ["-h"], False, root, "Linux/WSL2 only"),
        check_external_tool("RNAhybrid", ["-h"], False, root, "Linux/macOS Bioconda only"),
    ]
    df = pd.DataFrame(rows)
    tsv_path = v31 / "tool_availability_v3_1.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    mandatory_missing = df[(df["mandatory"]) & (df["status"] != "AVAILABLE")]
    md = f"""# BRD v3.1 Environment & Dependency Validation Report

**Checked:** {now_utc}  
**OS / Platform:** `{platform.platform()}` (`{os.name}`)  
**Python:** `{platform.python_version()}` (`{sys.executable}`)  
**Default Encoding / Temp UTF-8 LF Safety:** `{sys.getdefaultencoding()}` / `tmp_utf8_lf_ok={tmp_ok}`  
**Mandatory Python Packages Missing:** `{len(mandatory_missing)}`

## Complete Dependency Matrix (`tool_availability_v3_1.tsv`)
{df.to_markdown(index=False)}
"""
    md_path = v31 / "V3_1_Environment_Report.md"
    md_path.write_text(md, encoding="utf-8")
    return {"tsv": tsv_path, "report": md_path, "mandatory_ok": len(mandatory_missing) == 0}


if __name__ == "__main__":
    out = run_environment_validation()
    print(f"ENV_CHECK: {out}")
