"""QUAST integration wrapper with deterministic Python fallback for Windows."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict

from brd.sequence.fasta_parser import compute_assembly_stats
from brd.utils.subprocess_runner import run_tool


def is_quast_available() -> bool:
    return shutil.which("quast") is not None or shutil.which("quast.py") is not None


def assess_assembly(fasta_path: Path, quast_out_dir: Path, log_path: Path | None = None) -> Dict[str, object]:
    stats = compute_assembly_stats(fasta_path)
    if is_quast_available():
        cmd = [shutil.which("quast") or shutil.which("quast.py"), str(fasta_path), "-o", str(quast_out_dir), "--no-plots"]
        run_tool("quast", cmd, log_path=log_path, check=False)
        stats["qc_engine"] = "QUAST_CLI+Biopython"
    else:
        stats["qc_engine"] = "Biopython_QUAST_compatible (QUAST CLI Linux-only; unavailable natively on Windows)"
    return stats
