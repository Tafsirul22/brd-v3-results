"""BUSCO integration wrapper (Linux-only external tool check, never fabricates scores)."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict


def is_busco_available() -> bool:
    return shutil.which("busco") is not None


def get_busco_status(assembly_acc: str, genome_fasta: Path, lineage: str = "fungi_odb10") -> Dict[str, object]:
    if is_busco_available():
        return {
            "assembly_accession": assembly_acc,
            "lineage": lineage,
            "busco_version": "CLI_detected",
            "status": "AVAILABLE",
            "command": f"busco -i {genome_fasta} -l {lineage} -m genome -o {assembly_acc}_busco",
        }
    return {
        "assembly_accession": assembly_acc,
        "lineage": lineage,
        "busco_version": "Not installed (v6.1.0 Linux/macOS Bioconda only)",
        "complete_pct": "NA",
        "single_copy_pct": "NA",
        "duplicated_pct": "NA",
        "fragmented_pct": "NA",
        "missing_pct": "NA",
        "status": "NOT_RUN_LINUX_ONLY_NO_WSL_DISTRO",
        "command": f"busco -i '{genome_fasta}' -l {lineage} -m genome -o '{assembly_acc}_busco'",
        "notes": "Official BUSCO requires Linux/macOS (MetaEuk/HMMER/Augustus). WSL2 has no Linux distro installed; scores left NA (not fabricated).",
    }
