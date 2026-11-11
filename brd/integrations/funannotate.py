"""Funannotate and BRAKER3 availability wrappers (Linux-only genome annotation pipelines)."""

from __future__ import annotations

import shutil
from typing import Dict


def check_funannotate() -> Dict[str, str]:
    exe = shutil.which("funannotate")
    return {
        "tool": "Funannotate",
        "repository": "https://github.com/nextgenusfs/funannotate",
        "available": str(exe is not None),
        "status": "AVAILABLE" if exe else "NOT_RUN_LINUX_ONLY (use existing NCBI RefSeq/GenBank annotations; WSL2/Conda required for de novo annotation)",
    }


def check_braker() -> Dict[str, str]:
    exe = shutil.which("braker.pl")
    return {
        "tool": "BRAKER3",
        "repository": "https://github.com/Gaius-Augustus/BRAKER",
        "available": str(exe is not None),
        "status": "AVAILABLE" if exe else "NOT_RUN_LINUX_ONLY (requires Linux + GeneMark-ES/ETP license + AUGUSTUS)",
    }
