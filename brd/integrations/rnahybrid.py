"""RNAhybrid availability wrapper and ViennaRNA duplexfold thermodynamic helper."""

from __future__ import annotations

import shutil
from typing import Dict

try:
    import RNA  # Optional thermodynamic dependency (ViennaRNA)
except ImportError:
    RNA = None


def is_rnahybrid_available() -> bool:
    return shutil.which("RNAhybrid") is not None


def is_viennarna_available() -> bool:
    return RNA is not None


def compute_duplex_mfe(guide_rna: str, target_rna: str) -> Dict[str, object]:
    """Compute thermodynamic RNA-RNA duplex MFE via ViennaRNA RNA.duplexfold."""
    if RNA is None:
        raise RuntimeError("Optional dependency ViennaRNA ('RNA' module) is not installed in this environment.")
    df = RNA.duplexfold(guide_rna.upper().replace("T", "U"), target_rna.upper().replace("T", "U"))
    return {
        "mfe_kcal_mol": round(float(df.energy), 2),
        "structure": df.structure,
        "engine": "ViennaRNA_duplexfold_2.8.0 (RNAhybrid CLI Linux-only)",
    }
