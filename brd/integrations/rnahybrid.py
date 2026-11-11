"""RNAhybrid availability wrapper and ViennaRNA duplexfold thermodynamic helper."""

from __future__ import annotations

import shutil
from typing import Dict
import RNA


def is_rnahybrid_available() -> bool:
    return shutil.which("RNAhybrid") is not None


def compute_duplex_mfe(guide_rna: str, target_rna: str) -> Dict[str, object]:
    """Compute thermodynamic RNA-RNA duplex MFE via ViennaRNA RNA.duplexfold."""
    df = RNA.duplexfold(guide_rna.upper().replace("T", "U"), target_rna.upper().replace("T", "U"))
    return {
        "mfe_kcal_mol": round(float(df.energy), 2),
        "structure": df.structure,
        "engine": "ViennaRNA_duplexfold_2.8.0 (RNAhybrid CLI Linux-only)",
    }
