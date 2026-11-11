"""MAFFT and Biopython sequence alignment helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from brd.integrations.blast import _win_short_path
from brd.utils.subprocess_runner import run_tool


def run_mafft(mafft_bin: str, in_fasta: Path, out_fasta: Path, log_path: Optional[Path] = None) -> MultipleSeqAlignment:
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    in_p = _win_short_path(in_fasta)
    res = run_tool("mafft", [mafft_bin, "--auto", "--quiet", in_p], cwd=Path(in_p).parent, log_path=log_path)
    out_fasta.write_text(res.stdout, encoding="utf-8")
    return AlignIO.read(str(out_fasta), "fasta")
