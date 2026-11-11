"""DIAMOND protein alignment integration."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd

from brd.utils.subprocess_runner import run_tool

DIAMOND_COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovhsp",
    "scovhsp",
]


def build_diamond_db(diamond_bin: str, fasta_path: Path, db_prefix: Path, log_path: Optional[Path] = None, force: bool = False) -> Path:
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    dmnd = Path(str(db_prefix) + ".dmnd")
    needs_build = force or (not dmnd.exists()) or (fasta_path.exists() and dmnd.stat().st_mtime < fasta_path.stat().st_mtime)
    if needs_build:
        run_tool(
            "diamond",
            [diamond_bin, "makedb", "--in", str(fasta_path), "-d", str(db_prefix)],
            log_path=log_path,
        )
    return dmnd


def run_diamond_blastp(
    diamond_bin: str,
    query_faa: Path,
    db_prefix: Path,
    out_tsv: Path,
    evalue: float = 1e-5,
    max_target_seqs: int = 25,
    threads: int = 4,
    log_path: Optional[Path] = None,
) -> pd.DataFrame:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        diamond_bin,
        "blastp",
        "-q",
        str(query_faa),
        "-d",
        str(db_prefix),
        "-o",
        str(out_tsv),
        "--evalue",
        str(evalue),
        "--max-target-seqs",
        str(max_target_seqs),
        "--threads",
        str(threads),
        "--sensitive",
        "--outfmt",
        "6",
        *DIAMOND_COLS,
    ]
    run_tool("diamond", cmd, log_path=log_path)
    if out_tsv.stat().st_size == 0:
        return pd.DataFrame(columns=DIAMOND_COLS)
    return pd.read_csv(out_tsv, sep="\t", names=DIAMOND_COLS)
