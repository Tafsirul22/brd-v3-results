"""NCBI BLAST+ integration (makeblastdb, blastn, blastp, tblastn, blastx)."""

from __future__ import annotations

import ctypes
import os
from pathlib import Path
from typing import List, Optional

import pandas as pd

from brd.utils.subprocess_runner import run_tool

BLAST_COLS = [
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
    "qlen",
    "slen",
    "qcovs",
    "sstrand",
]


def _win_short_path(path: Path) -> str:
    """Convert Windows paths with spaces to 8.3 short paths when supported by NTFS."""
    s = str(path)
    if os.name != "nt":
        return s
    parent = path.parent
    parent.mkdir(parents=True, exist_ok=True)
    buf = ctypes.create_unicode_buffer(2048)
    # GetShortPathNameW requires target path (or parent for output prefixes) to exist
    target = str(path if path.exists() else parent)
    res = ctypes.windll.kernel32.GetShortPathNameW(target, buf, 2048)
    if res > 0 and buf.value:
        if path.exists():
            return buf.value
        return str(Path(buf.value) / path.name)
    return s


def make_blast_db(
    makeblastdb_bin: str,
    fasta_path: Path,
    db_prefix: Path,
    dbtype: str = "nucl",
    log_path: Optional[Path] = None,
) -> Path:
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    ext = ".nsq" if dbtype == "nucl" else ".psq"
    nal = Path(str(db_prefix) + ".nal")
    sq = Path(str(db_prefix) + ext)
    ndb = Path(str(db_prefix) + ".ndb")
    if not sq.exists() and not nal.exists() and not ndb.exists():
        in_p = _win_short_path(fasta_path)
        out_p = _win_short_path(db_prefix)
        run_tool(
            "makeblastdb",
            [makeblastdb_bin, "-in", in_p, "-dbtype", dbtype, "-out", out_p],
            cwd=Path(out_p).parent,
            log_path=log_path,
        )
    return db_prefix


def run_blast(
    blast_bin: str,
    query_fasta: Path,
    db_prefix: Path,
    out_tsv: Path,
    evalue: float = 1e-3,
    task: Optional[str] = None,
    word_size: Optional[int] = None,
    max_target_seqs: int = 50,
    num_threads: int = 4,
    extra_args: Optional[List[str]] = None,
    log_path: Optional[Path] = None,
) -> pd.DataFrame:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fmt_str = "6 " + " ".join(BLAST_COLS)
    q_p = _win_short_path(query_fasta)
    db_p = _win_short_path(db_prefix)
    o_p = _win_short_path(out_tsv)
    cmd = [
        blast_bin,
        "-query",
        q_p,
        "-db",
        db_p,
        "-out",
        o_p,
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(num_threads),
        "-outfmt",
        fmt_str,
    ]
    if task:
        cmd.extend(["-task", task])
    if word_size:
        cmd.extend(["-word_size", str(word_size)])
    if extra_args:
        cmd.extend(extra_args)
    run_tool(Path(blast_bin).stem, cmd, cwd=Path(o_p).parent, log_path=log_path)
    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLS)
    return pd.read_csv(out_tsv, sep="\t", names=BLAST_COLS)
