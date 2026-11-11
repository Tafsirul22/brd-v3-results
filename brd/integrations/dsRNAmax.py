"""dsRNAmax native Windows executable wrapper."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional
import pandas as pd

from brd.utils.subprocess_runner import run_tool


def run_dsrnamax(
    dsrnamax_bin: str,
    targets_fa: Path,
    out_csv: Path,
    out_log_txt: Path,
    off_targets: Optional[List[Path]] = None,
    construct_len: int = 300,
    kmer_len: int = 21,
    ot_kmer_len: int = 21,
    iterations: int = 100,
    log_path: Optional[Path] = None,
) -> dict:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        dsrnamax_bin,
        "-targets",
        str(targets_fa),
        "-constructLen",
        str(construct_len),
        "-kmerLen",
        str(kmer_len),
        "-otKmerLen",
        str(ot_kmer_len),
        "-iterations",
        str(iterations),
        "-csv",
        str(out_csv),
    ]
    if off_targets:
        cmd.extend(["-offTargets", ",".join(str(p) for p in off_targets)])

    res = run_tool("dsRNAmax", cmd, log_path=log_path, check=False)
    out_log_txt.write_text(res.stdout + "\n--- STDERR ---\n" + res.stderr, encoding="utf-8")

    # Parse generated sense-arm sequence from stdout or csv
    sense_seq = ""
    if out_csv.exists():
        lines = out_csv.read_text(encoding="utf-8", errors="replace").splitlines()
        for i, ln in enumerate(lines):
            if ln.startswith("dsRNA sense-arm sequence:") and i + 1 < len(lines):
                sense_seq = lines[i + 1].strip()
                break

    return {
        "returncode": res.returncode,
        "csv_path": str(out_csv),
        "log_txt": str(out_log_txt),
        "sense_arm_sequence": sense_seq,
        "construct_len": len(sense_seq),
        "stdout": res.stdout,
    }
