"""SIREN integration (siren_prefilter + thermodynamic ViennaRNA duplexfold + siren_designVIII CLI)."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from Bio import SeqIO

try:
    import RNA
except ImportError:
    RNA = None

from brd.utils.subprocess_runner import run_tool


def run_siren_prefilter_cli(
    python_bin: str,
    siren_repo: Path,
    targets_fa: Path,
    gene: str,
    outdir: Path,
    log_path: Optional[Path] = None,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    script = siren_repo / "src" / "siren_rnai" / "siren_prefilter.py"
    run_tool(
        "siren_prefilter",
        [python_bin, str(script), "--targets", str(targets_fa), "--gene", gene, "--outdir", str(outdir), "--write_log"],
        log_path=log_path,
    )
    return outdir / "targets_prefiltered.fa"


def evaluate_sirna_duplexfold_offtargets(
    target_name: str,
    target_seq: str,
    nontarget_fa: Path,
    sirna_size: int = 21,
    step: int = 2,
    mfe_threshold: float = -25.0,
) -> pd.DataFrame:
    """
    Thermodynamic off-target hybridization evaluation for SIREN-generated siRNAs
    against prefiltered sequences using ViennaRNA RNA.duplexfold (native Windows alternative
    to Linux-only RNAhybrid binary; explicitly annotated in tool provenance).
    """
    nt_recs = [r for r in SeqIO.parse(str(nontarget_fa), "fasta") if target_name not in r.id]
    off_target_map: Dict[str, List[str]] = {}
    seq_u = target_seq.upper().replace("T", "U")

    for i in range(0, len(seq_u) - sirna_size + 1, step):
        sirna = seq_u[i : i + sirna_size]
        if "N" in sirna:
            continue
        comp = str.maketrans("ACGU", "UGCA")
        guide_rna = sirna.translate(comp)[::-1]
        sname = f"sirna_{i+1}-{i+sirna_size}"
        for rec in nt_recs:
            t_rna = str(rec.seq).upper().replace("T", "U")
            if len(t_rna) < sirna_size:
                continue
            df = RNA.duplexfold(guide_rna, t_rna)
            if df.energy <= mfe_threshold:
                off_target_map.setdefault(rec.id, []).append(sname)

    rows = []
    for tid, snames in sorted(off_target_map.items(), key=lambda x: -len(x[1])):
        rows.append({"Off target": tid, "siRNA number": len(snames), "siRNA names": ", ".join(snames)})
    return pd.DataFrame(rows, columns=["Off target", "siRNA number", "siRNA names"])


def run_siren_design_cli(
    python_bin: str,
    siren_repo: Path,
    target_name: str,
    target_fa: Path,
    off_targets_summary_tsv: Path,
    out_tsv: Path,
    base_rnai_length: int = 250,
    threads: int = 2,
    log_path: Optional[Path] = None,
) -> pd.DataFrame:
    """Execute SIREN's unmodified siren_designVIII.py CLI script and return scored candidates."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    recs = list(SeqIO.parse(str(target_fa), "fasta"))
    tlen = len(recs[0].seq) if recs else 0
    # Ensure base_rnai_length - 50 <= base_rnai_length + 50 works cleanly on target length
    usable_base = base_rnai_length
    if tlen < base_rnai_length + 50:
        usable_base = 200 if tlen >= 250 else max(150, tlen - 50)
    script = siren_repo / "src" / "siren_rnai" / "siren_designVIII.py"
    run_tool(
        "siren_designVIII",
        [
            python_bin,
            str(script),
            "--target_path",
            str(target_fa),
            "--off_targets_summary_path",
            str(off_targets_summary_tsv),
            "--rnai_seq_length",
            str(usable_base),
            "--threads",
            str(threads),
            "--out",
            str(out_tsv),
        ],
        log_path=log_path,
    )
    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(out_tsv, sep="\t")
    df["target"] = target_name
    return df
