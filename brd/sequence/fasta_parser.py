"""FASTA parsing and assembly metric calculation."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List
from Bio import SeqIO


def compute_assembly_stats(fasta_path: Path) -> Dict[str, object]:
    """Compute QUAST-compatible assembly statistics deterministically."""
    lengths: List[int] = []
    gc_count = 0
    n_count = 0
    ambig_count = 0
    headers: List[str] = []
    seq_hashes: List[str] = []
    invalid_chars = set()

    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(rec.seq).upper()
        l = len(seq)
        lengths.append(l)
        headers.append(rec.id)
        seq_hashes.append(hashlib.md5(seq.encode("ascii", errors="ignore")).hexdigest())
        gc_count += seq.count("G") + seq.count("C")
        n_curr = seq.count("N")
        n_count += n_curr
        atgcn = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C") + n_curr
        ambig_count += (l - atgcn)
        bad = set(seq) - set("ACGTNRYSWKMBDHV")
        if bad:
            invalid_chars.update(bad)

    if not lengths:
        raise ValueError(f"Empty FASTA: {fasta_path}")

    lengths_sorted = sorted(lengths, reverse=True)
    total_len = sum(lengths_sorted)
    ungapped_len = total_len - n_count

    def nx_lx(frac: float):
        target = total_len * frac
        acc = 0
        for idx, l in enumerate(lengths_sorted, 1):
            acc += l
            if acc >= target:
                return l, idx
        return 0, 0

    n50, l50 = nx_lx(0.50)
    n90, l90 = nx_lx(0.90)
    dup_headers = len(headers) - len(set(headers))
    dup_seqs = len(seq_hashes) - len(set(seq_hashes))
    gc_pct = (gc_count / ungapped_len * 100.0) if ungapped_len > 0 else 0.0
    n_pct = (n_count / total_len * 100.0) if total_len > 0 else 0.0

    return {
        "scaffold_count": len(lengths_sorted),
        "total_length_bp": total_len,
        "ungapped_length_bp": ungapped_len,
        "largest_scaffold_bp": lengths_sorted[0],
        "smallest_scaffold_bp": lengths_sorted[-1],
        "n50_bp": n50,
        "l50": l50,
        "n90_bp": n90,
        "l90": l90,
        "gc_percent": round(gc_pct, 3),
        "n_count": n_count,
        "n_percent": round(n_pct, 4),
        "non_atgcn_ambiguous": ambig_count,
        "duplicate_headers": dup_headers,
        "duplicate_sequences": dup_seqs,
        "invalid_chars": "".join(sorted(invalid_chars)),
    }
