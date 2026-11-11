"""FASTA, GFF3, and file validation helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO

VALID_DNA = set("ACGTNacgtn")
VALID_STRICT_DNA = set("ACGTacgt")


def validate_fasta(path: Path, sequence_type: str = "nucl") -> Dict[str, object]:
    """Validate FASTA readability, non-emptiness, duplicate IDs, and alphabet."""
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"FASTA file missing or empty: {path}")
    ids: List[str] = []
    seen = set()
    dup_ids: List[str] = []
    invalid_chars = set()
    total_len = 0
    n_count = 0
    valid_set = VALID_DNA if sequence_type == "nucl" else set("ACDEFGHIKLMNPQRSTVWYXBZJUO*acdefghiklmnpqrstvwyxbzjuo")
    for rec in SeqIO.parse(str(path), "fasta"):
        if rec.id in seen:
            dup_ids.append(rec.id)
        seen.add(rec.id)
        ids.append(rec.id)
        s = str(rec.seq)
        total_len += len(s)
        if sequence_type == "nucl":
            n_count += s.upper().count("N")
        bad = set(s) - valid_set
        if bad:
            invalid_chars.update(bad)
    if not ids:
        raise ValueError(f"No FASTA records parsed in {path}")
    return {
        "path": str(path),
        "records": len(ids),
        "total_length": total_len,
        "duplicate_ids": len(dup_ids),
        "n_count": n_count,
        "invalid_chars": sorted(invalid_chars),
        "valid": len(dup_ids) == 0 and len(invalid_chars) == 0,
    }


def validate_nucleotide_seq(seq: str) -> Tuple[bool, Dict[str, object]]:
    """Validate an individual nucleotide sequence for alphabet and ambiguous bases."""
    s = seq.upper().strip()
    bad = sorted(set(s) - set("ACGTN"))
    n_bases = s.count("N")
    gc = (s.count("G") + s.count("C")) / len(s) * 100.0 if s else 0.0
    return (len(s) > 0 and len(bad) == 0), {
        "length": len(s),
        "n_bases": n_bases,
        "gc_percent": round(gc, 2),
        "invalid_chars": bad,
    }
