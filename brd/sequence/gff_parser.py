"""GFF3 parsing, coordinate validation, and harmonisation utilities."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


def parse_attributes(attr_str: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in attr_str.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def validate_gff3(gff_path: Path, seq_lengths: Dict[str, int]) -> Dict[str, object]:
    """Validate GFF3 coordinates, strands, phases, and parent-child structure."""
    genes: Dict[str, dict] = {}
    mrnas: Dict[str, dict] = {}
    cds_by_parent: Dict[str, List[dict]] = defaultdict(list)
    errors: List[str] = []
    feature_counts: Dict[str, int] = defaultdict(int)

    with open(gff_path, "r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                errors.append(f"Line {lineno}: expected 9 columns, got {len(parts)}")
                continue
            seqid, source, ftype, start_s, end_s, score, strand, phase, attrs_s = parts
            feature_counts[ftype] += 1
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                errors.append(f"Line {lineno}: non-integer coordinates {start_s}..{end_s}")
                continue
            if start < 1 or end < start:
                errors.append(f"Line {lineno}: invalid coordinate interval {start}..{end}")
            if seqid in seq_lengths and end > seq_lengths[seqid]:
                errors.append(f"Line {lineno}: end {end} > scaffold {seqid} len {seq_lengths[seqid]}")
            if strand not in ("+", "-", "."):
                errors.append(f"Line {lineno}: invalid strand {strand}")
            attrs = parse_attributes(attrs_s)
            fid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")

            if ftype == "gene":
                genes[fid] = {
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "attrs": attrs,
                }
            elif ftype in ("mRNA", "transcript"):
                mrnas[fid] = {
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "parent": parent,
                    "attrs": attrs,
                }
            elif ftype == "CDS":
                if phase not in ("0", "1", "2"):
                    errors.append(f"Line {lineno}: CDS missing valid phase ({phase})")
                for p in parent.split(","):
                    cds_by_parent[p].append(
                        {"seqid": seqid, "start": start, "end": end, "strand": strand, "phase": int(phase) if phase in ("0", "1", "2") else 0, "attrs": attrs}
                    )

    orphan_mrnas = sum(1 for m in mrnas.values() if m["parent"] and m["parent"] not in genes)
    return {
        "genes": genes,
        "mrnas": mrnas,
        "cds_by_parent": dict(cds_by_parent),
        "feature_counts": dict(feature_counts),
        "orphan_mrnas": orphan_mrnas,
        "error_count": len(errors),
        "sample_errors": errors[:20],
    }
