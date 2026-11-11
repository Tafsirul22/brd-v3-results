"""CDS completeness, translation verification, and sequence QC."""

from __future__ import annotations

from typing import Dict
from Bio.Seq import Seq


def validate_cds_translation(cds_seq: str, expected_protein: str | None = None) -> Dict[str, object]:
    """Validate CDS length mod 3, start/stop codons, internal stops, and translation agreement."""
    s = cds_seq.upper().strip()
    length = len(s)
    mod3_ok = (length % 3 == 0) and length >= 3
    has_start = s.startswith("ATG")
    has_stop = s[-3:] in ("TAA", "TAG", "TGA") if length >= 3 else False
    has_ambig = any(c not in "ACGT" for c in s)

    internal_stops = 0
    translated = ""
    if length >= 3 and not has_ambig:
        trim_len = length - (length % 3)
        translated = str(Seq(s[:trim_len]).translate(table=1))
        body = translated[:-1] if translated.endswith("*") else translated
        internal_stops = body.count("*")
    prot_clean = translated.rstrip("*") if translated else ""
    matches_expected = None
    if expected_protein is not None and prot_clean:
        matches_expected = (prot_clean == expected_protein.strip().rstrip("*"))

    return {
        "cds_length_nt": length,
        "mod3_ok": mod3_ok,
        "has_start_atg": has_start,
        "has_terminal_stop": has_stop,
        "has_ambiguous_bases": has_ambig,
        "internal_stop_count": internal_stops,
        "translated_aa_len": len(prot_clean),
        "matches_reference_protein": matches_expected,
        "is_complete_cds": bool(mod3_ok and has_start and has_stop and internal_stops == 0 and not has_ambig),
    }
