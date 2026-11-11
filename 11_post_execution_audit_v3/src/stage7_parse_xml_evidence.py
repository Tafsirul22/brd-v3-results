"""Parse exact tables, construct lengths, primer pairs, accessions, and methods from Stage 7 JATS XMLs."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

from paths_v3 import resolve_output_root, resolve_v3_root


FOCUS_PMCIDS = [
    "PMC5943259", "PMC7202221", "PMC6859480", "PMC5063301", "PMC8428832",
    "PMC5040644", "PMC7072719", "PMC8710829", "PMC7680546", "PMC7214474",
    "PMC7391711", "PMC6465682", "PMC8743017", "PMC3391263", "PMC4793245",
]


def _txt(el: ET.Element | None) -> str:
    return " ".join("".join(el.itertext()).split()) if el is not None else ""


def main() -> None:
    v3 = resolve_v3_root()
    ft_dir = v3 / "stage7_literature_benchmark" / "fulltexts"
    out_txt = v3 / "stage7_literature_benchmark" / "xml_extracted_details.txt"
    lines: list[str] = []

    for pmcid in FOCUS_PMCIDS:
        xp = ft_dir / f"{pmcid}.xml"
        if not xp.exists():
            continue
        root = ET.parse(xp).getroot()
        title = _txt(root.find(".//article-title"))
        doi = ""
        for aid in root.findall(".//article-id"):
            if aid.attrib.get("pub-id-type") == "doi":
                doi = _txt(aid)
                break
        authors = []
        for contrib in root.findall(".//contrib[@contrib-type='author']")[:2]:
            sn = _txt(contrib.find(".//surname"))
            if sn:
                authors.append(sn)
        year = _txt(root.find(".//pub-date/year"))
        full_txt = _txt(root)

        lines.append(f"\n==================== {pmcid} ({', '.join(authors)} et al., {year}; DOI: {doi}) ====================")
        lines.append(f"TITLE: {title}")
        abs_txt = _txt(root.find(".//abstract"))
        lines.append(f"ABSTRACT: {abs_txt[:900]}")

        # Search sentences mentioning dsRNA length (bp/nt), silencing %, accession, or uptake
        sents = re.split(r"(?<=[.!?])\s+", full_txt)
        hits = [
            s for s in sents
            if re.search(r"\b(\d{2,4}\s?(?:bp|nt)|dsRNA|SIGS|HIGS|silenc|lesion|inhibit|uptake|FGSG_|SS1G_|BCIN_|CYP51|DCL|AGO|CHI|MAPK|PKS)\b", s, re.I)
        ]
        for h in hits[:14]:
            if any(k in h.lower() for k in ("dsrna", "bp", "nt ", "silenc", "inhibit", "lesion", "uptake", "accession", "fgsg", "ss1g", "cyp51", "dcl")):
                lines.append(f"  SENT: {h[:320]}")

        # Extract tables completely (bounded to 35 rows/table)
        for ti, tw in enumerate(root.findall(".//table-wrap"), 1):
            cap = _txt(tw.find(".//caption"))
            lines.append(f"  -- TABLE {ti}: {cap[:220]}")
            for tr in tw.findall(".//tr")[:35]:
                cells = [_txt(td) for td in list(tr) if td.tag in ("th", "td")]
                if cells:
                    lines.append("    | " + " | ".join(cells) + " |")

    out_txt.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out_txt} ({len(lines)} lines)")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
