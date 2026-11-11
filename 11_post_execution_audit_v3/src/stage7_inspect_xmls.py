"""Inspect exact construct details, primer contexts, and verified sequences in Stage 7 XMLs."""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from paths_v3 import resolve_v3_root


def main() -> None:
    import urllib.request
    ft_dir = resolve_v3_root() / "stage7_literature_benchmark" / "fulltexts"
    for extra in ("PMC6859489", "PMC4515635", "PMC6764512"):
        xp = ft_dir / f"{extra}.xml"
        if not xp.exists():
            try:
                xb = urllib.request.urlopen(urllib.request.Request(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{extra}/fullTextXML", headers={"User-Agent": "BRD-v3/1.0"}), timeout=15).read()
                if b"<article" in xb:
                    xp.write_bytes(xb)
                    print(f"FETCHED {extra}.xml ({len(xb)} B)")
            except Exception as e:
                print(f"SKIP {extra}: {e}")
    for xp in sorted(ft_dir.glob("*.xml")):
        pid = xp.stem
        root = ET.parse(xp).getroot()
        txt = " ".join("".join(root.itertext()).split())
        sn = next(("".join(e.itertext()).strip() for e in root.findall(".//contrib/name/surname")), "-")
        yr = next(("".join(e.itertext()).strip() for e in root.findall(".//pub-date/year")), "-")
        doi = ""
        for a in root.findall(".//article-id"):
            if a.attrib.get("pub-id-type") == "doi":
                doi = "".join(a.itertext()).strip()
        nt_hits = [m.group(0) for m in re.finditer(r"\b[ACGTU]{19,600}\b", txt) if len(set(m.group(0))) >= 3]
        if nt_hits or pid in ("PMC6859489", "PMC13608320", "PMC13205770"):
            print(f"=== {pid} | {sn} {yr} | {doi} | nts={len(nt_hits)} ===")
            for s in nt_hits[:4]:
                pos = txt.find(s)
                print(f"  NT({len(s)}): {s}  <-- {txt[max(0, pos - 40) : min(len(txt), pos + len(s) + 30)]}")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
