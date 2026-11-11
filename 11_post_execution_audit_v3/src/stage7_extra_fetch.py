"""Fetch Leng 2015, Chen 2022, Wang 2026, LUCID, Hao 2021, Mukherjee 2024 XMLs & inspect dsRNA constructs."""

from __future__ import annotations

import json
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from paths_v3 import resolve_v3_root


def main() -> None:
    v3 = resolve_v3_root()
    ft = v3 / "stage7_literature_benchmark" / "fulltexts"
    raw = v3 / "stage7_literature_benchmark" / "raw_queries"
    queries = [
        'DOI:"10.1371/journal.pone.0128291"',
        'DOI:"10.3390/ijms232112922"',
        'DOI:"10.3390/jof12090656"',
        'DOI:"10.1111/pbi.70526"',
        'SsAGO2 AND Sclerotinia AND dsRNA',
        'TRI6 AND "Fusarium graminearum" AND dsRNA AND Hao',
    ]
    found = {}
    for q in queries:
        try:
            url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
                {"query": q, "resultType": "lite", "pageSize": 3, "format": "json"}
            )
            data = json.loads(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "BRD-v3/1.0"}), timeout=20).read().decode("utf-8"))
            res = data.get("resultList", {}).get("result", [])
        except Exception as exc:
            print("Q_ERR", q, exc, flush=True)
            res = []
        found[q] = res[:2]
        for r in res[:1]:
            pmcid = str(r.get("pmcid", "") or "")
            print(f"{q[:35]:35s} -> {pmcid or '-':11s} {r.get('doi','-')} | {(r.get('title','') or '')[:60]}", flush=True)
            if pmcid and not (ft / f"{pmcid}.xml").exists():
                for u in (
                    f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
                    f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid}&retmode=xml",
                ):
                    try:
                        xb = urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "BRD-v3/1.0"}), timeout=20).read()
                        if b"<article" in xb and len(xb) > 2000:
                            (ft / f"{pmcid}.xml").write_bytes(xb)
                            print(f"   DL_XML {pmcid} ({len(xb)} B)", flush=True)
                            break
                    except Exception:
                        continue
    (raw / "epmc_cited_resolutions.json").write_text(json.dumps(found, indent=2), encoding="utf-8")

    # Also extract exact Methods sentences about dsRNA construct sequence / size / primers in PMC6465682, PMC5063301, PMC7202221, PMC8710829, PMC7072719
    for pmcid in ("PMC6465682", "PMC5063301", "PMC7202221", "PMC8710829", "PMC7072719", "PMC5943259"):
        xp = ft / f"{pmcid}.xml"
        if not xp.exists():
            continue
        root = ET.parse(xp).getroot()
        txt = " ".join("".join(root.itertext()).split())
        print(f"--- CONSTRUCT DETAILS IN {pmcid} ---", flush=True)
        for m in re.finditer(r"[^.?!]{0,120}(?:primer|amplif|CYP3RNA|CYP3-dsRNA|FgAGO|FgDCL|BcDCL|ergosterol|asiR\d+|376 bp|791|SS1G_)[^.?!]{0,160}[.?!]", txt, re.I):
            print("  *", m.group(0).strip()[:220], flush=True)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
