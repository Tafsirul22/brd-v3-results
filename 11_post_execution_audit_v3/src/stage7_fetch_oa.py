"""Stage 7 Helper: Cleanly query Europe PMC OA studies & fetch JATS full-text XMLs."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

from paths_v3 import resolve_output_root, resolve_v3_root


EPMC_DEFAULT = Path.home() / ".gemini" / "config" / "plugins" / "science" / "skills" / "literature_search_europepmc" / "scripts" / "europepmc_api.py"


def _find_uv() -> str | None:
    found = shutil.which("uv")
    if found:
        return found
    local_uv = Path.home() / ".local" / "bin" / "uv.exe"
    return str(local_uv) if local_uv.exists() else None


def _epmc_search_direct(query: str, dest: Path, max_results: int = 15) -> int:
    if dest.exists() and dest.stat().st_size > 200:
        return len(json.loads(dest.read_text(encoding="utf-8")).get("results", []))
    params = urllib.parse.urlencode({"query": query, "resultType": "lite", "pageSize": max_results, "sort": "CITED desc", "format": "json"})
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{params}"
    for attempt in range(2):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "BRD-v3-Research/1.0"})
            with urllib.request.urlopen(req, timeout=25) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            results = data.get("resultList", {}).get("result", [])
            dest.write_text(json.dumps({"hitCount": data.get("hitCount", len(results)), "results": results}, indent=2), encoding="utf-8")
            return len(results)
        except Exception as exc:
            print(f"WARN query retry {attempt+1} for {dest.name}: {exc}", flush=True)
    return 0


def _fetch_xml_direct(pmcid: str, dest: Path) -> str:
    if dest.exists() and dest.stat().st_size > 2000:
        return f"CACHED ({dest.stat().st_size} B)"
    urls = [
        ("EPMC_REST", f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"),
        ("NCBI_EFETCH", f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid}&retmode=xml"),
    ]
    for src_tag, url in urls:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "BRD-v3-Research/1.0"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = resp.read()
            if (b"<article" in data or b"<pmc-articleset" in data) and len(data) > 2000:
                dest.write_bytes(data)
                return f"{src_tag}_OK ({len(data)} B)"
        except Exception:
            continue
    return "UNAVAILABLE_OR_CLOSED"


def main() -> None:
    out = resolve_output_root()
    v3 = resolve_v3_root(out)
    lit_dir = v3 / "stage7_literature_benchmark"
    raw_dir = lit_dir / "raw_queries"
    ft_dir = lit_dir / "fulltexts"
    raw_dir.mkdir(parents=True, exist_ok=True)
    ft_dir.mkdir(parents=True, exist_ok=True)

    # Targeted queries specifically for fungal dsRNA / SIGS / target KO primary papers
    queries = {
        "epmc_alternaria_sigs_specific.json": (
            'TITLE:"Alternaria alternata" AND (RNAi OR dsRNA OR SIGS OR silencing) AND OPEN_ACCESS:y'
        ),
        "epmc_aspergillus_rnai_specific.json": (
            'TITLE:"Aspergillus flavus" AND (RNAi OR dsRNA OR SIGS OR silencing) AND OPEN_ACCESS:y'
        ),
        "epmc_verticillium_magnaporthe_sigs.json": (
            '(TITLE:"Verticillium dahliae" OR TITLE:"Magnaporthe oryzae") AND (dsRNA OR SIGS OR "spray-induced") AND OPEN_ACCESS:y'
        ),
        "epmc_sclerotinia_sigs_specific.json": (
            'TITLE:"Sclerotinia sclerotiorum" AND (dsRNA OR SIGS OR "spray-induced") AND OPEN_ACCESS:y'
        ),
        "epmc_bsor_specific_papers.json": (
            '("Bipolaris sorokiniana" OR "Cochliobolus sativus") AND (CsFUS3 OR CsSTE11 OR CsSp1 OR CsPKS1 OR "spot blotch" dsRNA)'
        ),
    }
    for fname, q in queries.items():
        n = _epmc_search_direct(q, raw_dir / fname, max_results=15)
        print(f"QUERY {fname}: {n} results", flush=True)

    # High-priority OA PMCIDs directly verified across primary SIGS / B. sorokiniana literature
    pmcids = [
        "PMC8428832",  # Qiao et al. 2021 (PNAS/PBJ: dsRNA uptake across fungal pathogens)
        "PMC5063301",  # Koch et al. 2016 (PLoS Pathog: CYP3-dsRNA SIGS vs F. graminearum)
        "PMC7202221",  # Werner et al. 2020 (Front Plant Sci: AGO/DCL SIGS in F. graminearum)
        "PMC6859480",  # Koch et al. 2019 (F. graminearum SIGS)
        "PMC5943259",  # McLoughlin et al. 2018 (Sci Rep: exogenous dsRNA vs S. sclerotiorum & B. cinerea)
        "PMC8710829",  # Ergosterol dsRNA SIGS (Plant Biotechnol J 2022)
        "PMC7680546",  # Multiplex HIGS in F. graminearum (Plant Biotechnol J 2020)
        "PMC8743017",  # Zhang et al. 2022 (CsSp1 effector in B. sorokiniana!)
        "PMC7072719",  # Biomolecules 2020 dsRNA vs Botrytis cinerea
        "PMC5040644",  # Wang et al. 2016 (Nature Plants, B. cinerea DCL1/DCL2)
        "PMC7214474",  # Hu et al. 2020
    ]

    # Add top targeted hits strictly from species-specific RNAi/SIGS title searches (max 3 per file)
    for qname in ("epmc_alternaria_sigs_specific.json", "epmc_aspergillus_rnai_specific.json", "epmc_verticillium_magnaporthe_sigs.json", "epmc_sclerotinia_sigs_specific.json"):
        qf = raw_dir / qname
        if not qf.exists():
            continue
        res = json.loads(qf.read_text(encoding="utf-8")).get("results", [])
        print(f"=== {qname} ({len(res)}) ===", flush=True)
        added = 0
        for r in res[:6]:
            pmcid = str(r.get("pmcid", "") or "")
            oa = str(r.get("isOpenAccess", "") or "")
            title = (r.get("title", "") or "")[:80]
            print(f"  {pmcid or '-':11s} OA={oa:1s} yr={r.get('pubYear','-')} cit={r.get('citedByCount',0):4} doi={r.get('doi','-')} | {title}", flush=True)
            if pmcid and oa == "Y" and pmcid not in pmcids and added < 3:
                pmcids.append(pmcid)
                added += 1

    print("=== DOWNLOADING FULL-TEXT XMLS ===", flush=True)
    for pmcid in pmcids:
        st = _fetch_xml_direct(pmcid, ft_dir / f"{pmcid}.xml")
        print(f"  {pmcid}: {st}", flush=True)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
