"""Search Stage 7 XMLs for inline 5'-3' RNA/dsRNA/siRNA sequences & test A. alternata / M. oryzae / P. pachyrhizi GenBank reconstructions."""

from __future__ import annotations

import json
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio.Seq import Seq
from paths_v3 import resolve_v3_root


def main() -> None:
    v3 = resolve_v3_root()
    ft_dir = v3 / "stage7_literature_benchmark" / "fulltexts"
    ref_dir = v3 / "stage7_literature_benchmark" / "reconstructed_refs"
    ref_dir.mkdir(parents=True, exist_ok=True)

    # 1. Check for any inline siRNA/dsRNA sequences across all XMLs (allow RNA U and unicode dashes)
    for pid in ("PMC6465682", "PMC5040644", "PMC8428832", "PMC7072719", "PMC13608320", "PMC13205770", "PMC8710829", "PMC7202221"):
        xp = ft_dir / f"{pid}.xml"
        if not xp.exists():
            continue
        txt = " ".join("".join(ET.parse(xp).getroot().itertext()).split())
        hits = []
        for m in re.finditer(r"5['′\u2032\s-]*([ACGUTacgut\s-]{18,90})['′\u2032\s-]*3", txt):
            seq = re.sub(r"[^ACGUTacgut]", "", m.group(1)).upper()
            if len(seq) >= 18 and len(set(seq)) >= 3:
                ctx = txt[max(0, m.start() - 55) : min(len(txt), m.end() + 35)]
                hits.append((seq, ctx))
        if hits:
            print(f"INLINE {pid}: {len(hits)}")
            for s, c in hits[:6]:
                print(f"  {len(s)}nt {s} <-- {c}")

    # 2. Check Alternaria alternata pksJ (PMC3391263) + Botrytis cinerea erg11/erg13/erg1 (PMC8710829) + M. oryzae MoAP1 (PMC6465682) in NCBI nuccore
    lookups = [
        ("A_alternata_pksJ", "Alternaria alternata[Organism] AND pksJ[Gene]", "GACTCGATACAGCCATTGCAG", "GTGCCTCCATACCCGAATGAG"),
        ("M_oryzae_MoAP1", "Magnaporthe oryzae[Organism] AND MoAP1[Title]", "", ""),
    ]
    for tag, term, fwd, rev in lookups:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "term": term, "retmode": "json", "retmax": 3})
        try:
            ids = json.loads(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "BRD-v3/1.0"}), timeout=15).read().decode("utf-8"))["esearchresult"]["idlist"]
            print(f"ESEARCH {tag}: {ids}")
            if ids and fwd and rev:
                fa = urllib.request.urlopen(urllib.request.Request(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(ids)}&rettype=fasta&retmode=text", headers={"User-Agent": "BRD-v3/1.0"}), timeout=15).read().decode("utf-8")
                (ref_dir / f"{tag}.fasta").write_text(fa, encoding="utf-8")
                rc = str(Seq(rev).reverse_complement())
                i1 = fa.upper().find(fwd)
                # Search per sequence clean
                from Bio import SeqIO
                for rec in SeqIO.parse(str(ref_dir / f"{tag}.fasta"), "fasta"):
                    sq = str(rec.seq).upper()
                    p1, p2 = sq.find(fwd), sq.find(rc)
                    if p1 != -1 and p2 != -1 and p2 > p1:
                        print(f"  MATCH {tag} on {rec.id}: {p1+1}..{p2+len(rc)} ({p2+len(rc)-p1} bp)")
        except Exception as exc:
            print(f"LOOKUP_ERR {tag}: {exc}")

    # 3. Also print verbatim sentences containing asiR1115 in PMC6465682 or check if PMC4515635/PMC6764512 have tables
    for pid in ("PMC6465682",):
        txt = " ".join("".join(ET.parse(ft_dir / f"{pid}.xml").getroot().itertext()).split())
        for m in re.finditer(r"asiR1115|asiR1245|asiR1362", txt):
            print("CONTEXT:", txt[max(0, m.start() - 40) : m.end() + 90])


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
