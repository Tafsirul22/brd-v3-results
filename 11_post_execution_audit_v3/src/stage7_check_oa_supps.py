"""Check Europe PMC supplementaryFiles ZIP archive for primary benchmark papers."""

from __future__ import annotations

import io
import sys
import urllib.request
import zipfile
from pathlib import Path

from paths_v3 import resolve_v3_root


def main() -> None:
    supp_dir = resolve_v3_root() / "stage7_literature_benchmark" / "supplements"
    supp_dir.mkdir(parents=True, exist_ok=True)
    for pid in ("PMC5943259", "PMC5063301", "PMC7202221", "PMC6465682", "PMC13608320", "PMC13205770"):
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pid}/supplementaryFiles"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "BRD-v3/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = resp.read()
            zf = zipfile.ZipFile(io.BytesIO(data))
            names = zf.namelist()
            out_sub = supp_dir / pid
            out_sub.mkdir(exist_ok=True)
            zf.extractall(out_sub)
            print(f"{pid}: OK ({len(data)} B) files={names}")
        except Exception as exc:
            print(f"{pid}: UNAVAIL ({exc})")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
