"""NCBI Datasets CLI, Datasets REST API, and Entrez E-utilities integration."""

from __future__ import annotations

import json
import time
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests
from Bio import Entrez

from brd.utils.subprocess_runner import run_tool


class NCBIClient:
    """Rate-limited NCBI Entrez + Datasets API/CLI wrapper."""

    def __init__(
        self,
        email: str,
        api_key: Optional[str] = None,
        datasets_bin: str = "datasets",
        rate_limit_sec: float = 0.35,
        log_path: Optional[Path] = None,
    ) -> None:
        self.email = email
        self.api_key = api_key
        self.datasets_bin = datasets_bin
        self.rate_limit_sec = 0.12 if api_key else rate_limit_sec
        self.log_path = log_path
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

    def _sleep(self) -> None:
        time.sleep(self.rate_limit_sec)

    def fetch_taxon_assemblies(self, taxid: int) -> List[Dict[str, Any]]:
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report?page_size=100"
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        return resp.json().get("reports", [])

    def download_genome_package(
        self,
        accession: str,
        include: List[str],
        dest_dir: Path,
        overwrite: bool = False,
    ) -> Path:
        dest_dir.mkdir(parents=True, exist_ok=True)
        marker = dest_dir / f".downloaded_{accession}"
        zip_path = dest_dir / f"{accession}.zip"
        if marker.exists() and not overwrite:
            return dest_dir
        inc_str = ",".join(include)
        cmd = [
            self.datasets_bin,
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            inc_str,
            "--no-progressbar",
            "--filename",
            str(zip_path),
        ]
        run_tool("datasets", cmd, log_path=self.log_path, timeout=1800)
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(dest_dir)
        zip_path.unlink(missing_ok=True)
        marker.write_text("ok\n", encoding="utf-8")
        return dest_dir

    def efetch_text(self, db: str, ids: List[str], rettype: str = "gb", retmode: str = "text") -> str:
        self._sleep()
        handle = Entrez.efetch(db=db, id=",".join(ids), rettype=rettype, retmode=retmode)
        data = handle.read()
        handle.close()
        return data.decode("utf-8", errors="replace") if isinstance(data, bytes) else data

    def esummary(self, db: str, ids: List[str]):
        self._sleep()
        handle = Entrez.esummary(db=db, id=",".join(ids))
        recs = Entrez.read(handle)
        handle.close()
        return recs if isinstance(recs, dict) else list(recs)

    def esearch(self, db: str, term: str, retmax: int = 20) -> Dict[str, Any]:
        self._sleep()
        handle = Entrez.esearch(db=db, term=term, retmax=retmax)
        res = Entrez.read(handle)
        handle.close()
        return dict(res)
