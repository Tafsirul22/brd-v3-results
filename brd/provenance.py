"""Data and execution provenance recorder for BRD."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

from brd.utils.checksums import compute_sha256


class ProvenanceTracker:
    """Append-only provenance ledger storing inputs, outputs, checksums, and tool versions."""

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir
        self.meta_dir = output_dir / "00_project_metadata"
        self.meta_dir.mkdir(parents=True, exist_ok=True)
        self.tsv_path = self.meta_dir / "data_provenance.tsv"
        self.json_path = self.meta_dir / "provenance.json"
        self.records: List[Dict[str, str]] = []
        if self.json_path.exists():
            try:
                self.records = json.loads(self.json_path.read_text(encoding="utf-8"))
            except Exception:
                self.records = []

    def record(
        self,
        module: str,
        step: str,
        source: str,
        output_file: Path,
        tool: str = "python",
        tool_version: str = "",
        accession_or_version: str = "",
        notes: str = "",
    ) -> None:
        sha = compute_sha256(output_file) if output_file.exists() and output_file.is_file() else ""
        size = str(output_file.stat().st_size) if output_file.exists() and output_file.is_file() else "0"
        entry = {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "module": module,
            "step": step,
            "source": source,
            "output_path": str(output_file),
            "size_bytes": size,
            "sha256": sha,
            "tool": tool,
            "tool_version": tool_version,
            "accession_or_version": accession_or_version,
            "notes": notes,
        }
        self.records.append(entry)
        self.flush()

    def flush(self) -> None:
        fields = [
            "timestamp_utc",
            "module",
            "step",
            "source",
            "output_path",
            "size_bytes",
            "sha256",
            "tool",
            "tool_version",
            "accession_or_version",
            "notes",
        ]
        with open(self.tsv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(self.records)
        self.json_path.write_text(json.dumps(self.records, indent=2), encoding="utf-8")
