"""Checksums, file validation, subprocess runner, and reporting utilities."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Iterable


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def compute_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, header: Iterable[str], rows: Iterable[Iterable[object]]) -> Path:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(list(header))
        for r in rows:
            writer.writerow(list(r))
    return path
