"""Subprocess wrapper capturing command, stdout, stderr, return code, and runtime."""

from __future__ import annotations

import json
import subprocess
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class CommandResult:
    tool: str
    command: List[str]
    returncode: int
    duration_sec: float
    stdout: str
    stderr: str
    timestamp_utc: str
    cwd: str


def run_tool(
    tool: str,
    cmd: List[str],
    cwd: Optional[Path] = None,
    log_path: Optional[Path] = None,
    check: bool = True,
    timeout: int = 3600,
) -> CommandResult:
    t0 = time.perf_counter()
    ts = datetime.now(timezone.utc).isoformat()
    proc = subprocess.run(
        [str(c) for c in cmd],
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )
    dt = round(time.perf_counter() - t0, 3)
    res = CommandResult(
        tool=tool,
        command=[str(c) for c in cmd],
        returncode=proc.returncode,
        duration_sec=dt,
        stdout=proc.stdout,
        stderr=proc.stderr,
        timestamp_utc=ts,
        cwd=str(cwd or Path.cwd()),
    )
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(asdict(res)) + "\n")
    if check and proc.returncode != 0:
        raise RuntimeError(
            f"Command failed ({tool}, rc={proc.returncode}): {' '.join(res.command)}\nSTDERR: {proc.stderr[:1000]}"
        )
    return res
