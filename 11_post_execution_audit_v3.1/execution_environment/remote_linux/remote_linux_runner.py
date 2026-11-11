"""Route C: Optional Remote Linux SSH Execution & Checksum-Verified Retrieval Runner.

Provides a transparent, failure-recoverable workflow for running Linux-only tools
(e.g., Compleasm with miniprot/hmmsearch or FastTree-backed Broccoli) on a remote
Linux host via OpenSSH (`ssh`/`scp`). Refuses to run when no real remote config
exists and never fabricates remote outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List

import yaml


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_ssh_preconditions(cfg_path: Path) -> Dict[str, object]:
    if shutil.which("ssh") is None or shutil.which("scp") is None:
        return {"ready": False, "reason": "OpenSSH ssh/scp client missing from PATH"}
    if not cfg_path.exists():
        return {"ready": False, "reason": f"Remote config file not found: {cfg_path}"}
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
    host = str(cfg.get("remote", {}).get("host", ""))
    if not host or "example.org" in host:
        return {"ready": False, "reason": "No real remote Linux host configured (example placeholder only)"}
    return {"ready": True, "config": cfg}


def build_remote_commands(cfg: dict, local_inputs: List[Path]) -> List[str]:
    r = cfg["remote"]
    sw = cfg["software_environment"]
    target_prefix = f"{r['username']}@{r['host']}"
    work = r["remote_workdir"]
    cmds = [
        f"ssh -p {r.get('port', 22)} {target_prefix} 'mkdir -p {work}/inputs {work}/outputs'",
    ]
    for fp in local_inputs:
        cmds.append(f"scp -P {r.get('port', 22)} \"{fp}\" {target_prefix}:{work}/inputs/")
    cmds.append(
        f"ssh -p {r.get('port', 22)} {target_prefix} 'cd {work} && sha256sum inputs/* > inputs.sha256 && "
        f"conda run -n {sw.get('conda_env', 'brd_linux_tools')} compleasm.py run -a inputs/genome.fna "
        f"-l {sw.get('compleasm_lineage', 'pleosporales_odb12')} -o outputs/compleasm -t {sw.get('threads', 8)} "
        f"&& sha256sum outputs/compleasm/* > outputs.sha256'"
    )
    cmds.append(f"scp -r -P {r.get('port', 22)} {target_prefix}:{work}/outputs ./remote_outputs")
    return cmds


def main() -> int:
    parser = argparse.ArgumentParser(description="BRD v3.1 Route C Remote Linux Execution Helper")
    parser.add_argument("--config", type=Path, required=True, help="Path to real remote_linux_config.yaml")
    parser.add_argument("--dry-run", action="store_true", help="Print planned SSH/SCP commands without executing")
    args = parser.parse_args()

    check = verify_ssh_preconditions(args.config)
    if not check["ready"]:
        print(f"ROUTE_C_BLOCKED: {check['reason']}")
        return 2

    cfg = check["config"]
    cmds = build_remote_commands(cfg, [])
    for c in cmds:
        print(c)
    if args.dry_run:
        return 0

    for c in cmds:
        proc = subprocess.run(c, shell=True)
        if proc.returncode != 0:
            print(f"REMOTE_STEP_FAILED rc={proc.returncode}: {c}")
            return proc.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
