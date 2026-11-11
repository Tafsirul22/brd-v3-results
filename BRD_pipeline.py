"""Bipolaris-RNAi Designer (BRD) — Master Command-Line Interface."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List

from brd import __version__
from brd.config import BRDConfig
from brd.logger import get_logger
from brd.modules.m01_acquisition import run_m01
from brd.modules.m02_genome_qc import run_m02
from brd.modules.m03_annotation import run_m03
from brd.modules.m04_comparative_genomics import run_m04
from brd.modules.m05_functional_evidence import run_m05
from brd.modules.m06_target_prioritisation import run_m06
from brd.modules.m07_dsRNA_design import run_m07
from brd.modules.m08_off_target import run_m08
from brd.modules.m09_ranking import run_m09
from brd.modules.m10_reporting import run_m10

MODULES: Dict[str, tuple[str, Callable[[BRDConfig], dict], List[str]]] = {
    "m01": ("M01 Data acquisition", run_m01, ["01_data_acquisition/isolate_metadata.tsv", "01_data_acquisition/download_manifest.tsv"]),
    "m02": ("M02 Genome QC", run_m02, ["02_genome_qc/assembly_quality_summary.tsv"]),
    "m03": ("M03 Annotation", run_m03, ["03_annotation/harmonised/gene_id_mapping.tsv"]),
    "m04": ("M04 Comparative genomics", run_m04, ["04_comparative_genomics/nucleotide_conservation/target_nucleotide_conservation.tsv"]),
    "m05": ("M05 Functional evidence", run_m05, ["05_functional_evidence/functional_evidence_summary/functional_evidence_summary.tsv"]),
    "m06": ("M06 Target prioritisation", run_m06, ["06_target_prioritisation/target_scoring/target_prioritisation.tsv"]),
    "m07": ("M07 dsRNA design", run_m07, ["07_dsRNA_design/candidate_sequences/dsRNA_candidates.tsv"]),
    "m08": ("M08 Off-target screening", run_m08, ["08_off_target_screening/summaries/off_target_summary.tsv"]),
    "m09": ("M09 Candidate ranking", run_m09, ["09_candidate_ranking/integrated_scores/final_candidate_ranking.tsv"]),
    "m10": ("M10 Reporting", run_m10, ["10_final_reports/BRD_Final_Research_Report.html", "10_final_reports/BRD_Execution_Summary.txt"]),
}


def _state_file(cfg: BRDConfig) -> Path:
    return cfg.output_dir / "logs" / "pipeline_state.json"


def _load_state(cfg: BRDConfig) -> dict:
    sf = _state_file(cfg)
    if sf.exists():
        return json.loads(sf.read_text(encoding="utf-8"))
    return {}


def _save_state(cfg: BRDConfig, state: dict) -> None:
    sf = _state_file(cfg)
    sf.parent.mkdir(parents=True, exist_ok=True)
    sf.write_text(json.dumps(state, indent=2), encoding="utf-8")


def print_status(cfg: BRDConfig) -> None:
    state = _load_state(cfg)
    print(f"Bipolaris-RNAi Designer (BRD) v{__version__}")
    print(f"Output directory: {cfg.output_dir}")
    print("\n--- Dependency Status ---")
    for key in ["blastn", "makeblastdb", "diamond", "datasets", "mafft", "dsrnamax"]:
        tp = cfg.tool_path(key)
        ok = Path(tp).exists() or shutil.which(tp) is not None
        print(f"  {key:<12}: {'FOUND' if ok else 'MISSING'} ({tp})")
    for ltool in ["quast", "busco", "orthofinder", "funannotate", "braker.pl", "RNAhybrid"]:
        print(f"  {ltool:<12}: {'FOUND' if shutil.which(ltool) else 'NOT_FOUND (Linux-only; handled via documented fallback/NA)'}")

    print("\n--- Module Status ---")
    for mkey, (label, _, outs) in MODULES.items():
        all_exist = all((cfg.output_dir / o).exists() and (cfg.output_dir / o).stat().st_size > 0 for o in outs)
        st = state.get(mkey, {}).get("status", "COMPLETED" if all_exist else "PENDING")
        print(f"  {label:<28}: {st:<10} -> {', '.join(outs)}")


def run_modules(cfg: BRDConfig, keys: List[str], resume: bool = False) -> int:
    logger = get_logger("BRD.CLI", cfg.output_dir / "logs" / "brd_execution.log")
    state = _load_state(cfg)
    for mkey in keys:
        label, fn, outs = MODULES[mkey]
        all_exist = all((cfg.output_dir / o).exists() and (cfg.output_dir / o).stat().st_size > 0 for o in outs)
        if resume and state.get(mkey, {}).get("status") == "COMPLETED" and all_exist:
            logger.info("Skipping %s (already COMPLETED with verified outputs).", label)
            continue
        logger.info("Starting %s...", label)
        try:
            fn(cfg)
            state[mkey] = {"status": "COMPLETED", "timestamp_utc": datetime.now(timezone.utc).isoformat()}
            _save_state(cfg, state)
        except Exception as exc:
            logger.exception("%s FAILED: %s", label, exc)
            state[mkey] = {"status": "FAILED", "error": str(exc), "timestamp_utc": datetime.now(timezone.utc).isoformat()}
            _save_state(cfg, state)
            return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="BRD_pipeline.py",
        description="Bipolaris-RNAi Designer (BRD) — Complete M01–M10 Fungal RNAi Computational Pipeline",
    )
    parser.add_argument("--config", type=str, default=None, help="Path to project_configuration.yaml")
    parser.add_argument("--run", type=str, default=None, help="Run 'all' or a specific module: m01..m10")
    parser.add_argument("--resume", action="store_true", help="Resume pipeline, skipping completed modules")
    parser.add_argument("--status", action="store_true", help="Report module status, dependencies, and outputs")
    args = parser.parse_args()

    cfg = BRDConfig.load(Path(args.config) if args.config else None)
    if args.status:
        print_status(cfg)
        return 0
    if args.resume:
        return run_modules(cfg, list(MODULES.keys()), resume=True)
    if args.run:
        target = args.run.lower().strip()
        if target == "all":
            return run_modules(cfg, list(MODULES.keys()), resume=False)
        if target in MODULES:
            return run_modules(cfg, [target], resume=False)
        parser.error(f"Unknown module '{args.run}'. Valid choices: all, {', '.join(MODULES.keys())}")
    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
