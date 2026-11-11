"""Module M02 — Genome Quality Control."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from brd.config import BRDConfig
from brd.integrations.busco import get_busco_status
from brd.integrations.quast import assess_assembly
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker
from brd.utils.report_utils import df_to_html_table, wrap_html_report, write_text


def run_m02(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    cmd_log = cfg.output_dir / "logs" / "commands.jsonl"
    logger = get_logger("BRD.M02", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    m01_iso = cfg.output_dir / "01_data_acquisition" / "isolate_metadata.tsv"
    if not m01_iso.exists():
        raise FileNotFoundError(f"M01 output missing: {m01_iso}")

    m02_dir = cfg.module_dir("02_genome_qc")
    quast_dir = m02_dir / "quast"
    busco_dir = m02_dir / "busco"
    fig_dir = m02_dir / "figures"
    for d in (quast_dir, busco_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    df_iso = pd.read_csv(m01_iso, sep="\t")
    qc_rows: List[dict] = []
    busco_rows: List[dict] = []

    for _, row in df_iso.iterrows():
        acc = str(row["assembly_accession"])
        iso_id = str(row["isolate_id"])
        fna = Path(str(row["genome_fasta"]))
        if not fna.exists():
            logger.warning("Missing FASTA for %s: %s", acc, fna)
            continue

        stats = assess_assembly(fna, quast_dir / acc, log_path=cmd_log)
        busco_st = get_busco_status(acc, fna, lineage="fungi_odb10")
        busco_rows.append({"isolate_id": iso_id, **busco_st})

        # Quality decision criteria
        flags = []
        if stats["total_length_bp"] < 28_000_000 or stats["total_length_bp"] > 45_000_000:
            flags.append("unusual_genome_size")
        if stats["gc_percent"] < 46.0 or stats["gc_percent"] > 54.0:
            flags.append("unusual_gc_content")
        if stats["n_percent"] > 5.0:
            flags.append("high_ambiguous_bases")
        if stats["duplicate_headers"] > 0:
            flags.append("duplicate_headers")
        if stats["invalid_chars"]:
            flags.append(f"invalid_chars:{stats['invalid_chars']}")

        qc_status = "RETAIN" if not flags else "FLAG"
        qc_reason = "Passed FASTA validity, size (28-45 Mb), GC (46-54%), and N-content (<5%) checks" if not flags else "; ".join(flags)

        qc_rows.append(
            {
                "isolate_id": iso_id,
                "strain": row["strain"],
                "assembly_accession": acc,
                "assembly_level": row["assembly_level"],
                "geographic_origin": row["geographic_origin"],
                **stats,
                "qc_status": qc_status,
                "qc_reason": qc_reason,
                "busco_status": busco_st["status"],
            }
        )

    df_qc = pd.DataFrame(qc_rows)
    df_busco = pd.DataFrame(busco_rows)

    asm_stats_p = m02_dir / "assembly_statistics.tsv"
    asm_qual_p = m02_dir / "assembly_quality_summary.tsv"
    quast_sum_p = m02_dir / "quast_summary.tsv"
    busco_sum_p = m02_dir / "busco_summary.tsv"

    df_qc.to_csv(asm_stats_p, sep="\t", index=False)
    df_qc.to_csv(asm_qual_p, sep="\t", index=False)
    df_qc.to_csv(quast_sum_p, sep="\t", index=False)
    df_busco.to_csv(busco_sum_p, sep="\t", index=False)
    df_qc.to_csv(quast_dir / "quast_summary.tsv", sep="\t", index=False)
    df_busco.to_csv(busco_dir / "busco_summary.tsv", sep="\t", index=False)

    # Generate Genome QC Matplotlib Figures (PNG + SVG)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    df_sorted = df_qc.sort_values("n50_bp", ascending=True)
    y_labels = [f"{r['strain']} ({r['assembly_accession']})" for _, r in df_sorted.iterrows()]

    axes[0].barh(y_labels, df_sorted["total_length_bp"] / 1e6, color="#2563eb", alpha=0.88)
    axes[0].set_xlabel("Total Assembly Length (Mb)")
    axes[0].set_title("Assembly Size Across Isolates")
    axes[0].tick_params(axis="y", labelsize=8)

    axes[1].barh(y_labels, df_sorted["n50_bp"] / 1e6, color="#059669", alpha=0.88)
    axes[1].set_xlabel("Contig/Scaffold N50 (Mb)")
    axes[1].set_title("Contiguity (N50)")
    axes[1].set_yticks([])

    axes[2].scatter(df_sorted["gc_percent"], df_sorted["scaffold_count"], c="#d97706", s=60, alpha=0.85, edgecolors="#1e293b")
    axes[2].set_xlabel("GC Content (%)")
    axes[2].set_ylabel("Scaffold / Contig Count")
    axes[2].set_title("GC Content vs Fragmentation")
    axes[2].grid(True, linestyle="--", alpha=0.4)

    fig.suptitle(f"Bipolaris sorokiniana Isolate Panel Assembly QC (n = {len(df_qc)} biological isolates)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    png_path = fig_dir / "genome_quality_visualisation.png"
    svg_path = fig_dir / "genome_quality_visualisation.svg"
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path)
    plt.close(fig)

    # GC distribution plot
    fig2, ax2 = plt.subplots(figsize=(7, 4.5))
    ax2.hist(df_qc["gc_percent"], bins=10, color="#0284c7", edgecolor="#0f172a", alpha=0.85)
    ax2.axvline(df_qc["gc_percent"].mean(), color="#dc2626", linestyle="--", label=f"Mean GC = {df_qc['gc_percent'].mean():.2f}%")
    ax2.set_xlabel("GC Percentage (%)")
    ax2.set_ylabel("Number of Isolate Assemblies")
    ax2.set_title("GC% Distribution Across B. sorokiniana Assemblies")
    ax2.legend()
    fig2.tight_layout()
    fig2.savefig(fig_dir / "gc_distribution.png", dpi=300)
    fig2.savefig(fig_dir / "gc_distribution.svg")
    plt.close(fig2)

    body = (
        f"<div class='callout'><strong>Scope & Tool Provenance:</strong> Computed for {len(df_qc)} deduplicated biological "
        f"isolates using deterministic Biopython QUAST-compatible sequence inspection. Native Linux QUAST/BUSCO CLI binaries "
        f"were unavailable natively on Windows 11 without an installed WSL2 Linux distribution; BUSCO completeness fields are "
        f"explicitly marked as <code>NOT_RUN_LINUX_ONLY_NO_WSL_DISTRO</code> (not fabricated).</div>"
        f"<h2>Assembly Quality Metrics</h2>{df_to_html_table(df_qc[['isolate_id', 'assembly_accession', 'assembly_level', 'total_length_bp', 'scaffold_count', 'n50_bp', 'l50', 'gc_percent', 'n_percent', 'qc_status']])}"
        f"<h2>Visualisation</h2><figure><img src='figures/genome_quality_visualisation.png' alt='Genome QC'/><figcaption>Assembly length, N50 contiguity, and GC% vs scaffold fragmentation.</figcaption></figure>"
        f"<h2>BUSCO Execution Status</h2>{df_to_html_table(df_busco[['isolate_id', 'assembly_accession', 'lineage', 'status', 'notes']])}"
    )
    rep_html = m02_dir / "genome_quality_report.html"
    write_text(rep_html, wrap_html_report("BRD M02 — Genome Quality Control Report", "Bipolaris sorokiniana Multi-Isolate Panel", body))

    prov.record("M02", "Genome QC", str(m01_iso), asm_qual_p, "Biopython+QUAST-compat", "1.86", f"{len(df_qc)}_isolates")
    logger.info("M02 complete: %d assemblies quality-assessed.", len(df_qc))
    return {
        "assembly_quality_summary": asm_qual_p,
        "quast_summary": quast_sum_p,
        "busco_summary": busco_sum_p,
        "report": rep_html,
    }
