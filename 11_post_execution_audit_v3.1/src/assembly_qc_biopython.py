"""Stage 3: Dedicated Biopython-Based Assembly Quality Control Module for BRD v3.1.

Computes deterministic per-contig and per-assembly metrics using Biopython SeqIO
without invoking or imitating QUAST. Includes synthetic-assembly N50/L50 calculation,
FASTA integrity checks (empty files, duplicate headers, duplicate contig sequences,
invalid nucleotide characters, short-assembly warnings), and visualization plots.
"""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

try:
    from .paths_v3_1 import resolve_repo_root, resolve_v31_root
except ImportError:
    from paths_v3_1 import resolve_repo_root, resolve_v31_root

VALID_IUPAC_NUCL = set("ACGTNRYSWKMBDHV")
STANDARD_ACGT = set("ACGT")
DEFAULT_SHORT_ASSEMBLY_BP = 25_000_000


def compute_n50_l50(lengths: Sequence[int]) -> Tuple[int, int]:
    """Compute exact N50 and L50 from contig lengths (sorted descending).

    Returns (0, 0) for empty input or zero total length.
    """
    pos_lengths = sorted((int(x) for x in lengths if int(x) > 0), reverse=True)
    if not pos_lengths:
        return 0, 0
    total = sum(pos_lengths)
    threshold = 0.5 * total
    acc = 0
    for idx, ln in enumerate(pos_lengths, start=1):
        acc += ln
        if acc >= threshold:
            return ln, idx
    return pos_lengths[-1], len(pos_lengths)


def analyze_fasta_assembly(
    fasta_path: Path,
    isolate_id: str = "",
    assembly_accession: str = "",
    min_expected_bp: int = DEFAULT_SHORT_ASSEMBLY_BP,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Analyze a single assembly FASTA file and return summary, contig stats, and warnings."""
    warnings_list: List[Dict[str, Any]] = []
    contig_rows: List[Dict[str, Any]] = []

    if not fasta_path.exists():
        warnings_list.append(
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "warning_type": "MISSING_ASSEMBLY",
                "severity": "CRITICAL",
                "detail": f"Expected assembly FASTA missing at {fasta_path}",
            }
        )
        return (
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "fasta_file": fasta_path.name,
                "status": "MISSING",
                "total_length_bp": 0,
                "contig_count": 0,
                "contigs_gt_500bp": 0,
                "contigs_gt_1kb": 0,
                "contigs_gt_5kb": 0,
                "contigs_gt_10kb": 0,
                "n50_bp": 0,
                "l50": 0,
                "largest_contig_bp": 0,
                "smallest_contig_bp": 0,
                "mean_contig_bp": 0.0,
                "median_contig_bp": 0.0,
                "gc_percent": 0.0,
                "gc_percent_ungapped_acgt": 0.0,
                "n_base_count": 0,
                "other_iupac_ambiguous_count": 0,
                "ambiguous_nucleotide_count": 0,
                "ambiguous_nucleotide_percent": 0.0,
                "duplicate_id_count": 0,
                "duplicate_contig_seq_count": 0,
                "invalid_char_count": 0,
                "qc_engine": "Biopython_SeqIO",
            },
            contig_rows,
            warnings_list,
        )

    lengths: List[int] = []
    seen_ids: Dict[str, int] = {}
    seen_seq_hashes: Dict[str, List[str]] = {}
    dup_id_count = 0
    invalid_char_count = 0
    invalid_chars_seen: set[str] = set()
    total_gc = 0
    total_acgt = 0
    total_n = 0
    total_other_iupac = 0

    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(rec.seq).upper()
        ln = len(seq)
        lengths.append(ln)
        seen_ids[rec.id] = seen_ids.get(rec.id, 0) + 1
        if seen_ids[rec.id] > 1:
            dup_id_count += 1

        seq_md5 = hashlib.md5(seq.encode("ascii", errors="replace")).hexdigest()
        seen_seq_hashes.setdefault(seq_md5, []).append(rec.id)

        a_cnt = seq.count("A")
        t_cnt = seq.count("T")
        g_cnt = seq.count("G")
        c_cnt = seq.count("C")
        n_cnt = seq.count("N")
        acgt_cnt = a_cnt + t_cnt + g_cnt + c_cnt
        gc_cnt = g_cnt + c_cnt

        char_set = set(seq)
        bad_chars = char_set - VALID_IUPAC_NUCL
        rec_invalid_cnt = sum(seq.count(ch) for ch in bad_chars) if bad_chars else 0
        if bad_chars:
            invalid_chars_seen.update(bad_chars)
            invalid_char_count += rec_invalid_cnt

        other_iupac = max(0, ln - acgt_cnt - n_cnt - rec_invalid_cnt)
        rec_ambig = n_cnt + other_iupac

        total_gc += gc_cnt
        total_acgt += acgt_cnt
        total_n += n_cnt
        total_other_iupac += other_iupac

        contig_rows.append(
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "contig_id": rec.id,
                "length_bp": ln,
                "gc_percent_total": round(gc_cnt / ln * 100.0, 4) if ln > 0 else 0.0,
                "gc_percent_acgt": round(gc_cnt / acgt_cnt * 100.0, 4) if acgt_cnt > 0 else 0.0,
                "n_base_count": n_cnt,
                "ambiguous_nucleotide_count": rec_ambig,
                "invalid_char_count": rec_invalid_cnt,
                "seq_md5": seq_md5,
            }
        )

    if not lengths:
        warnings_list.append(
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "warning_type": "EMPTY_FASTA",
                "severity": "CRITICAL",
                "detail": f"Zero FASTA records parsed in {fasta_path.name}",
            }
        )
        return (
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "fasta_file": fasta_path.name,
                "status": "EMPTY",
                "total_length_bp": 0,
                "contig_count": 0,
                "contigs_gt_500bp": 0,
                "contigs_gt_1kb": 0,
                "contigs_gt_5kb": 0,
                "contigs_gt_10kb": 0,
                "n50_bp": 0,
                "l50": 0,
                "largest_contig_bp": 0,
                "smallest_contig_bp": 0,
                "mean_contig_bp": 0.0,
                "median_contig_bp": 0.0,
                "gc_percent": 0.0,
                "gc_percent_ungapped_acgt": 0.0,
                "n_base_count": 0,
                "other_iupac_ambiguous_count": 0,
                "ambiguous_nucleotide_count": 0,
                "ambiguous_nucleotide_percent": 0.0,
                "duplicate_id_count": 0,
                "duplicate_contig_seq_count": 0,
                "invalid_char_count": 0,
                "qc_engine": "Biopython_SeqIO",
            },
            contig_rows,
            warnings_list,
        )

    total_len = int(sum(lengths))
    ambig_total = total_n + total_other_iupac
    n50, l50 = compute_n50_l50(lengths)
    dup_seq_count = sum(len(ids) - 1 for ids in seen_seq_hashes.values() if len(ids) > 1)

    if dup_id_count > 0:
        warnings_list.append(
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "warning_type": "DUPLICATE_SEQUENCE_IDS",
                "severity": "HIGH",
                "detail": f"{dup_id_count} duplicate FASTA identifier(s) detected",
            }
        )
    if dup_seq_count > 0:
        warnings_list.append(
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "warning_type": "DUPLICATE_CONTIG_SEQUENCES",
                "severity": "MEDIUM",
                "detail": f"{dup_seq_count} identical contig sequence(s) found by MD5",
            }
        )
    if invalid_char_count > 0:
        warnings_list.append(
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "warning_type": "INVALID_NUCLEOTIDE_CHARACTERS",
                "severity": "HIGH",
                "detail": f"{invalid_char_count} invalid chars ({''.join(sorted(invalid_chars_seen))})",
            }
        )
    if total_len < min_expected_bp:
        warnings_list.append(
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "warning_type": "SUSPICIOUSLY_SHORT_ASSEMBLY",
                "severity": "HIGH",
                "detail": f"Total length {total_len:,} bp is below threshold {min_expected_bp:,} bp",
            }
        )
    ambig_pct = (ambig_total / total_len * 100.0) if total_len > 0 else 0.0
    if ambig_pct > 1.0:
        warnings_list.append(
            {
                "isolate_id": isolate_id or fasta_path.stem,
                "assembly_accession": assembly_accession,
                "warning_type": "ELEVATED_AMBIGUOUS_N_PCT",
                "severity": "MEDIUM",
                "detail": f"Ambiguous nucleotides represent {ambig_pct:.4f}% ({ambig_total:,} bp) of assembly (scaffold gaps)",
            }
        )

    arr = np.asarray(lengths, dtype=np.int64)
    summary = {
        "isolate_id": isolate_id or fasta_path.stem,
        "assembly_accession": assembly_accession,
        "fasta_file": fasta_path.name,
        "status": "OK",
        "total_length_bp": total_len,
        "contig_count": len(lengths),
        "contigs_gt_500bp": int((arr > 500).sum()),
        "contigs_gt_1kb": int((arr > 1000).sum()),
        "contigs_gt_5kb": int((arr > 5000).sum()),
        "contigs_gt_10kb": int((arr > 10000).sum()),
        "n50_bp": int(n50),
        "l50": int(l50),
        "largest_contig_bp": int(arr.max()),
        "smallest_contig_bp": int(arr.min()),
        "mean_contig_bp": round(float(arr.mean()), 2),
        "median_contig_bp": round(float(np.median(arr)), 2),
        "gc_percent": round(total_gc / total_len * 100.0, 4) if total_len > 0 else 0.0,
        "gc_percent_ungapped_acgt": round(total_gc / total_acgt * 100.0, 4) if total_acgt > 0 else 0.0,
        "n_base_count": int(total_n),
        "other_iupac_ambiguous_count": int(total_other_iupac),
        "ambiguous_nucleotide_count": int(ambig_total),
        "ambiguous_nucleotide_percent": round(ambig_pct, 4),
        "duplicate_id_count": int(dup_id_count),
        "duplicate_contig_seq_count": int(dup_seq_count),
        "invalid_char_count": int(invalid_char_count),
        "qc_engine": f"Biopython_SeqIO_{SeqIO.__package__}",
    }
    return summary, contig_rows, warnings_list


def run_stage3_assembly_qc(repo_root: Path | None = None) -> Dict[str, Path]:
    """Run Biopython assembly QC across all 19 authoritative B. sorokiniana isolates."""
    root = resolve_repo_root(repo_root)
    v31 = resolve_v31_root(root)
    s3 = v31 / "stage3_assembly_qc"
    plot_dir = s3 / "assembly_qc_plots_v3_1"
    s3.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    auth_iso = pd.read_csv(root / "11_post_execution_audit_v3" / "authoritative_data" / "auth_isolates.tsv", sep="\t")
    summaries: List[Dict[str, Any]] = []
    all_contigs: List[Dict[str, Any]] = []
    all_warnings: List[Dict[str, Any]] = []

    for _, row in auth_iso.iterrows():
        rel_path = str(row["relative_genome_path"])
        fpath = root / rel_path
        summ, c_rows, w_rows = analyze_fasta_assembly(
            fpath,
            isolate_id=str(row["isolate_id"]),
            assembly_accession=str(row["assembly_accession"]),
        )
        summ["assembly_level"] = str(row.get("assembly_level", ""))
        summ["v3_quality_flag"] = str(row.get("assembly_quality_flag", ""))
        summaries.append(summ)
        all_contigs.extend(c_rows)
        all_warnings.extend(w_rows)

    sum_df = pd.DataFrame(summaries).sort_values(["n50_bp", "total_length_bp"], ascending=[False, False]).reset_index(drop=True)
    contig_df = pd.DataFrame(all_contigs)
    warn_df = pd.DataFrame(all_warnings) if all_warnings else pd.DataFrame(
        columns=["isolate_id", "assembly_accession", "warning_type", "severity", "detail"]
    )

    sum_tsv = s3 / "assembly_qc_summary_v3_1.tsv"
    contig_tsv = s3 / "assembly_contig_statistics_v3_1.tsv"
    warn_tsv = s3 / "assembly_qc_warnings_v3_1.tsv"
    sum_df.to_csv(sum_tsv, sep="\t", index=False)
    contig_df.to_csv(contig_tsv, sep="\t", index=False)
    warn_df.to_csv(warn_tsv, sep="\t", index=False)

    # Visualization plots
    if not sum_df.empty and (sum_df["status"] == "OK").any():
        ok_df = sum_df[sum_df["status"] == "OK"].copy()
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

        y_pos = np.arange(len(ok_df))
        axes[0].barh(y_pos, ok_df["total_length_bp"] / 1e6, color="#0284c7")
        axes[0].set_yticks(y_pos)
        axes[0].set_yticklabels(ok_df["isolate_id"], fontsize=8.5)
        axes[0].invert_yaxis()
        axes[0].set_xlabel("Total Assembly Length (Mb)")
        axes[0].set_title("Assembly Span (19 B. sorokiniana Isolates)")
        axes[0].grid(axis="x", alpha=0.3)

        axes[1].barh(y_pos, ok_df["n50_bp"] / 1e6, color="#0d9488")
        axes[1].set_yticks(y_pos)
        axes[1].set_yticklabels([])
        axes[1].invert_yaxis()
        axes[1].set_xlabel("N50 (Mb)")
        axes[1].set_title("Contiguity (N50 & Contig Count)")
        for idx, (_, rr) in enumerate(ok_df.iterrows()):
            axes[1].text(rr["n50_bp"] / 1e6 + 0.03, idx, f"n={int(rr['contig_count'])}, L50={int(rr['l50'])}", va="center", fontsize=7.5)
        axes[1].grid(axis="x", alpha=0.3)

        axes[2].barh(y_pos, ok_df["ambiguous_nucleotide_percent"], color="#d97706")
        axes[2].set_yticks(y_pos)
        axes[2].set_yticklabels([])
        axes[2].invert_yaxis()
        axes[2].set_xlabel("Ambiguous Bases (% N + IUPAC)")
        axes[2].set_title("Ambiguous Gap Fraction")
        axes[2].grid(axis="x", alpha=0.3)

        fig.suptitle("BRD v3.1 Biopython Assembly Quality Control (Not QUAST)", fontsize=12, fontweight="bold")
        fig.tight_layout()
        fig.savefig(plot_dir / "assembly_qc_overview_v3_1.png", dpi=200)
        plt.close(fig)

        # Contig length distribution (log10 bp boxplot across isolates)
        fig2, ax2 = plt.subplots(figsize=(13, 5))
        order = ok_df["isolate_id"].tolist()
        box_data = [np.log10(np.maximum(1, contig_df[contig_df["isolate_id"] == iso]["length_bp"].to_numpy())) for iso in order]
        ax2.boxplot(box_data, tick_labels=order, showfliers=True)
        ax2.axhline(np.log10(10_000), color="#64748b", linestyle="--", linewidth=0.9, label="10 kb")
        ax2.axhline(np.log10(1_000), color="#94a3b8", linestyle=":", linewidth=0.9, label="1 kb")
        ax2.set_ylabel("Contig Length (log10 bp)")
        ax2.set_title("Per-Isolate Contig Length Distribution (Biopython SeqIO)")
        ax2.tick_params(axis="x", rotation=45, labelsize=8.5)
        ax2.legend(loc="lower left", fontsize=8)
        fig2.tight_layout()
        fig2.savefig(plot_dir / "contig_length_distribution_v3_1.png", dpi=200)
        plt.close(fig2)

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    methods_md = f"""# Stage 3 — Biopython Assembly Quality Control Methodology (BRD v3.1)

**Executed:** {now_utc}  
**Execution Route:** Route A — Native Windows 11 (`Biopython SeqIO`)  
**Tool Replacement Notice:** This module is a dedicated **Biopython** Python assembly statistics and FASTA integrity analyzer. **QUAST was NOT executed** (`Linux/WSL2` absent) and these outputs are never reported as QUAST results.

## Metrics & Integrity Checks Computed
- **Per-assembly size & count metrics:** `total_length_bp`, `contig_count`, strict length thresholds `contigs_gt_500bp`, `contigs_gt_1kb`, `contigs_gt_5kb`, `contigs_gt_10kb`, `largest_contig_bp`, `smallest_contig_bp`, `mean_contig_bp`, `median_contig_bp`.
- **Contiguity metrics:** Exact `N50` and `L50` computed via `compute_n50_l50()` (validated on synthetic unit test assemblies).
- **Base composition & ambiguity:** `gc_percent` (total span), `gc_percent_ungapped_acgt` (`G+C / (A+C+G+T)` excluding `N` gaps), `n_base_count`, `other_iupac_ambiguous_count`, `ambiguous_nucleotide_count`, `ambiguous_nucleotide_percent`.
- **Integrity alerts (`assembly_qc_warnings_v3_1.tsv`):** Detects missing files, empty FASTAs, duplicate sequence IDs, MD5-identical duplicate contig records, invalid non-IUPAC characters, assemblies `<25 Mb`, and scaffold gap fractions `>1.0%`.

## Summary Snapshot Across All {len(sum_df)} Isolates
- **Available & parsed:** `{int((sum_df['status'] == 'OK').sum())}/{len(sum_df)}`
- **Warnings raised:** `{len(warn_df)}` (`{'; '.join(sorted(warn_df['warning_type'].unique())) if not warn_df.empty else 'None'}`)

## Explicit Scientific Limitations (Not QUAST)
Biopython assembly metrics measure only intrinsic FASTA composition and contiguity. They **do not** independently establish:
1. Reference-based structural misassemblies, translocations, or relocations.
2. Biological accuracy of local consensus bases (e.g., long-read homopolymer indels in `WAI2406` / `WAI2411`).
3. Chromosome-scale completeness or absence of chimeric scaffolding.
4. Complete gene space representation or annotation completeness.
"""
    methods_path = s3 / "assembly_qc_methods.md"
    methods_path.write_text(methods_md, encoding="utf-8")

    return {
        "summary": sum_tsv,
        "contigs": contig_tsv,
        "warnings": warn_tsv,
        "methods": methods_path,
        "plots_dir": plot_dir,
    }


if __name__ == "__main__":
    for k, v in run_stage3_assembly_qc().items():
        print(f"STAGE3_ASM_OUT {k}: {v}")
