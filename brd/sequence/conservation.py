"""Nucleotide-level multi-isolate conservation, SNP/indel detection, and 21-mer preservation."""

from __future__ import annotations

from typing import Dict, List, Tuple
import pandas as pd
from Bio.Align import MultipleSeqAlignment


def compute_pairwise_identity_and_coverage(seq_a: str, seq_b: str) -> Tuple[float, float]:
    """Compute pairwise nucleotide identity (%) and alignment coverage (%) from aligned strings."""
    aligned_cols = 0
    matches = 0
    non_gap_a = 0
    non_gap_b = 0
    for a, b in zip(seq_a.upper(), seq_b.upper()):
        if a != "-":
            non_gap_a += 1
        if b != "-":
            non_gap_b += 1
        if a != "-" and b != "-":
            aligned_cols += 1
            if a == b:
                matches += 1
    ident = (matches / aligned_cols * 100.0) if aligned_cols > 0 else 0.0
    denom = max(non_gap_a, non_gap_b, 1)
    cov = aligned_cols / denom * 100.0
    return round(ident, 3), round(cov, 3)


def analyze_msa_conservation(
    msa: MultipleSeqAlignment,
    ref_id: str,
    window_size: int = 21,
) -> Dict[str, object]:
    """
    Map MSA columns onto reference sequence coordinates (1-based),
    detect SNPs, indels, per-position conservation fraction, conserved blocks,
    and 21-mer exact conservation across all aligned isolates.
    """
    seqs = {rec.id: str(rec.seq).upper() for rec in msa}
    if ref_id not in seqs:
        ref_id = msa[0].id
    ref_aln = seqs[ref_id]
    other_ids = [rid for rid in seqs if rid != ref_id]
    n_isolates = len(seqs)

    # Map MSA index -> reference ungapped index (0-based)
    pos_rows = []
    ref_pos = -1
    ungapped_chars: Dict[str, List[str]] = {rid: [] for rid in seqs}

    snp_count = 0
    indel_events = 0
    prev_indel = False
    pairwise_stats = []

    for rid in other_ids:
        pid, cov = compute_pairwise_identity_and_coverage(ref_aln, seqs[rid])
        pairwise_stats.append({"isolate": rid, "pairwise_identity_pct": pid, "alignment_coverage_pct": cov})

    ref_ungapped_cols = []
    for col_idx, ref_base in enumerate(ref_aln):
        col_bases = [seqs[rid][col_idx] for rid in seqs]
        if ref_base == "-":
            if not prev_indel:
                indel_events += 1
                prev_indel = True
            continue
        ref_pos += 1
        prev_indel = False
        bases_non_gap = [b for b in col_bases if b != "-"]
        gaps = sum(1 for b in col_bases if b == "-")
        if gaps > 0:
            indel_events += 1
        matches_ref = sum(1 for b in col_bases if b == ref_base)
        cons_frac = matches_ref / n_isolates
        unique_non_gap = set(bases_non_gap)
        is_snp = len(unique_non_gap) > 1
        if is_snp:
            snp_count += 1
        pos_rows.append(
            {
                "ref_pos_1based": ref_pos + 1,
                "ref_base": ref_base,
                "conservation_frac": round(cons_frac, 4),
                "n_isolates_matching_ref": matches_ref,
                "n_gaps": gaps,
                "is_snp": is_snp,
                "observed_alleles": ",".join(sorted(unique_non_gap)),
            }
        )
        ref_ungapped_cols.append(col_idx)

    # Identify contiguous 100%-conserved blocks (>= window_size)
    conserved_blocks = []
    in_block = False
    b_start = 0
    for i, pr in enumerate(pos_rows):
        fully_cons = (pr["conservation_frac"] == 1.0 and pr["n_gaps"] == 0)
        if fully_cons and not in_block:
            in_block = True
            b_start = i + 1
        elif not fully_cons and in_block:
            in_block = False
            b_end = i
            if (b_end - b_start + 1) >= window_size:
                conserved_blocks.append((b_start, b_end, b_end - b_start + 1))
    if in_block:
        b_end = len(pos_rows)
        if (b_end - b_start + 1) >= window_size:
            conserved_blocks.append((b_start, b_end, b_end - b_start + 1))

    # Compute exact k-mer conservation fraction for every window_size-mer along ref
    ref_ungapped_seq = ref_aln.replace("-", "")
    other_ungapped = {rid: seqs[rid].replace("-", "") for rid in other_ids}
    kmer_exact_fractions = []
    for i in range(max(0, len(ref_ungapped_seq) - window_size + 1)):
        kmer = ref_ungapped_seq[i : i + window_size]
        hits = 1 + sum(1 for rid in other_ids if kmer in other_ungapped[rid])
        kmer_exact_fractions.append(hits / n_isolates)

    mean_pairwise_id = (
        sum(p["pairwise_identity_pct"] for p in pairwise_stats) / len(pairwise_stats) if pairwise_stats else 100.0
    )
    min_pairwise_id = min((p["pairwise_identity_pct"] for p in pairwise_stats), default=100.0)
    mean_kmer_cons = (
        sum(kmer_exact_fractions) / len(kmer_exact_fractions) * 100.0 if kmer_exact_fractions else 0.0
    )

    return {
        "n_isolates": n_isolates,
        "ref_length_nt": len(ref_ungapped_seq),
        "mean_pairwise_identity_pct": round(mean_pairwise_id, 3),
        "min_pairwise_identity_pct": round(min_pairwise_id, 3),
        "snp_sites": snp_count,
        "indel_columns": indel_events,
        "conserved_blocks_ge21nt": conserved_blocks,
        "mean_21mer_isolate_conservation_pct": round(mean_kmer_cons, 3),
        "pct_100pct_conserved_21mers": round(
            sum(1 for f in kmer_exact_fractions if f == 1.0) / max(len(kmer_exact_fractions), 1) * 100.0, 2
        ),
        "pairwise_stats": pairwise_stats,
        "per_position": pd.DataFrame(pos_rows),
        "kmer_conservation_fractions": kmer_exact_fractions,
    }
