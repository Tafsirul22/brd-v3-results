"""Stage 6: Multi-Isolate Conservation Audit Across 4 Target Genes, 50 Candidates, and 19 Isolates (BRD v3.1)."""

from __future__ import annotations

import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from Bio import AlignIO, SeqIO

try:
    from .paths_v3_1 import resolve_repo_root, resolve_v31_root
except ImportError:
    from paths_v3_1 import resolve_repo_root, resolve_v31_root


def run_stage6_conservation(repo_root: Path | None = None) -> Dict[str, Path]:
    root = resolve_repo_root(repo_root)
    v3 = root / "11_post_execution_audit_v3"
    v31 = resolve_v31_root(root)
    s6 = v31 / "stage6_conservation"
    aln_out = s6 / "conservation_alignment_files"
    s6.mkdir(parents=True, exist_ok=True)
    aln_out.mkdir(parents=True, exist_ok=True)
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    auth_t = pd.read_csv(v3 / "authoritative_data" / "auth_targets.tsv", sep="\t")
    auth_c = pd.read_csv(v3 / "authoritative_data" / "auth_candidates.tsv", sep="\t")
    auth_iso = pd.read_csv(v3 / "authoritative_data" / "auth_isolates.tsv", sep="\t")
    auth_g = pd.read_csv(v3 / "authoritative_data" / "auth_guides.tsv", sep="\t")
    cons_v3 = pd.read_csv(v3 / "stage3_isolate_conservation" / "candidate_conservation_multimetric_v3.tsv", sep="\t").set_index("candidate_id")
    pair_v3 = pd.read_csv(v3 / "stage3_isolate_conservation" / "candidate_by_isolate_alignments_v3.tsv", sep="\t")
    t_cons_v1 = pd.read_csv(root / "04_comparative_genomics" / "nucleotide_conservation" / "target_nucleotide_conservation.tsv", sep="\t").set_index("target")
    broc_t = pd.read_csv(v31 / "stage5_broccoli_orthology" / "broccoli_vs_diamond_target_genes_v3_1.tsv", sep="\t").set_index("target")

    # 1. Copy and directly audit multi-isolate CDS alignment files
    aln_src_dir = root / "04_comparative_genomics" / "target_orthologs" / "target_cds_alignments"
    gene_rows: List[Dict[str, Any]] = []
    for _, tr in auth_t.iterrows():
        tname = str(tr["target"])
        ref_cds = str(tr["cds_sequence"]).upper() if "cds_sequence" in tr else ""
        aln_fa = aln_src_dir / f"{tname}_cds_aligned.fasta"
        un_fa = aln_src_dir / f"{tname}_isolates_cds.fasta"
        if aln_fa.exists():
            shutil.copy2(aln_fa, aln_out / aln_fa.name)
        if un_fa.exists():
            shutil.copy2(un_fa, aln_out / un_fa.name)

        n_records = 0
        ambig_n_total = 0
        gap_cols_total = 0
        aln_len = 0
        seq_identities_vs_ref: List[float] = []
        coverages_vs_ref: List[float] = []

        if aln_fa.exists():
            aln = AlignIO.read(str(aln_fa), "fasta")
            n_records = len(aln)
            aln_len = aln.get_alignment_length()
            ref_rec = next((r for r in aln if "ND90Pr" in r.id), aln[0])
            ref_str = str(ref_rec.seq).upper()
            for r in aln:
                s = str(r.seq).upper()
                ambig_n_total += s.count("N")
                gaps = s.count("-")
                gap_cols_total += gaps
                comp_pairs = [(a, b) for a, b in zip(ref_str, s) if a != "-" and b != "-"]
                if comp_pairs:
                    same = sum(1 for a, b in comp_pairs if a == b)
                    seq_identities_vs_ref.append(same / len(comp_pairs) * 100.0)
                ref_ungapped = sum(1 for a in ref_str if a != "-")
                aligned_on_ref = sum(1 for a, b in zip(ref_str, s) if a != "-" and b != "-")
                coverages_vs_ref.append(aligned_on_ref / max(ref_ungapped, 1) * 100.0)

        tc = t_cons_v1.loc[tname]
        bt = broc_t.loc[tname]
        sub_p = pair_v3[pair_v3["target"] == tname]
        sub_hq = sub_p[sub_p["is_high_quality_non_artifact_isolate"].astype(bool)]

        gene_rows.append(
            {
                "target": tname,
                "locus_tag": str(tr["locus_tag"]),
                "refseq_protein_id": str(tr["refseq_protein_id"]),
                "ref_cds_length_nt": int(tr["cds_length_nt"]),
                "isolates_evaluated": int(tc["isolates_evaluated"]),
                "isolates_with_gene": int(tc["isolates_with_gene"]),
                "isolates_complete_cds": int(tc["isolates_complete_cds"]),
                "missing_or_incomplete_cds_isolates": int(tc["isolates_evaluated"]) - int(tc["isolates_complete_cds"]),
                "alignment_length_nt": aln_len,
                "msa_mean_nt_identity_vs_nd90pr_pct": round(sum(seq_identities_vs_ref) / max(len(seq_identities_vs_ref), 1), 4),
                "msa_min_nt_identity_vs_nd90pr_pct": round(min(seq_identities_vs_ref), 4) if seq_identities_vs_ref else 0.0,
                "msa_max_nt_identity_vs_nd90pr_pct": round(max(seq_identities_vs_ref), 4) if seq_identities_vs_ref else 0.0,
                "msa_mean_coverage_pct": round(sum(coverages_vs_ref) / max(len(coverages_vs_ref), 1), 4),
                "msa_min_coverage_pct": round(min(coverages_vs_ref), 4) if coverages_vs_ref else 0.0,
                "msa_gap_pct_across_matrix": round(gap_cols_total / max(n_records * aln_len, 1) * 100.0, 4),
                "msa_ambiguous_n_bases_total": int(ambig_n_total),
                "snp_sites_msa": int(tc["snp_sites"]),
                "indel_columns_msa": int(tc["indel_columns"]),
                "mean_21mer_isolate_conservation_pct": float(tc["mean_21mer_isolate_conservation_pct"]),
                "pct_100pct_conserved_21mers": float(tc["pct_100pct_conserved_21mers"]),
                "candidate_level_min_nt_id_hq17_pct": round(float(sub_hq["nucleotide_identity_pct"].min()), 4),
                "candidate_level_min_21mer_hq17_pct": round(float(sub_hq["exact_21mer_preservation_pct"].min()), 4),
                "broccoli_v2_orthogroup": str(bt["broccoli_og"]),
                "diamond_rbh_orthogroup": str(bt["diamond_rbh_og"]),
                "artifact_isolates_flagged": "WAI2406,WAI2411 (unpolished long-read homopolymer indels)",
            }
        )
    gene_tsv = s6 / "gene_conservation_v3_1.tsv"
    pd.DataFrame(gene_rows).to_csv(gene_tsv, sep="\t", index=False)

    # 2. Candidate conservation summary (explicitly prevents high mean hiding low-identity isolate)
    guide_groups = auth_g.groupby("candidate_id")
    cand_rows: List[Dict[str, Any]] = []
    for _, cr in auth_c.iterrows():
        cid = str(cr["candidate_id"])
        cv = cons_v3.loc[cid]
        sub_p = pair_v3[pair_v3["candidate_id"] == cid].copy()
        sub_hq = sub_p[sub_p["is_high_quality_non_artifact_isolate"].astype(bool)]

        min_all_row = sub_p.sort_values(["exact_21mer_preservation_pct", "nucleotide_identity_pct"], ascending=[True, True]).iloc[0]
        min_hq_row = sub_hq.sort_values(["exact_21mer_preservation_pct", "nucleotide_identity_pct"], ascending=[True, True]).iloc[0]
        max_nt_all = float(sub_p["nucleotide_identity_pct"].max())
        mean_cov_hq = float(sub_hq["alignment_coverage_pct"].mean())
        min_cov_hq = float(sub_hq["alignment_coverage_pct"].min())
        gap_bases_all = int(sub_p["gap_bases"].sum())
        gap_bases_hq = int(sub_hq["gap_bases"].sum())
        cand_len = int(cr["candidate_length_nt"])
        gap_pct_all = round(gap_bases_all / max(len(sub_p) * cand_len, 1) * 100.0, 4)
        gap_pct_hq = round(gap_bases_hq / max(len(sub_hq) * cand_len, 1) * 100.0, 4)

        gg = guide_groups.get_group(cid)
        n_windows = len(gg)
        n_windows_100_all19 = int((gg["isolates_conserved_100pct"].astype(int) == 19).sum()) if "isolates_conserved_100pct" in gg.columns else 0

        cand_rows.append(
            {
                "candidate_id": cid,
                "target": str(cr["target"]),
                "candidate_length_nt": cand_len,
                "total_21nt_windows": n_windows,
                "isolates_evaluated": int(cv["isolates_evaluated"]),
                "hq_non_artifact_isolates": int(cv["hq_non_artifact_isolates"]),
                "mean_nt_identity_all19_pct": float(cv["mean_nt_identity_all19_pct"]),
                "min_nt_identity_all19_pct": float(cv["min_nt_identity_all19_pct"]),
                "max_nt_identity_all19_pct": round(max_nt_all, 4),
                "mean_nt_identity_hq17_pct": float(cv["mean_nt_identity_hq17_pct"]),
                "min_nt_identity_hq17_pct": float(cv["min_nt_identity_hq17_pct"]),
                "lowest_identity_isolate_all19": f"{min_all_row['isolate_id']} ({min_all_row['nucleotide_identity_pct']}% nt / {min_all_row['exact_21mer_preservation_pct']}% 21-mer)",
                "lowest_identity_isolate_hq17": f"{min_hq_row['isolate_id']} ({min_hq_row['nucleotide_identity_pct']}% nt / {min_hq_row['exact_21mer_preservation_pct']}% 21-mer)",
                "mean_alignment_coverage_all19_pct": float(cv["mean_alignment_coverage_all19_pct"]),
                "min_alignment_coverage_all19_pct": float(cv["min_alignment_coverage_all19_pct"]),
                "mean_alignment_coverage_hq17_pct": round(mean_cov_hq, 4),
                "min_alignment_coverage_hq17_pct": round(min_cov_hq, 4),
                "gap_bases_all19": gap_bases_all,
                "gap_pct_all19": gap_pct_all,
                "gap_bases_hq17": gap_bases_hq,
                "gap_pct_hq17": gap_pct_hq,
                "missing_or_unresolved_isolates_n": int(cv["n_isolates_missing_or_unresolved"]),
                "ambiguous_n_bases_in_candidate": int(str(cr["sequence"]).upper().count("N")),
                "mean_21mer_preservation_all19_pct": float(cv["mean_21mer_preservation_all19_pct"]),
                "min_21mer_preservation_all19_pct": float(cv["min_21mer_preservation_all19_pct"]),
                "mean_21mer_preservation_hq17_pct": float(cv["mean_21mer_preservation_hq17_pct"]),
                "min_21mer_preservation_hq17_pct": float(cv["min_21mer_preservation_hq17_pct"]),
                "guides_pct_preserved_all19": float(cv["guides_pct_preserved_all19"]),
                "guides_pct_preserved_hq17": float(cv["guides_pct_preserved_hq17"]),
                "windows_21nt_100pct_all19_n": n_windows_100_all19,
                "conservation_rank_v3": int(cv["conservation_rank_v3"]),
            }
        )
    cand_tsv = s6 / "candidate_conservation_v3_1.tsv"
    pd.DataFrame(cand_rows).to_csv(cand_tsv, sep="\t", index=False)

    # 3. Wide per-isolate conservation matrix (50 candidates x 19 isolates for identity & 21-mer preservation)
    iso_order = auth_iso["isolate_id"].tolist()
    piv_id = pair_v3.pivot(index="candidate_id", columns="isolate_id", values="nucleotide_identity_pct")[iso_order]
    piv_id.columns = [f"{c}_nt_id_pct" for c in piv_id.columns]
    piv_k = pair_v3.pivot(index="candidate_id", columns="isolate_id", values="exact_21mer_preservation_pct")[iso_order]
    piv_k.columns = [f"{c}_21mer_pct" for c in piv_k.columns]
    mat_df = pd.concat([piv_id, piv_k], axis=1).reset_index()
    mat_tsv = s6 / "isolate_conservation_matrix_v3_1.tsv"
    mat_df.to_csv(mat_tsv, sep="\t", index=False)

    methods_md = f"""# Stage 6 — Multi-Isolate Conservation Methodology (BRD v3.1)

**Executed:** {now_utc}

## Audit Principles Enforced
1. **Per-gene and per-candidate separation:** Target-level MAFFT CDS alignments (`conservation_alignment_files/*_cds_aligned.fasta`) and candidate-level semi-global alignments (`candidate_by_isolate_alignments_v3.tsv`, `950` candidate-by-isolate pairs) are evaluated independently across all `19` *B. sorokiniana* isolates.
2. **Minimum identity & 21-mer window preservation:** Reporting enforces `min_nt_identity_hq17_pct`, `min_21mer_preservation_hq17_pct`, `lowest_identity_isolate_hq17`, and `lowest_identity_isolate_all19` alongside means so localized SNPs (e.g., `LK93` in `CsSp1`, dropping 21-mer preservation to `79.38%` in `BRD_CsSp1_CAND008` despite `98.89%` nucleotide identity) are never masked.
3. **Assembly polishing artifact separation:** `WAI2406` (`+137` excess internal stops in 50 core markers) and `WAI2411` (`+182` excess stops) are preserved in `ALL-19` summaries and cleanly separated from the `17` high-quality (`HQ-17`) assemblies.
4. **21-nt sliding window assessment:** Every possible 21-nt guide (`candidate_length_nt - 20`) within all `50` candidates (`9,301` guide windows total) is tracked across `ALL-19` and `HQ-17`.
"""
    (s6 / "conservation_methods.md").write_text(methods_md, encoding="utf-8")

    return {
        "gene_conservation": gene_tsv,
        "candidate_conservation": cand_tsv,
        "isolate_matrix": mat_tsv,
        "alignments_dir": aln_out,
        "methods": s6 / "conservation_methods.md",
    }


if __name__ == "__main__":
    for k, v in run_stage6_conservation().items():
        print(f"STAGE6_CONS_OUT {k}: {v}")
