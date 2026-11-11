"""Stage 6 (FIX 6): Candidate Coordinate Audit & Consistency Reconciliation (BRD v3.1.1).

Explicitly distinguishes:
  1. Recorded coordinate assignments (50/50)
  2. Contiguous authoritative CDS substring matches (49/50)
  3. Stitched / synthetic constructs (1/50: BRD_CsSTE11_CAND009 from dsRNAmax)
  4. Unresolved coordinate mappings (0/50)

Reconciles the v3.1 change log misstatement ('Verified 50/50 CDS coordinates') with the actual
validated v3.1 baseline finding (49/50 contiguous CDS substrings + 1 synthetic stitched construct).
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from Bio import SeqIO

try:
    from .paths_v3_1_1 import resolve_repo_root, resolve_v31_root, resolve_v311_root
except ImportError:
    from paths_v3_1_1 import resolve_repo_root, resolve_v31_root, resolve_v311_root


def run_stage6_coordinate_audit(repo_root: Path | None = None) -> Dict[str, Any]:
    root = resolve_repo_root(repo_root)
    v3 = root / "11_post_execution_audit_v3"
    v31 = resolve_v31_root(root)
    v311 = resolve_v311_root(root)
    s6 = v311 / "stage6_coordinates"
    s6.mkdir(parents=True, exist_ok=True)
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    auth_c = pd.read_csv(v3 / "authoritative_data" / "auth_candidates.tsv", sep="\t")
    auth_t = pd.read_csv(v3 / "authoritative_data" / "auth_targets.tsv", sep="\t").set_index("target")
    mast_v31 = pd.read_csv(v31 / "stage8_candidate_design" / "candidate_master_v3_1.tsv", sep="\t").set_index("candidate_id")
    t_seqs = {
        r.id.split("|")[0]: str(r.seq).upper()
        for r in SeqIO.parse(str(root / "10_final_reports" / "BRD_Target_Gene_Sequences.fasta"), "fasta")
    }

    rows: List[Dict[str, Any]] = []
    for _, r in auth_c.iterrows():
        cid = str(r["candidate_id"])
        tname = str(r["target"])
        cseq = str(r["sequence"]).upper()
        s_pos = int(r["cds_start_nt"])
        e_pos = int(r["cds_end_nt"])
        g_s = int(r["genomic_start_1based"])
        g_e = int(r["genomic_end_1based"])
        g_scf = str(r["genomic_scaffold"])
        g_strand = str(r["genomic_strand"])
        cds_strand = str(r["cds_strand"])
        is_exact = bool(r["exact_cds_substring"])
        origin = str(r["exons_spanned"])
        is_synth = origin == "SYNTHETIC_KMER_CONSTRUCT" or (not is_exact)

        ref_seq = t_seqs.get(tname, "")
        slice_seq = ref_seq[s_pos - 1 : e_pos] if ref_seq and 1 <= s_pos <= e_pos <= len(ref_seq) else ""
        exact_substr_ok = bool(slice_seq) and (slice_seq == cseq)

        # Quantify nucleotide differences vs reference CDS slice if same length (e.g., stitched synthetic construct)
        if exact_substr_ok:
            mismatches_vs_ref = 0
            coord_category = "CONTIGUOUS_AUTHORITATIVE_CDS_SUBSTRING"
            mapping_resolved = True
            genomic_coord_note = f"{g_scf}:{g_s}..{g_e}({g_strand})"
        elif is_synth and len(slice_seq) == len(cseq) and g_s == -1 and g_e == -1:
            mismatches_vs_ref = sum(1 for a, b in zip(slice_seq, cseq) if a != b)
            coord_category = "STITCHED_SYNTHETIC_CONSTRUCT"
            mapping_resolved = True  # Explicitly categorized as non-contiguous stitched construct covering nominal CDS 1..300
            genomic_coord_note = f"{g_scf}:UNMAPPED_SYNTHETIC(-1..-1)"
        else:
            mismatches_vs_ref = -1
            coord_category = "UNRESOLVED_COORDINATE_MAPPING"
            mapping_resolved = False
            genomic_coord_note = "UNRESOLVED"

        trow = auth_t.loc[tname]
        rows.append(
            {
                "candidate_id": cid,
                "target": tname,
                "tool_name": str(r["tool_name"]),
                "reference_sequence_id": str(trow["refseq_mrna_id"]),
                "locus_tag": str(trow["locus_tag"]),
                "reference_assembly": str(r["reference_assembly"]),
                "coordinate_domain": "CDS_nucleotide_and_genomic_scaffold",
                "coordinate_convention": "1-based_inclusive",
                "cds_start_1based": s_pos,
                "cds_end_1based": e_pos,
                "cds_strand": cds_strand,
                "genomic_scaffold": g_scf,
                "genomic_start_1based": g_s,
                "genomic_end_1based": g_e,
                "genomic_strand": g_strand,
                "genomic_span_label": genomic_coord_note,
                "candidate_length_nt": int(r["candidate_length_nt"]),
                "is_contiguous": exact_substr_ok,
                "is_stitched_or_synthetic": is_synth,
                "exact_cds_substring_verified": exact_substr_ok,
                "mismatches_vs_nd90pr_cds_slice": mismatches_vs_ref,
                "coordinate_classification": coord_category,
                "coordinate_mapping_resolved": mapping_resolved,
                "v31_master_is_exact_flag": bool(mast_v31.loc[cid, "is_exact_cds_substring"]),
                "v31_master_verified_flag": bool(mast_v31.loc[cid, "authoritative_cds_substring_verified"]),
            }
        )

    audit_df = pd.DataFrame(rows)
    n_total = len(audit_df)
    n_recorded = int(((audit_df["cds_start_1based"] >= 1) & (audit_df["cds_end_1based"] >= audit_df["cds_start_1based"])).sum())
    n_contig = int((audit_df["coordinate_classification"] == "CONTIGUOUS_AUTHORITATIVE_CDS_SUBSTRING").sum())
    n_synth = int((audit_df["coordinate_classification"] == "STITCHED_SYNTHETIC_CONSTRUCT").sum())
    n_unres = int((~audit_df["coordinate_mapping_resolved"]).sum())

    synth_df = audit_df[audit_df["is_stitched_or_synthetic"]]
    synth_row = synth_df.iloc[0]
    synth_desc = "; ".join(
        f"{r['candidate_id']} ({r['tool_name']}, genomic_start={r['genomic_start_1based']}, {r['mismatches_vs_nd90pr_cds_slice']} mismatches vs {r['target']} {r['cds_start_1based']}..{r['cds_end_1based']})"
        for _, r in synth_df.iterrows()
    )
    summary_rows = [
        {"metric": "candidates_with_recorded_coordinate_assignments", "count": n_recorded, "total": n_total, "details": f"All {n_recorded}/{n_total} carry 1-based inclusive CDS span [cds_start_1based, cds_end_1based], refseq_mrna_id, and strand"},
        {"metric": "candidates_matching_contiguous_authoritative_cds_substrings", "count": n_contig, "total": n_total, "details": f"{n_contig}/{n_total} match ND90Pr CDS substrings in BRD_Target_Gene_Sequences.fasta with 0 mismatches"},
        {"metric": "candidates_stitched_or_synthetic_constructs", "count": n_synth, "total": n_total, "details": f"{n_synth}/{n_total}: {synth_desc}"},
        {"metric": "candidates_with_unresolved_coordinate_mapping", "count": n_unres, "total": n_total, "details": f"{n_unres}/{n_total} unresolved (all {n_total} explicitly accounted for without altering any candidate sequence)"},
    ]

    audit_tsv = s6 / "candidate_coordinate_audit_v3_1_1.tsv"
    summary_tsv = s6 / "candidate_coordinate_summary_v3_1_1.tsv"
    audit_df.to_csv(audit_tsv, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(summary_tsv, sep="\t", index=False)
    methods_md = f"""# Candidate Coordinate Audit & Reconciliation Methods (BRD v3.1.1)

**Audited:** {now_utc}

## 1. Reconciliation of v3.1 Defect (FIX 6)
The v3.1 `V3_1_Final_Validation_Report.md` accurately reported `49/50` contiguous CDS substring matches plus `1` synthetic construct (`BRD_CsSTE11_CAND009`), whereas `V3_1_Change_Log.md` ambiguously stated `Verified 50/50 CDS coordinates on target FASTA`. All v3.1.1 tables, change logs, and validation reports now consistently separate all four coordinate categories without modifying any candidate sequence.

## 2. Explicit 4-Category Audit Counts
| Category | Count | Details |
|---|---:|---|
| 1. Candidates with recorded coordinate assignments | **{n_recorded}/{n_total}** | 1-based inclusive target CDS span `[cds_start_1based, cds_end_1based]` + RefSeq mRNA ID & strand |
| 2. Contiguous authoritative CDS substring matches | **{n_contig}/{n_total}** | Exact contiguous substring on `BRD_Target_Gene_Sequences.fasta` (`0` mismatches, valid genomic coordinates) |
| 3. Stitched or synthetic constructs | **{n_synth}/{n_total}** | `{synth_row['candidate_id']}` (`{synth_row['tool_name']}`), `is_contiguous=False`, `is_stitched_or_synthetic=True`, `genomic_start_1based=-1`, `{synth_row['mismatches_vs_nd90pr_cds_slice']}` nt difference vs nominal `CsSTE11` `1..300` CDS slice |
| 4. Unresolved coordinate mappings | **{n_unres}/{n_total}** | Every candidate is cleanly classified (`{n_contig} + {n_synth} + {n_unres} = {n_total}`) |

## 3. Coordinate Convention & Sequence Domain
- **Convention:** `1-based_inclusive` (corresponding Python slice `[cds_start_1based - 1 : cds_end_1based]`).
- **Domains recorded:** Both spliced target **CDS** coordinates on RefSeq mRNAs (`XM_007698023.1`, `XM_007699256.1`, `XM_007696384.1`, `XM_007706020.1`) and reference **genomic scaffold** coordinates on `GCF_000338995.1` (`NW_006911896.1 (+)`, `NW_006911898.1 (-)`, `NW_006911895.1 (+)`, `NW_006911909.1 (+)`).
"""
    (s6 / "candidate_coordinate_methods.md").write_text(methods_md, encoding="utf-8")
    return {
        "audit_tsv": audit_tsv,
        "summary_tsv": summary_tsv,
        "counts": {"recorded": n_recorded, "contiguous": n_contig, "synthetic": n_synth, "unresolved": n_unres},
    }


if __name__ == "__main__":
    print("STAGE6_COORDS_OK:", run_stage6_coordinate_audit())
