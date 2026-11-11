"""Automated verification suite for Bipolaris-RNAi Designer v3 (BRD v3)."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pandas as pd
import pytest


def _v3_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _output_root() -> Path:
    return _v3_root().parent


def test_portability_no_hardcoded_machine_paths() -> None:
    v3 = _v3_root()
    out = _output_root()
    forbidden = ("ASUS " + "BANGLADESH", "New " + "folder (2)", r"C:\Us" + "ers")
    files_to_check = list((v3 / "src").glob("*.py")) + list((v3 / "tests").glob("*.py")) + list((out / "brd").rglob("*.py"))
    assert len(files_to_check) >= 8
    for p in files_to_check:
        txt = p.read_text(encoding="utf-8", errors="replace")
        for pat in forbidden:
            assert pat not in txt, f"Hardcoded machine path {pat!r} found in {p}"


def test_stage1_authoritative_invariants() -> None:
    v3 = _v3_root()
    out = _output_root()
    auth = v3 / "authoritative_data"
    targets = pd.read_csv(auth / "auth_targets.tsv", sep="\t")
    isolates = pd.read_csv(auth / "auth_isolates.tsv", sep="\t")
    cands = pd.read_csv(auth / "auth_candidates.tsv", sep="\t")
    guides = pd.read_csv(auth / "auth_guides.tsv", sep="\t")

    assert set(targets["target"]) == {"CsFUS3", "CsSTE11", "CsSp1", "CsPKS1"}
    assert len(isolates) == 19
    for rel in isolates["relative_genome_path"]:
        assert (out / str(rel)).exists(), f"Missing relative genome path: {rel}"
    artifact_isos = set(isolates[isolates["assembly_quality_flag"].str.contains("HOMOPOLYMER|ARTIFACT", case=False, na=False)]["isolate_id"])
    assert {"WAI2406", "WAI2411"}.issubset(artifact_isos)

    assert len(cands) == 50
    for _, r in cands.iterrows():
        seq = str(r["sequence"]).upper()
        assert len(seq) == int(r["candidate_length_nt"])
        assert hashlib.sha256(seq.encode("utf-8")).hexdigest() == str(r["sequence_sha256"])
        gc = round((seq.count("G") + seq.count("C")) / len(seq) * 100.0, 2)
        assert abs(gc - float(r["gc_percent"])) <= 0.05
        assert int(r["cds_end_nt"]) > int(r["cds_start_nt"])
        assert int(r["genomic_end_1based"]) >= int(r["genomic_start_1based"])

    cand_seq_map = cands.set_index("candidate_id")["sequence"].to_dict()
    assert len(guides) > 500
    for cid, g_seq in zip(guides["candidate_id"], guides["sense_dna_21nt"]):
        assert len(str(g_seq)) == 21
        assert str(g_seq) in str(cand_seq_map[cid]).upper()

    rec = pd.read_csv(v3 / "stage1_reconciliation" / "reconciliation_log_v3.tsv", sep="\t")
    assert len(rec) >= 4
    joined = " ".join(rec.astype(str).to_numpy().ravel())
    assert "1438.44" in joined and "0.9170" in joined


def test_stage2_honest_blockers_and_repaired_rbh() -> None:
    v3 = _v3_root()
    s2 = v3 / "stage2_genome_qc_comparative"
    bl = pd.read_csv(s2 / "unavailable_failed_or_incomplete_analyses_v3.tsv", sep="\t").set_index("tool")
    for tool in ("QUAST", "BUSCO", "OrthoFinder"):
        assert tool in bl.index
        assert bl.loc[tool, "status"] == "NOT_RUN_DEPENDENCY_BLOCKED"

    # Must NEVER invent QUAST/BUSCO/OrthoFinder outputs
    assert not list(s2.glob("*quast_report*"))
    assert not list(s2.glob("*busco_summary*"))

    core = pd.read_csv(s2 / "orthology_diamond_rbh" / "core_accessory_summary_v3.tsv", sep="\t")
    n_core_both = int((core["bsor_annotated_present"] == 2).sum())
    assert n_core_both > 5000, f"Expected repaired core orthogroups >5000, got {n_core_both}"

    iso_screen = pd.read_csv(s2 / "core_integrity_screen" / "isolates_artifact_vs_quality_v3.tsv", sep="\t")
    assert len(iso_screen) == 19
    low_q = set(iso_screen[iso_screen["total_internal_stops_across_50_markers"] >= 100]["isolate_id"])
    assert low_q == {"WAI2406", "WAI2411"}


def test_stage3_conservation_math_and_sensitivity() -> None:
    s3 = _v3_root() / "stage3_isolate_conservation"
    cm = pd.read_csv(s3 / "candidate_conservation_multimetric_v3.tsv", sep="\t")
    assert len(cm) == 50
    for _, r in cm.iterrows():
        assert r["min_nt_identity_hq17_pct"] <= r["mean_nt_identity_hq17_pct"] + 1e-6
        assert r["min_21mer_preservation_all19_pct"] <= r["mean_21mer_preservation_all19_pct"] + 1e-6
        assert int(r["hq_non_artifact_isolates"]) == 17  # WAI2406 & WAI2411 excluded from HQ-17
        assert 0 <= int(r["n_isolates_missing_or_unresolved"]) <= 2

    sens = pd.read_csv(s3 / "conservation_threshold_sensitivity_v3.tsv", sep="\t")
    assert len(sens) >= 9
    assert set(sens["threshold_pct"]).issuperset({85.0, 90.0, 95.0})


def test_stage4_expanded_panel_and_selectivity_invariants() -> None:
    s4 = _v3_root() / "stage4_off_target_expanded"
    manifest = pd.read_csv(s4 / "expanded_reference_manifest_v3.tsv", sep="\t")
    assert len(manifest) >= 15
    assert {"crop_wheat", "crop_barley", "crop_rice", "crop_maize", "beneficial_fungi", "bacteria"}.issubset(set(manifest["group"]))

    gaps = pd.read_csv(s4 / "reference_panel_coverage_and_gaps_v3.tsv", sep="\t")
    assert len(gaps) >= 17
    assert "screening_scope_note" in gaps.columns
    assert int((gaps["transcript_screened"].astype(bool) | gaps["genome_screened"].astype(bool)).sum()) >= 17

    sel = pd.read_csv(s4 / "candidate_selectivity_expanded_v3.tsv", sep="\t")
    assert len(sel) == 50
    assert "heuristic_sequence_selectivity_index" in sel.columns
    assert "safety_score" not in " ".join(sel.columns).lower()
    for _, r in sel.iterrows():
        any_exact = (
            int(r["wheat_exact_21"]) + int(r["barley_exact_21"]) + int(r["rice_exact_21"]) + int(r["maize_exact_21"])
            + int(r["beneficial_fungi_tx_exact_21"]) + int(r["beneficial_fungi_genome_exact_21"])
            + int(r["bacteria_exact_21"]) + int(r["self_paralog_or_genome_exact_21"])
        ) > 0
        assert bool(r["passes_zero_exact_21mer_filter"]) is (not any_exact)

    # Unit check dual-strand 21-mer matching helper + exact equivalence to naive sliding window across all offsets
    import sys as _sys
    _sys.path.insert(0, str(_v3_root() / "src"))
    from stage4_off_target_expanded import _scan_exact_21mers  # noqa: PLC0415

    tmp_fa = _v3_root() / "_tmp_unit_test_strand.fa"
    try:
        guide = "ATCGATCGATCGATCGATCGA"
        # Place matches at EVERY offset mod 11 (0..10) + end boundary on both strands
        seq1 = "".join("G" * off + guide + "T" * 7 for off in range(12))
        rc_guide = "TCGATCGATCGATCGATCGAT"
        seq2 = "".join("C" * off + rc_guide + "A" * 5 for off in range(12))
        tmp_fa.write_text(f">tx1\n{seq1}\n>tx2\n{seq2}\n", encoding="utf-8")
        hits = _scan_exact_21mers(tmp_fa, {guide})
        strands = {h[2] for h in hits}
        assert strands == {"+", "-"}
        # Compare strictly to brute-force naive scan
        naive = set()
        for rid, sq in (("tx1", seq1), ("tx2", seq2)):
            for i in range(len(sq) - 20):
                sub = sq[i : i + 21]
                if sub == guide:
                    naive.add((rid, i, "+"))
                elif sub == rc_guide:
                    naive.add((rid, i, "-"))
        assert set(hits) == naive
    finally:
        if tmp_fa.exists():
            tmp_fa.unlink()


def test_stage5_clustering_and_non_overlap() -> None:
    s5 = _v3_root() / "stage5_redundancy_design"
    ov = pd.read_csv(s5 / "candidate_pairwise_overlaps_v3.tsv", sep="\t")
    cl = pd.read_csv(s5 / "candidate_to_region_cluster_map_v3.tsv", sep="\t").set_index("candidate_id")
    reg = pd.read_csv(s5 / "independent_target_regions_v3.tsv", sep="\t")

    assert len(cl) == 50
    assert 1 <= len(reg) < 50
    for _, r in ov.iterrows():
        if bool(r["clustered_as_same_region"]):
            assert cl.loc[r["candidate_a"], "region_id"] == cl.loc[r["candidate_b"], "region_id"]
    sweep = pd.read_csv(s5 / "exhaustive_cds_180nt_sweep_v3.tsv", sep="\t")
    assert set(sweep["target"]) == {"CsFUS3", "CsSTE11", "CsSp1", "CsPKS1"}


def test_stage6_pareto_separation_and_reproducibility() -> None:
    s6 = _v3_root() / "stage6_multidim_ranking"
    bio = pd.read_csv(s6 / "biological_target_relevance_v3.tsv", sep="\t").set_index("target")
    # PHI:2260 unaffected_pathogenicity must prevent CsPKS1 from outranking CsFUS3/CsSTE11!
    assert bio.loc["CsFUS3", "biological_target_relevance_score"] > bio.loc["CsPKS1", "biological_target_relevance_score"]
    assert bio.loc["CsSTE11", "biological_target_relevance_score"] > bio.loc["CsPKS1", "biological_target_relevance_score"]

    rk = pd.read_csv(s6 / "candidate_multidim_ranking_v3.tsv", sep="\t")
    assert len(rk) == 50
    for _, r in rk.iterrows():
        recalc = round(
            0.30 * float(r["biological_target_relevance_score"])
            + 0.35 * float(r["dsrna_targetability_score"])
            + 0.35 * float(r["sequence_selectivity_score"]),
            4,
        )
        assert abs(recalc - float(r["heuristic_composite_prioritization_index"])) <= 1e-4

    short = pd.read_csv(s6 / "independent_region_shortlist_v3.tsv", sep="\t")
    assert short["region_id"].is_unique
    assert len(short) == int(pd.read_csv(_v3_root() / "stage5_redundancy_design" / "independent_target_regions_v3.tsv", sep="\t").shape[0])

    # Unit check Pareto solver
    import sys as _sys
    _sys.path.insert(0, str(_v3_root() / "src"))
    from stage6_multidim_ranking import _compute_pareto_tiers  # noqa: PLC0415

    df_u = pd.DataFrame(
        [
            {"a": 0.9, "b": 0.5, "c": 0.9},
            {"a": 0.5, "b": 0.9, "c": 0.9},
            {"a": 0.4, "b": 0.4, "c": 0.4},
        ]
    )
    assert _compute_pareto_tiers(df_u, ["a", "b", "c"]) == [1, 1, 2]

    sens = pd.read_csv(s6 / "stage6_weight_sensitivity_v3.tsv", sep="\t")
    assert sens["scenario"].nunique() == 5
    boot = pd.read_csv(s6 / "stage6_bootstrap_stability_v3.tsv", sep="\t")
    assert len(boot) >= 1


def test_stage7_literature_benchmark_categories() -> None:
    s7 = _v3_root() / "stage7_literature_benchmark"
    ds = pd.read_csv(s7 / "literature_benchmark_dataset_v3.tsv", sep="\t")
    assert len(ds) >= 10
    assert ds["doi"].notna().all()
    allowed_cats = {
        "1_Direct_dsRNA_Sequence",
        "2_Reconstructed_Sequence",
        "3_Gene_Level_Only",
        "4_Computational_Only",
        "5_Ambiguous_or_Conflicting_Outcome",
    }
    assert set(ds["benchmark_category"]).issubset(allowed_cats)
    assert "3_Gene_Level_Only" in set(ds["benchmark_category"])
    assert "5_Ambiguous_or_Conflicting_Outcome" in set(ds["benchmark_category"])


def test_stage8_report_18_sections_and_hashes() -> None:
    v3 = _v3_root()
    rep = v3 / "V3_Validation_and_Execution_Report.md"
    assert rep.exists()
    txt = rep.read_text(encoding="utf-8")
    for sec_num in range(1, 19):
        assert re.search(rf"^##\s+{sec_num}\.", txt, re.MULTILINE), f"Missing mandatory Section {sec_num} in V3 report"

    # Must NOT assert unproven biological efficacy or safety without explicit negation
    cleaned_txt = txt.lower().replace("**", "")
    assert not re.search(r"(?<!not )proven fungicide", cleaned_txt)
    assert "not proven fungicides" in cleaned_txt
    assert any(ph in cleaned_txt for ph in ("environmental biosafety", "ecological safety", "biological safety"))
    hashes = pd.read_csv(v3 / "hashes_v3.tsv", sep="\t")
    assert len(hashes) >= 15
    for _, hr in hashes.iterrows():
        fp = _output_root() / str(hr["relative_path"])
        assert fp.exists(), f"Hashed file missing: {fp}"
        assert hashlib.sha256(fp.read_bytes()).hexdigest() == str(hr["sha256"])
