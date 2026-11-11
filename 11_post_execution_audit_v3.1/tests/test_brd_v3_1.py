"""Comprehensive verification suite for Bipolaris-RNAi Designer v3.1 (clean-extraction portable)."""

from __future__ import annotations

import hashlib
import importlib
import sys
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO


def _v31_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _repo_root() -> Path:
    return _v31_root().parent


sys.path.insert(0, str(_v31_root() / "src"))

from assembly_qc_biopython import analyze_fasta_assembly, compute_n50_l50  # noqa: E402
from stage4_compleasm_runner import parse_compleasm_full_table, parse_compleasm_summary  # noqa: E402
from stage5_broccoli_orthology import parse_broccoli_ortholog_pairs, parse_broccoli_orthologous_groups  # noqa: E402


def test_01_portability_no_hardcoded_user_paths() -> None:
    v31 = _v31_root()
    forbidden = ("ASUS " + "BANGLADESH", "New " + "folder (2)", r"C:\Us" + "ers")
    py_files = list((v31 / "src").glob("*.py")) + list((v31 / "tests").glob("*.py"))
    assert len(py_files) >= 7
    for fp in py_files:
        txt = fp.read_text(encoding="utf-8", errors="replace")
        for pat in forbidden:
            assert pat not in txt, f"Hardcoded path fragment {pat!r} found in {fp}"


def test_02_optional_viennarna_import_guard() -> None:
    """Optional ViennaRNA import must not break non-folding BRD integration imports."""
    for mod_name in ("brd.integrations.rnahybrid", "brd.integrations.siren", "brd.integrations.lucid"):
        mod = importlib.import_module(mod_name)
        assert mod is not None


def test_03_tool_investigation_and_routes() -> None:
    v31 = _v31_root()
    ti = v31 / "tool_investigation"
    ee = v31 / "execution_environment"
    comp = pd.read_csv(ti / "github_tool_comparison.tsv", sep="\t")
    assert len(comp) == 3
    repos = set(comp["Repository"])
    assert repos == {
        "https://github.com/biopython/biopython",
        "https://github.com/huangnengCSU/compleasm",
        "https://github.com/rderelle/Broccoli",
    }
    for req_col in ("Version", "License", "Operating_system", "Dependencies", "Installation", "Execution", "Input", "Output", "Status"):
        assert req_col in comp.columns
        assert comp[req_col].notna().all()

    routes = pd.read_csv(ee / "execution_routes_status.tsv", sep="\t")
    assert len(routes) == 3
    assert "NOT_EXECUTED" in " ".join(routes["availability"]) or "UNAVAILABLE" in " ".join(routes["availability"])


def test_04_n50_l50_synthetic_and_fasta_qc_edge_cases(tmp_path: Path) -> None:
    # Known deterministic N50/L50 check: lengths [50, 30, 20] -> total 100 -> 50% = 50 -> N50=50, L50=1
    assert compute_n50_l50([20, 50, 30]) == (50, 1)
    # [40, 30, 20, 10] -> total 100 -> 50% = 50 -> 40+30=70 >= 50 -> N50=30, L50=2
    assert compute_n50_l50([10, 40, 20, 30]) == (30, 2)
    assert compute_n50_l50([]) == (0, 0)

    # Empty FASTA warning
    empty_fa = tmp_path / "empty.fna"
    empty_fa.write_text("", encoding="utf-8")
    summ_e, _, warns_e = analyze_fasta_assembly(empty_fa, isolate_id="EMPTY_TEST")
    assert summ_e["status"] == "EMPTY"
    assert any(w["warning_type"] == "EMPTY_FASTA" for w in warns_e)

    # Duplicate header + duplicate sequence + invalid char + short assembly warnings
    bad_fa = tmp_path / "bad.fna"
    bad_fa.write_text(">c1\nATGCATGC\n>c1\nATGCATGC\n>c3\nATGX\n", encoding="utf-8")
    summ_b, _, warns_b = analyze_fasta_assembly(bad_fa, isolate_id="BAD_TEST", min_expected_bp=100)
    w_types = {w["warning_type"] for w in warns_b}
    assert {"DUPLICATE_SEQUENCE_IDS", "DUPLICATE_CONTIG_SEQUENCES", "INVALID_NUCLEOTIDE_CHARACTERS", "SUSPICIOUSLY_SHORT_ASSEMBLY"}.issubset(w_types)
    assert summ_b["n50_bp"] == 8
    assert summ_b["l50"] == 2
    assert summ_b["invalid_char_count"] == 1

    # Missing assembly warning
    _, _, warns_m = analyze_fasta_assembly(tmp_path / "nonexistent.fna", isolate_id="MISS")
    assert warns_m[0]["warning_type"] == "MISSING_ASSEMBLY"


def test_05_stage3_biopython_assembly_qc_19_isolates() -> None:
    s3 = _v31_root() / "stage3_assembly_qc"
    df = pd.read_csv(s3 / "assembly_qc_summary_v3_1.tsv", sep="\t")
    assert len(df) == 19
    assert (df["status"] == "OK").all()
    assert (df["duplicate_id_count"] == 0).all()
    assert (df["invalid_char_count"] == 0).all()
    assert "Not QUAST" in (s3 / "assembly_qc_methods.md").read_text(encoding="utf-8")
    # Never claim QUAST executed
    assert not any("QUAST_CLI" in str(x) for x in df["qc_engine"])


def test_06_stage4_compleasm_parsers_and_honest_blocker(tmp_path: Path) -> None:
    syn_sum = tmp_path / "summary.txt"
    syn_sum.write_text(
        "## lineage: pleosporales_odb12\nS:98.12%, 5952\nD:0.41%, 25\nF:0.35%, 21\nI:0.12%, 7\nM:1.00%, 61\nN:6066\n",
        encoding="utf-8",
    )
    parsed = parse_compleasm_summary(syn_sum, "pleosporales_odb12")
    assert parsed["total_markers_assessed"] == 6066
    assert parsed["complete_markers"] == 5977
    assert parsed["fragmented_markers_total"] == 28
    assert parsed["missing_markers"] == 61
    assert parsed["completeness_pct"] == 98.53

    syn_ft = tmp_path / "full_table.tsv"
    syn_ft.write_text("Gene\tStatus\tSequence\nM1\tSingle\tctg1\nM2\tMissing\t.\n", encoding="utf-8")
    ft_df = parse_compleasm_full_table(syn_ft)
    assert len(ft_df) == 2

    s4 = _v31_root() / "stage4_compleasm"
    csum = pd.read_csv(s4 / "compleasm_summary_v3_1.tsv", sep="\t")
    assert len(csum) == 19
    assert (csum["status"] == "NOT_RUN_DEPENDENCY_BLOCKED").all()
    lin = pd.read_csv(s4 / "compleasm_lineage_metadata.tsv", sep="\t")
    assert "pleosporales_odb12" in set(lin["lineage"])
    assert int(lin[lin["lineage"] == "pleosporales_odb12"]["marker_count"].iloc[0]) == 6066


def test_07_stage5_broccoli_v2_outputs_and_rbh_comparison(tmp_path: Path) -> None:
    # Synthetic parser validation
    syn_og = tmp_path / "table_OGs_protein_names.txt"
    syn_og.write_text("#OG_name\tspA.faa\tspB.faa\nOG_1\tA1\tB1 B2\n", encoding="utf-8")
    og_p = parse_broccoli_orthologous_groups(syn_og, ["spA", "spB"])
    assert len(og_p) == 1

    syn_pr = tmp_path / "orthologous_pairs.txt"
    syn_pr.write_text("#protA\tprotB\nA1\tB1\n", encoding="utf-8")
    pr_p = parse_broccoli_ortholog_pairs(syn_pr)
    assert len(pr_p) == 1

    s5 = _v31_root() / "stage5_broccoli_orthology"
    comp = pd.read_csv(s5 / "broccoli_vs_diamond_rbh_comparison_v3_1.tsv", sep="\t")
    targets = pd.read_csv(s5 / "broccoli_vs_diamond_target_genes_v3_1.tsv", sep="\t").set_index("target")
    assert set(targets.index) == {"CsSTE11", "CsPKS1", "CsFUS3", "CsSp1"}
    assert targets.loc["CsSTE11", "broccoli_og"] == "OG_1786"
    assert targets.loc["CsPKS1", "broccoli_og"] == "OG_4251"
    assert targets.loc["CsFUS3", "broccoli_og"] == "OG_1315"
    assert targets.loc["CsSp1", "broccoli_og"] == "OG_9092"
    assert int(comp.iloc[0]["analysis_b_broccoli_v2_nj"]) == 9589


def test_08_stage6_conservation_invariants() -> None:
    s6 = _v31_root() / "stage6_conservation"
    gc = pd.read_csv(s6 / "gene_conservation_v3_1.tsv", sep="\t")
    cc = pd.read_csv(s6 / "candidate_conservation_v3_1.tsv", sep="\t")
    mat = pd.read_csv(s6 / "isolate_conservation_matrix_v3_1.tsv", sep="\t")
    assert len(gc) == 4
    assert len(cc) == 50
    assert len(mat) == 50
    for _, r in cc.iterrows():
        assert r["min_nt_identity_hq17_pct"] <= r["mean_nt_identity_hq17_pct"] + 1e-6 <= r["max_nt_identity_all19_pct"] + 1e-6
        assert r["min_21mer_preservation_hq17_pct"] <= r["mean_21mer_preservation_hq17_pct"] + 1e-6
    for tname in ("CsFUS3", "CsPKS1", "CsSTE11", "CsSp1"):
        assert (s6 / "conservation_alignment_files" / f"{tname}_cds_aligned.fasta").exists()


def test_09_stage7_offtarget_checksums_and_near_positions() -> None:
    s7 = _v31_root() / "stage7_offtarget"
    meta = pd.read_csv(s7 / "offtarget_reference_metadata.tsv", sep="\t")
    ex = pd.read_csv(s7 / "offtarget_exact_hits_v3_1.tsv", sep="\t")
    near = pd.read_csv(s7 / "offtarget_near_hits_v3_1.tsv", sep="\t")
    summ = pd.read_csv(s7 / "offtarget_summary_v3_1.tsv", sep="\t")

    assert len(meta) == 18
    assert meta["transcript_fasta_sha256"].str.len().ge(64).all()
    assert len(ex) == 1214
    assert "plus(+)" in set(ex["alignment_orientation"]) or "minus(-)" in set(ex["alignment_orientation"])
    assert len(near) == 16073
    assert (near["total_21nt_mismatches"] >= 1).all()
    assert (near["mismatch_positions_1based"].astype(str) != "NONE").all()
    assert len(summ) == 50
    assert int(summ["passes_zero_exact_21mer_filter"].sum()) == 32


def test_10_stage8_to_11_coordinates_targets_ranking_and_literature() -> None:
    v31 = _v31_root()
    root = _repo_root()
    mast = pd.read_csv(v31 / "stage8_candidate_design" / "candidate_master_v3_1.tsv", sep="\t")
    fa_seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(str(v31 / "stage8_candidate_design" / "candidate_sequences_v3_1.fasta"), "fasta")}
    assert len(mast) == 50
    exact_mask = mast["is_exact_cds_substring"].astype(bool)
    assert int(exact_mask.sum()) == 49
    assert (mast.loc[exact_mask, "authoritative_cds_substring_verified"].astype(bool)).all()
    assert list(mast.loc[~exact_mask, "candidate_id"]) == ["BRD_CsSTE11_CAND009"]
    assert (mast["length_consistency_verified"].astype(bool)).all()
    for _, r in mast.iterrows():
        assert fa_seqs[r["candidate_id"]] == str(r["sequence"]).upper()

    tev = pd.read_csv(v31 / "stage9_biological_targets" / "target_evidence_table_v3_1.tsv", sep="\t").set_index("gene_id")
    assert "PHI:2260" in str(tev.loc["CsPKS1", "external_database_evidence"])
    assert "139.61" in str(tev.loc["CsSp1", "expression_evidence_separate"])

    rk = pd.read_csv(v31 / "stage10_ranking" / "candidate_ranking_v3_1.tsv", sep="\t")
    reg_rk = pd.read_csv(v31 / "stage10_ranking" / "independent_region_ranking_v3_1.tsv", sep="\t")
    assert len(rk) == 50
    assert len(reg_rk) == 12
    for _, r in rk.iterrows():
        recalc = round(
            0.30 * float(r["biological_target_relevance_score"])
            + 0.35 * float(r["dsrna_targetability_score"])
            + 0.35 * float(r["sequence_selectivity_score"]),
            4,
        )
        assert abs(recalc - float(r["heuristic_composite_prioritization_index"])) <= 1e-4

    # Near-match burden visible prominently for shortlist leaders
    clean_reg = reg_rk[reg_rk["passes_zero_exact_21mer_filter"]]
    assert (clean_reg["near_20of21_guides_total"] >= 3).all()
    assert str(clean_reg.iloc[0]["candidate_id"]) == "BRD_CsSTE11_CAND013"

    lit = pd.read_csv(v31 / "stage11_literature_benchmark" / "literature_benchmark_v3_1.tsv", sep="\t")
    assert len(lit) == 25
    assert "relevance_to_current_targets" in lit.columns


def test_11_stage13_14_17_hydration_manifest_reports_and_hashes() -> None:
    v31 = _v31_root()
    root = _repo_root()

    g_man = pd.read_csv(v31 / "authoritative_data" / "genome_checksums_manifest_v3_1.tsv", sep="\t")
    assert len(g_man) == 19
    assert "GCF_000338995.1_Cocsa1_genomic.fna" in set(g_man["filename"])
    assert (v31 / "src" / "hydrate_inputs_v3_1.py").exists()

    for req_doc in (
        "V3_1_Final_Validation_Report.md",
        "V3_1_Change_Log.md",
        "V3_1_Environment_Report.md",
        "README.md",
    ):
        p = v31 / req_doc
        assert p.exists() and p.stat().st_size > 200, f"Missing/empty doc: {p}"

    hashes_df = pd.read_csv(v31 / "hashes_v3_1.tsv", sep="\t")
    assert len(hashes_df) >= 40
    assert not any(str(p).endswith("hashes_v3_1.tsv") for p in hashes_df["relative_path"])
    for _, hr in hashes_df.iterrows():
        fp = root / str(hr["relative_path"])
        assert fp.exists(), f"Manifest file missing on disk: {fp}"
        assert hashlib.sha256(fp.read_bytes()).hexdigest() == str(hr["sha256"]), f"SHA256 mismatch: {fp}"
