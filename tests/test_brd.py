"""Automated test suite for Bipolaris-RNAi Designer (BRD)."""

from __future__ import annotations

from pathlib import Path
import pandas as pd
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from brd.analysis.bootstrap import bootstrap_conservation_and_ranks
from brd.analysis.pareto import compute_pareto_ranks
from brd.analysis.sensitivity import run_sensitivity_scenarios
from brd.analysis.target_scoring import compute_weighted_target_scores
from brd.sequence.cds_validation import validate_cds_translation
from brd.sequence.conservation import analyze_msa_conservation
from brd.sequence.fasta_parser import compute_assembly_stats
from brd.sequence.gff_parser import validate_gff3
from brd.utils.file_validation import validate_fasta, validate_nucleotide_seq


def test_fasta_validation_and_stats(tmp_path: Path) -> None:
    fa = tmp_path / "test.fna"
    fa.write_text(">scf1\nATGCATGCATGC\n>scf2\nGGGGCCCC\n", encoding="utf-8")
    res = validate_fasta(fa, "nucl")
    assert res["valid"] is True
    assert res["records"] == 2
    stats = compute_assembly_stats(fa)
    assert stats["total_length_bp"] == 20
    assert stats["scaffold_count"] == 2
    assert stats["n50_bp"] == 12


def test_fasta_duplicate_header_fails(tmp_path: Path) -> None:
    fa = tmp_path / "dup.fna"
    fa.write_text(">scf1\nATGC\n>scf1\nGGGG\n", encoding="utf-8")
    res = validate_fasta(fa, "nucl")
    assert res["valid"] is False
    assert res["duplicate_ids"] == 1


def test_gff3_invalid_coords(tmp_path: Path) -> None:
    gff = tmp_path / "bad.gff3"
    gff.write_text("scf1\tJGI\tgene\t50\t20\t.\t+\t.\tID=g1\n", encoding="utf-8")
    info = validate_gff3(gff, {"scf1": 100})
    assert info["error_count"] >= 1


def test_cds_translation_validation() -> None:
    # ATG GCT TAA -> MA*
    ok = validate_cds_translation("ATGGCTTAA", expected_protein="MA")
    assert ok["is_complete_cds"] is True
    assert ok["matches_reference_protein"] is True
    bad = validate_cds_translation("ATGTAATAA")  # Internal stop
    assert bad["internal_stop_count"] == 1
    assert bad["is_complete_cds"] is False


def test_nucleotide_conservation_snps_and_kmers() -> None:
    s1 = "ATGGCTAGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
    s2 = "ATGGCTAGTCGTAGTTAGCTAGCTAGCTAGCTAGCTAGCT"  # 1 SNP at 1-based pos 15 -> 25 nt conserved right block
    msa = MultipleSeqAlignment([SeqRecord(Seq(s1), id="iso1"), SeqRecord(Seq(s2), id="iso2")])
    res = analyze_msa_conservation(msa, ref_id="iso1", window_size=21)
    assert res["snp_sites"] == 1
    assert res["n_isolates"] == 2
    assert res["mean_pairwise_identity_pct"] < 100.0
    assert max(b[2] for b in res["conserved_blocks_ge21nt"]) >= 21


def test_pareto_ranking_and_missing_data() -> None:
    df = pd.DataFrame(
        [
            {"candidate_id": "C1", "target": "T1", "conservation_score": 1.0, "off_target_safety_score": 0.9, "target_evidence_score": 0.8, "guide_quality_score": 0.7, "tool_agreement_score": 1.0},
            {"candidate_id": "C2", "target": "T1", "conservation_score": 0.8, "off_target_safety_score": 0.7, "target_evidence_score": 0.6, "guide_quality_score": 0.5, "tool_agreement_score": 0.7},
            {"candidate_id": "C3", "target": "T2", "conservation_score": 0.9, "off_target_safety_score": 1.0, "target_evidence_score": 0.8, "guide_quality_score": 0.7, "tool_agreement_score": 1.0},
        ]
    )
    p = compute_pareto_ranks(df, ["conservation_score", "off_target_safety_score", "target_evidence_score", "guide_quality_score"])
    assert set(p[p["is_pareto_optimal"]]["candidate_id"]) == {"C1", "C3"}
    assert int(p[p["candidate_id"] == "C2"]["pareto_rank"].iloc[0]) == 2

    # Missing evidence conservative scoring check
    tdf = pd.DataFrame(
        [{"target": "T1", "score_conservation": 1.0, "score_off_target": float("nan"), "score_functional": 0.9, "score_expression": 0.8, "score_designability": 1.0}]
    )
    scored = compute_weighted_target_scores(
        tdf,
        {
            "multi_isolate_nucleotide_conservation": 0.25,
            "off_target_safety": 0.25,
            "direct_functional_pathogenicity_evidence": 0.20,
            "infection_stage_expression": 0.15,
            "dsrna_designability": 0.15,
        },
    )
    assert float(scored["known_weight_fraction"].iloc[0]) == 0.75
    assert float(scored["weighted_score_conservative"].iloc[0]) < float(scored["normalized_known_score"].iloc[0])


def test_bootstrap_and_sensitivity() -> None:
    cdf = pd.DataFrame(
        [
            {"candidate_id": "C1", "target": "T1", "off_target_safety_score": 1.0, "target_evidence_score": 0.9, "guide_quality_score": 0.8, "tool_agreement_score": 1.0},
            {"candidate_id": "C2", "target": "T2", "off_target_safety_score": 0.5, "target_evidence_score": 0.5, "guide_quality_score": 0.5, "tool_agreement_score": 0.5},
        ]
    )
    mat = pd.DataFrame(
        [
            {"candidate_id": "C1", "iso1": 1.0, "iso2": 0.98, "iso3": 1.0},
            {"candidate_id": "C2", "iso1": 0.6, "iso2": 0.55, "iso3": 0.65},
        ]
    )
    w = {"nucleotide_conservation": 0.25, "off_target_safety": 0.25, "target_evidence_and_expression": 0.20, "guide_quality_and_accessibility": 0.15, "tool_agreement_and_seq_qc": 0.15}
    boot = bootstrap_conservation_and_ranks(cdf, mat, w, n_iterations=50, seed=42)
    assert len(boot) == 2
    assert boot.iloc[0]["candidate_id"] == "C1"

    cdf["conservation_score"] = [0.99, 0.60]
    sens = run_sensitivity_scenarios(
        cdf,
        "candidate_id",
        {"c": "conservation_score", "o": "off_target_safety_score"},
        {"s1": {"c": 0.8, "o": 0.2}, "s2": {"c": 0.2, "o": 0.8}},
        baseline_name="s1",
    )
    assert int(sens.iloc[0]["rank_s1"]) == 1


def test_nucleotide_seq_validator() -> None:
    ok, st = validate_nucleotide_seq("ATGCATGC")
    assert ok and st["gc_percent"] == 50.0
    bad, _ = validate_nucleotide_seq("ATGCX123")
    assert not bad
