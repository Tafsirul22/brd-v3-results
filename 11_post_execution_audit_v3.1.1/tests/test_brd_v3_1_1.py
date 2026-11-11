"""Comprehensive Verification Suite for Bipolaris-RNAi Designer v3.1.1 (clean-extraction portable)."""

from __future__ import annotations

import hashlib
import importlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO


def _v311_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _repo_root() -> Path:
    return _v311_root().parent


def _v31_root() -> Path:
    return _repo_root() / "11_post_execution_audit_v3.1"


sys.path.insert(0, str(_repo_root()))
sys.path.insert(0, str(_v311_root() / "src"))
sys.path.insert(0, str(_v31_root() / "src"))

from assembly_qc_biopython import analyze_fasta_assembly, compute_n50_l50  # noqa: E402
from hydrate_broccoli_proteomes_v3_1_1 import validate_protein_fasta  # noqa: E402
from hydrate_offtarget_references_v3_1_1 import (  # noqa: E402
    _scan_exact_21mers_fast,
    build_ncbi_download_url,
    validate_nucleotide_fasta,
)
from stage5_broccoli_orthology import parse_broccoli_ortholog_pairs, parse_broccoli_orthologous_groups  # noqa: E402
from verify_manifest_v3_1_1 import verify_manifest  # noqa: E402


def test_01_portability_no_hardcoded_machine_paths() -> None:
    """Fix 4 Action 10: Active code, brd/, and all test suites use relative paths only."""
    root = _repo_root()
    v311 = _v311_root()
    forbidden = ("ASUS " + "BANGLADESH", "New " + "folder (2)", r"C:\Us" + "ers")
    files = (
        list((root / "brd").rglob("*.py"))
        + list((root / "tests").glob("*.py"))
        + list((root / "11_post_execution_audit" / "tests").glob("*.py"))
        + list((v311 / "src").glob("*.py"))
        + list((v311 / "tests").glob("*.py"))
    )
    assert len(files) >= 15
    for fp in files:
        txt = fp.read_text(encoding="utf-8", errors="replace")
        for pat in forbidden:
            assert pat not in txt, f"Hardcoded machine path {pat!r} found in {fp}"


def test_02a_primer3_absent_optional_graceful_degradation() -> None:
    """Fix 1 (absent state): Clean import and unrelated execution succeed without `primer3`; Primer3 calls fail explicitly."""
    root = _repo_root()
    code_absent = r"""
import sys
class _BlockPrimer3:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "primer3" or fullname.startswith("primer3."):
            raise ModuleNotFoundError("No module named 'primer3'", name="primer3")
sys.meta_path.insert(0, _BlockPrimer3())

import brd
import brd.integrations.lucid as lucid
import brd.modules.m07_dsRNA_design as m07
import brd.modules.m10_reporting as m10
from brd.analysis.pareto import compute_pareto_ranks
from brd.utils.file_validation import validate_nucleotide_seq
import pandas as pd

assert lucid.PRIMER3_AVAILABLE is False
assert lucid.is_primer3_available() is False
ok, _ = validate_nucleotide_seq("ATGCATGCATGC")
assert ok
df = pd.DataFrame([{"id": "A", "a": 1.0, "b": 0.5}, {"id": "B", "a": 0.5, "b": 0.4}])
p = compute_pareto_ranks(df, ["a", "b"])
assert list(p["pareto_rank"]) == [1, 2]

try:
    lucid.design_amplification_primers("T", "ATGC" * 60, 0, 200, 200)
    raise AssertionError("Expected RuntimeError when invoking Primer3 while absent")
except RuntimeError as exc:
    assert "primer3-py" in str(exc)

if lucid.VIENNARNA_AVAILABLE and lucid.SIFI21_AVAILABLE:
    wdf, _ = lucid.design_lucid_windows("T", "ATGCATGC" * 40, min_len=200, max_len=250, design_primers=False)
    assert len(wdf) > 0 and (wdf["primer_left"] == "NA").all()
print("PRIMER3_ABSENT_ISOLATED_OK")
"""
    p_abs = subprocess.run([sys.executable, "-c", code_absent], cwd=str(root), capture_output=True, text=True, check=True)
    assert "PRIMER3_ABSENT_ISOLATED_OK" in p_abs.stdout


def test_02b_primer3_installed_live_design() -> None:
    """Fix 1 (installed state): Verify Primer3 PCR primer design executes cleanly when primer3-py is available."""
    import brd.integrations.lucid as lucid  # noqa: PLC0415

    if not lucid.PRIMER3_AVAILABLE:
        pytest.skip("Optional primer3-py not installed in current minimal environment")
    assert lucid.is_primer3_available() is True
    seq = (
        "ATGCCACCAGCAGGAAGCGGCAGCTCCAGGAAAATCTCCTTCAACGTCTCGGAGCAGTACGACATTCAGGATGTTGTAGGAGAGGGCGCATACGGTGTCG"
        "TATGCTCTGCCTTGCACAAGCCATCGGGCCAGAAGGTTGCCATTAAGAAGATCACGCCCTTTGATCACTCCATGTTCTGCTTGAGAACACTGCGAGAGAT"
        "GAAGCTGCTGCGCTACTTCAACCACGAAAACATCATCTCCATCCTGGACATCCAAAAGCCCCGCAACTACGAGACCTTCACCGAGGTCTACCTAATTCAG"
    )
    prim = lucid.design_amplification_primers("CsTest", seq, 0, 200, 200)
    assert len(prim["primer_left"]) >= 18 and prim["primer_left"] != "NA"


def test_03_optional_viennarna_absent_import_and_runtime_guard() -> None:
    """Verify core BRD imports and unrelated modules run cleanly even if RNA, sifi21_plus, and primer3 are all absent."""
    root = _repo_root()
    code_no_opt = r"""
import sys
class _BlockOptional:
    BLOCK = {"primer3", "RNA", "sifi21_plus", "snakemake"}
    def find_spec(self, fullname, path=None, target=None):
        top = fullname.split(".")[0]
        if top in self.BLOCK:
            raise ModuleNotFoundError(f"No module named '{top}'", name=top)
sys.meta_path.insert(0, _BlockOptional())

import brd.integrations.lucid as lucid
import brd.integrations.rnahybrid as rnahybrid
import brd.integrations.siren as siren
import brd.modules.m07_dsRNA_design as m07
import brd.modules.m10_reporting as m10
from brd.sequence.cds_validation import validate_cds_translation

assert lucid.PRIMER3_AVAILABLE is False
assert lucid.VIENNARNA_AVAILABLE is False
assert rnahybrid.is_viennarna_available() is False
assert validate_cds_translation("ATGGCTTAA", expected_protein="MA")["is_complete_cds"] is True

try:
    rnahybrid.compute_duplex_mfe("AUGCAUGCAUGCAUGCAUGCA", "AUGCAUGCAUGCAUGCAUGCA")
    raise AssertionError("Expected RuntimeError when ViennaRNA absent")
except RuntimeError as exc:
    assert "ViennaRNA" in str(exc)

try:
    lucid.compute_unpaired_probs("ATGC" * 60)
    raise AssertionError("Expected RuntimeError when ViennaRNA absent")
except RuntimeError as exc:
    assert "ViennaRNA" in str(exc)
print("ALL_OPTIONAL_ABSENT_ISOLATED_OK")
"""
    p_res = subprocess.run([sys.executable, "-c", code_no_opt], cwd=str(root), capture_output=True, text=True, check=True)
    assert "ALL_OPTIONAL_ABSENT_ISOLATED_OK" in p_res.stdout


def test_04_stage0_v31_baseline_preserved_and_snapshot_valid() -> None:
    """Confirm Stage 0 inspection snapshot, unchanged v3.1 scientific tables, and preserved rankings."""
    v311 = _v311_root()
    v31 = _v31_root()
    snap = json.loads((v311 / "stage0_audit" / "v3_1_baseline_snapshot.json").read_text(encoding="utf-8"))
    assert snap["v31_archive_crc_test"] == "PASS (None)"
    assert snap["v31_archive_file_count"] == 340

    reg_now = pd.read_csv(v31 / "stage10_ranking" / "independent_region_ranking_v3_1.tsv", sep="\t")
    clean_now = reg_now[reg_now["passes_zero_exact_21mer_filter"]]["candidate_id"].tolist()
    clean_snap = [r["candidate_id"] for r in snap["clean_shortlist_8"]]
    assert clean_now == clean_snap
    assert clean_now[0] == "BRD_CsSTE11_CAND013"


def test_05_stage2_proteome_recovery_validation_and_broccoli_rbh(tmp_path: Path) -> None:
    """Fix 2: Verify proteome recovery manifests, FASTA validator edge cases, and Broccoli vs RBH separation."""
    v311 = _v311_root()
    v31 = _v31_root()
    s2 = v311 / "stage2_proteome_hydration"

    pman = pd.read_csv(s2 / "broccoli_proteome_manifest_v3_1_1.tsv", sep="\t")
    pval = pd.read_csv(s2 / "broccoli_proteome_validation.tsv", sep="\t")
    assert len(pman) == 3
    assert set(pman["species_id"]) == {"Bcin_ASM83294v1", "Bsor_ND90Pr", "Bsor_Sacc._Shoemaker"}
    assert (pman["expected_sha256"] == pman["actual_sha256"]).all()
    assert (pval["matches_historical_sha256"].astype(bool)).all()
    comb_row = pval[pval["file"] == "all_proteomes_combined_v3.faa"].iloc[0]
    assert int(comb_row["protein_count"]) == 36714

    # Positive + negative FASTA validator unit tests
    good_faa = tmp_path / "good.faa"
    good_faa.write_text(">SpA||P1\nMAKT\n>SpA||P2\nMGSS\n", encoding="utf-8")
    v_ok = validate_protein_fasta(good_faa, expected_count=2, require_species_prefix="SpA")
    assert v_ok["protein_count"] == 2 and v_ok["duplicate_id_count"] == 0

    dup_faa = tmp_path / "dup.faa"
    dup_faa.write_text(">P1\nMAKT\n>P1\nMGSS\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="Duplicate"):
        validate_protein_fasta(dup_faa)

    # Broccoli v2 & DIAMOND RBH preserved separation
    targets = pd.read_csv(v31 / "stage5_broccoli_orthology" / "broccoli_vs_diamond_target_genes_v3_1.tsv", sep="\t").set_index("target")
    assert targets.loc["CsFUS3", "broccoli_og"] == "OG_1315"
    assert targets.loc["CsPKS1", "broccoli_og"] == "OG_4251"
    assert targets.loc["CsSTE11", "broccoli_og"] == "OG_1786"
    assert targets.loc["CsSp1", "broccoli_og"] == "OG_9092"


def test_06_stage3_offtarget_hydration_validation_and_guards(tmp_path: Path) -> None:
    """Fix 3: Verify off-target hydration manifest, SHA256 validator fail-fast behavior, and k-mer smoke scan."""
    v311 = _v311_root()
    s3 = v311 / "stage3_offtarget_hydration"
    oman = pd.read_csv(s3 / "offtarget_reference_manifest_v3_1_1.tsv", sep="\t")
    oval = pd.read_csv(s3 / "offtarget_reference_validation.tsv", sep="\t")
    assert len(oman) == 18
    assert (oval["validation_status"] == "PASS").all()
    assert len(oval) == 33  # 18 transcriptomes + 15 genomic assemblies (wheat/barley/maize multi-GB nuclear genomes excluded)
    log_txt = (s3 / "offtarget_reference_retrieval_log.txt").read_text(encoding="utf-8")
    assert "invariant_ok=True" in log_txt and "BUILT_OK" in log_txt

    assert "include_annotation_type=RNA_FASTA,GENOME_FASTA" in build_ncbi_download_url("GCF_034140825.1", "rna.fna", True)
    assert "include_annotation_type=CDS_FASTA" in build_ncbi_download_url("GCA_002897155.2", "cds_from_genomic.fna", False)

    # Nucleotide validator mismatch failure guard
    fa = tmp_path / "nt.fna"
    fa.write_text(">seq1\nATGCATGCATGCATGCATGCATGC\n", encoding="utf-8")
    good_sha = hashlib.sha256(fa.read_bytes()).hexdigest()
    res_ok = validate_nucleotide_fasta(fa, good_sha, fa.stat().st_size, 1)
    assert res_ok["sequence_count"] == 1
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        validate_nucleotide_fasta(fa, "0" * 64, fa.stat().st_size, 1)

    # Exact dual-strand 21-mer scan smoke check
    guide = "ATGGCCATGCTCGCCTCAAAG"
    fa2 = tmp_path / "scan.fna"
    fa2.write_text(f">rec1\nAAAA{guide}TTTT\n", encoding="utf-8")
    assert _scan_exact_21mers_fast(fa2, {guide}) == 1


def test_07_assembly_qc_n50_l50_and_conservation_preserved() -> None:
    """Preserve Biopython assembly QC & conservation results with deterministic unit tests."""
    v31 = _v31_root()
    assert compute_n50_l50([20, 50, 30]) == (50, 1)
    assert compute_n50_l50([10, 40, 20, 30]) == (30, 2)
    asm = pd.read_csv(v31 / "stage3_assembly_qc" / "assembly_qc_summary_v3_1.tsv", sep="\t")
    assert len(asm) == 19 and (asm["status"] == "OK").all()
    gc = pd.read_csv(v31 / "stage6_conservation" / "gene_conservation_v3_1.tsv", sep="\t")
    assert len(gc) == 4


def test_08_candidate_coordinate_audit_consistency() -> None:
    """Fix 6: Verify 50 recorded = 49 contiguous CDS + 1 synthetic construct (CAND009) + 0 unresolved consistently everywhere."""
    v311 = _v311_root()
    s6 = v311 / "stage6_coordinates"
    audit_df = pd.read_csv(s6 / "candidate_coordinate_audit_v3_1_1.tsv", sep="\t")
    summ_df = pd.read_csv(s6 / "candidate_coordinate_summary_v3_1_1.tsv", sep="\t").set_index("metric")

    assert len(audit_df) == 50
    assert int(summ_df.loc["candidates_with_recorded_coordinate_assignments", "count"]) == 50
    assert int(summ_df.loc["candidates_matching_contiguous_authoritative_cds_substrings", "count"]) == 49
    assert int(summ_df.loc["candidates_stitched_or_synthetic_constructs", "count"]) == 1
    assert int(summ_df.loc["candidates_with_unresolved_coordinate_mapping", "count"]) == 0
    assert (audit_df["coordinate_convention"] == "1-based_inclusive").all()
    assert list(audit_df.loc[audit_df["is_stitched_or_synthetic"], "candidate_id"]) == ["BRD_CsSTE11_CAND009"]

    for doc_name in ("V3_1_1_Change_Log.md", "V3_1_1_Final_Validation_Report.md"):
        doc_p = v311 / doc_name
        if doc_p.exists():
            txt = doc_p.read_text(encoding="utf-8")
            assert "50/50 CDS" not in txt and "50/50 contiguous" not in txt
            assert "49/50" in txt and "BRD_CsSTE11_CAND009" in txt


def test_09_ranking_pareto_and_bootstrap_preservation() -> None:
    """Preserve exact v3.1 candidate rankings, Pareto tier separation, sensitivity, and bootstrap invariants."""
    v31 = _v31_root()
    rk = pd.read_csv(v31 / "stage10_ranking" / "candidate_ranking_v3_1.tsv", sep="\t")
    reg_rk = pd.read_csv(v31 / "stage10_ranking" / "independent_region_ranking_v3_1.tsv", sep="\t")
    assert len(rk) == 50 and len(reg_rk) == 12
    for _, r in rk.iterrows():
        recalc = round(
            0.30 * float(r["biological_target_relevance_score"])
            + 0.35 * float(r["dsrna_targetability_score"])
            + 0.35 * float(r["sequence_selectivity_score"]),
            4,
        )
        assert abs(recalc - float(r["heuristic_composite_prioritization_index"])) <= 1e-4
    clean = reg_rk[reg_rk["passes_zero_exact_21mer_filter"]]
    expected = [
        ("BRD_CsSTE11_CAND013", 0.9418),
        ("BRD_CsSTE11_CAND014", 0.9367),
        ("BRD_CsSTE11_CAND010", 0.8861),
        ("BRD_CsPKS1_CAND012", 0.7979),
        ("BRD_CsPKS1_CAND010", 0.7762),
        ("BRD_CsPKS1_CAND013", 0.7841),
        ("BRD_CsPKS1_CAND006", 0.7836),
        ("BRD_CsSp1_CAND008", 0.7363),
    ]
    actual = list(zip(clean["candidate_id"], clean["heuristic_composite_prioritization_index"].round(4)))
    assert actual == expected


def test_10_manifest_completeness_and_verifier_script(tmp_path: Path) -> None:
    """Fix 5: Manifest verifies 100% of release files (including logs & test reports) and detects tampering/unlisted files."""
    root = _repo_root()
    v311 = _v311_root()
    mpath = v311 / "hashes_v3_1_1.tsv"
    assert mpath.exists()

    res = verify_manifest(repo_root=root, manifest_path=mpath, check_unlisted_on_disk=True)
    assert res["ok"], f"Manifest verification failed: {res}"
    assert res["self_included"] == []
    assert res["unlisted_on_disk"] == []

    listed = set(pd.read_csv(mpath, sep="\t")["relative_path"])
    for req_log in (
        "11_post_execution_audit_v3.1.1/V3_1_1_Execution_Log.txt",
        "11_post_execution_audit_v3.1.1/V3_1_1_Test_Results.txt",
        "11_post_execution_audit_v3.1.1/V3_1_1_Test_Summary.md",
        "11_post_execution_audit_v3.1.1/V3_1_1_Final_Validation_Report.md",
    ):
        assert req_log in listed, f"Expected release artifact omitted from hashes_v3_1_1.tsv: {req_log}"

    # Synthetic negative verification check: unlisted file must cause verify_manifest -> False and non-zero CLI exit
    fake_root = tmp_path / "mini_repo"
    (fake_root / "brd").mkdir(parents=True)
    (fake_root / "config").mkdir(parents=True)
    (fake_root / "config" / "project_configuration.yaml").write_text("project: {}\n", encoding="utf-8")
    (fake_root / "11_post_execution_audit_v3.1.1" / "src").mkdir(parents=True)
    (fake_root / "README.md").write_text("hi\n", encoding="utf-8")
    from verify_manifest_v3_1_1 import generate_hashes_v3_1_1  # noqa: PLC0415

    hm = generate_hashes_v3_1_1(fake_root)
    assert verify_manifest(fake_root, hm, check_unlisted_on_disk=True)["ok"] is True
    (fake_root / "tests").mkdir(exist_ok=True)
    (fake_root / "tests" / "unlisted_extra.py").write_text("x=1\n", encoding="utf-8")
    bad_res = verify_manifest(fake_root, hm, check_unlisted_on_disk=True)
    assert bad_res["ok"] is False and "tests/unlisted_extra.py" in bad_res["unlisted_on_disk"]


def test_11_historical_tests_and_clean_env_inventory() -> None:
    """Fix 4: Confirm historical v2, v3, core, v3.1, and v3.1.1 test files & dependency manifests are packaged."""
    root = _repo_root()
    v311 = _v311_root()
    for tf in (
        "tests/test_brd.py",
        "11_post_execution_audit/tests/test_audit_v2.py",
        "11_post_execution_audit_v3/tests/test_brd_v3.py",
        "11_post_execution_audit_v3.1/tests/test_brd_v3_1.py",
        "11_post_execution_audit_v3.1.1/tests/test_brd_v3_1_1.py",
        "11_post_execution_audit_v3.1.1/V3_1_1_Test_Inventory.tsv",
        "11_post_execution_audit_v3.1.1/requirements.txt",
        "11_post_execution_audit_v3.1.1/requirements-optional.txt",
    ):
        p = root / tf
        assert p.exists() and p.stat().st_size > 100, f"Missing required test or env file: {tf}"
    req_txt = (v311 / "requirements.txt").read_text(encoding="utf-8")
    assert "primer3-py" in (v311 / "requirements-optional.txt").read_text(encoding="utf-8")
    assert "\nprimer3-py==" not in req_txt
