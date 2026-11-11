"""Automated cross-file validation tests (Section 10) for BRD v2 post-execution audit."""

import json
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

OUT_DIR = Path(r"C:\Users\ASUS BANGLADESH\Desktop\New folder (2)\Output")
AUDIT_DIR = OUT_DIR / "11_post_execution_audit"
FIN_V1 = OUT_DIR / "10_final_reports"


def test_01_fasta_sequences_match_tsv_v2():
    cands = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidates_v2.tsv", sep="\t")
    fa = {r.id: str(r.seq).upper() for r in SeqIO.parse(str(AUDIT_DIR / "BRD_dsRNA_Candidates_v2.fasta"), "fasta")}
    assert len(cands) == 50
    assert len(fa) == 50
    assert set(cands["candidate_id"]) == set(fa.keys())
    for _, r in cands.iterrows():
        assert fa[r["candidate_id"]] == str(r["sequence"]).upper()


def test_02_candidate_lengths_and_gc_correct():
    cands = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidates_v2.tsv", sep="\t")
    for _, r in cands.iterrows():
        seq = str(r["sequence"]).upper()
        assert int(r["candidate_length_nt"]) == len(seq)
        expected_gc = round((seq.count("G") + seq.count("C")) / len(seq) * 100.0, 2)
        assert abs(float(r["gc_percent"]) - expected_gc) <= 0.01
        assert int(r["n_bases"]) == seq.count("N")


def test_03_target_coordinates_valid_and_substrings():
    cands = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidates_v2.tsv", sep="\t")
    tgs = pd.read_csv(FIN_V1 / "BRD_Target_Gene_Summary.tsv", sep="\t")
    for _, r in tgs.iterrows():
        assert int(r["start"]) < int(r["end"])
    refs = {r.id.split("|")[0]: str(r.seq).upper() for r in SeqIO.parse(str(FIN_V1 / "BRD_Target_Gene_Sequences.fasta"), "fasta")}
    for _, r in cands.iterrows():
        if bool(r["exact_cds_substring"]):
            s, e = int(r["candidate_start_nt"]), int(r["candidate_end_nt"])
            assert 1 <= s < e <= len(refs[r["target"]])
            assert refs[r["target"]][s - 1 : e] == str(r["sequence"]).upper()


def test_04_accessions_traceable_v2():
    ident = pd.read_csv(AUDIT_DIR / "BRD_Target_Gene_Identity_Evidence_v2.tsv", sep="\t")
    pks1 = pd.read_csv(AUDIT_DIR / "BRD_CsPKS1_Record_Audit_v2.tsv", sep="\t")
    prim = pd.read_csv(AUDIT_DIR / "BRD_Target_Primer_Alignments_v2.tsv", sep="\t")
    assert set(ident["target"]) == {"CsPKS1", "CsFUS3", "CsSTE11", "CsSp1"}
    assert len(pks1) == 9
    assert "HQ830033.1" in set(pks1["accession"])
    for k in range(134, 142):
        assert f"KF769{k}.1" in set(pks1["accession"])
    assert pks1["maps_to_COCSADRAFT_157080"].all()
    assert len(prim) == 12
    cssp1_prim = prim[prim["target"] == "CsSp1"]
    assert len(cssp1_prim) == 5
    assert (cssp1_prim["identity_pct"] == 100.0).all()


def test_05_candidate_ids_match_across_m07_m08_m09_v2():
    c_ids = set(pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidates_v2.tsv", sep="\t")["candidate_id"])
    o_ids = set(pd.read_csv(AUDIT_DIR / "BRD_Off_Target_Analysis_v2.tsv", sep="\t")["candidate_id"])
    r1_ids = set(pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidate_Ranking_v2.tsv", sep="\t")["candidate_id"])
    r2_ids = set(pd.read_csv(AUDIT_DIR / "final_candidate_ranking_v2.tsv", sep="\t")["candidate_id"])
    b_ids = set(pd.read_csv(AUDIT_DIR / "bootstrap_results_v2.tsv", sep="\t")["candidate_id"])
    s_ids = set(pd.read_csv(AUDIT_DIR / "sensitivity_analysis_v2.tsv", sep="\t")["candidate_id"])
    cf_ids = set(pd.read_csv(AUDIT_DIR / "candidate_confidence_v2.tsv", sep="\t")["candidate_id"])
    assert c_ids == o_ids == r1_ids == r2_ids == b_ids == s_ids == cf_ids


def test_06_no_stale_pending_m08():
    for fn in [
        "BRD_dsRNA_Candidates_v2.tsv",
        "BRD_Off_Target_Analysis_v2.tsv",
        "BRD_dsRNA_Candidate_Ranking_v2.tsv",
        "final_candidate_ranking_v2.tsv",
    ]:
        txt = (AUDIT_DIR / fn).read_text(encoding="utf-8")
        assert "PENDING_M08" not in txt, f"Found stale PENDING_M08 in {fn}"


def test_07_ranking_values_exact_reproducible():
    rk = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidate_Ranking_v2.tsv", sep="\t")
    fk = pd.read_csv(AUDIT_DIR / "final_candidate_ranking_v2.tsv", sep="\t")
    pd.testing.assert_frame_equal(rk, fk)
    recomp = (
        rk["conservation_score"] * 0.25
        + rk["off_target_safety_score_v2"] * 0.25
        + rk["target_evidence_score"] * 0.20
        + rk["guide_quality_score"] * 0.15
        + rk["tool_agreement_score"] * 0.15
    ).round(4)
    assert (rk["composite_score_v2"].round(4) == recomp).all()
    assert list(rk["final_rank_v2"]) == list(range(1, len(rk) + 1))


def test_08_exact_off_target_matches_reported_and_excluded():
    ex = pd.read_csv(AUDIT_DIR / "off_target_exact_21mer_matches_v2.tsv", sep="\t")
    ot = pd.read_csv(AUDIT_DIR / "BRD_Off_Target_Analysis_v2.tsv", sep="\t").set_index("candidate_id")
    rk = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidate_Ranking_v2.tsv", sep="\t").set_index("candidate_id")
    assert not ex.empty
    for cid in ex["candidate_id"].unique():
        total_ex = (
            int(ot.loc[cid, "wheat_exact_21mer_transcripts"])
            + int(ot.loc[cid, "barley_exact_21mer_transcripts"])
            + int(ot.loc[cid, "beneficial_fungi_exact_21mer_transcripts"])
            + int(ot.loc[cid, "fungal_paralog_exact_21mer_transcripts"])
        )
        assert total_ex > 0
        assert not bool(rk.loc[cid, "passes_computational_filters_v2"])
        assert "Excluded" in str(rk.loc[cid, "computational_priority_category_v2"])


def test_09_missing_linux_tools_clearly_labelled():
    sw = pd.read_csv(AUDIT_DIR / "BRD_Software_Integration_Audit_v2.tsv", sep="\t").set_index("software_or_component")
    for missing_tool in ["QUAST CLI", "BUSCO", "OrthoFinder (MCL/FastME)", "Funannotate / BRAKER3", "RNAhybrid"]:
        assert sw.loc[missing_tool, "audit_status_v2"] == "NOT_RUN"
    qc = pd.read_csv(OUT_DIR / "02_genome_qc" / "assembly_quality_summary.tsv", sep="\t")
    assert (qc["busco_status"] == "NOT_RUN_LINUX_ONLY_NO_WSL_DISTRO").all()


def test_10_final_report_agrees_with_v2_tsv_and_v1_preserved():
    # Confirm original v1 BRD_dsRNA_Candidates.tsv untouched (retains original v1 PENDING_M08 marker)
    v1_c = pd.read_csv(FIN_V1 / "BRD_dsRNA_Candidates.tsv", sep="\t")
    assert (v1_c["off_target_status"] == "PENDING_M08").all()
    for rfn in [
        "BRD_Final_Research_Report_v2.md",
        "BRD_Final_Research_Report_v2.html",
        "BRD_Final_Research_Report_v2.pdf",
        "BRD_Execution_Summary_v2.txt",
    ]:
        p = AUDIT_DIR / rfn
        assert p.exists() and p.stat().st_size > 1000
    rk = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidate_Ranking_v2.tsv", sep="\t")
    dual = pd.read_csv(AUDIT_DIR / "target_biological_vs_dsrna_prioritisation_v2.tsv", sep="\t")
    top1 = str(rk[rk["passes_computational_filters_v2"]].iloc[0]["candidate_id"])
    n_pass = int(rk["passes_computational_filters_v2"].sum())
    bio_expected = ", ".join(
        f"#{r['biological_target_rank']} {r['target']} ({r['biological_target_score']})"
        for _, r in dual.sort_values("biological_target_rank").iterrows()
    )
    dsrna_expected = ", ".join(
        f"#{r['dsrna_targetability_rank']} {r['target']} ({r['best_passing_candidate_score_v2']})"
        for _, r in dual.sort_values("dsrna_targetability_rank").iterrows()
    )
    md_txt = (AUDIT_DIR / "BRD_Final_Research_Report_v2.md").read_text(encoding="utf-8")
    ex_txt = (AUDIT_DIR / "BRD_Execution_Summary_v2.txt").read_text(encoding="utf-8")
    assert top1 in md_txt and top1 in ex_txt
    assert f"{n_pass}/50" in md_txt and f"{n_pass}/50" in ex_txt
    assert bio_expected in md_txt and bio_expected in ex_txt
    assert dsrna_expected in md_txt and dsrna_expected in ex_txt


def test_11_both_orientations_screened_and_all_exact_21mers_reported():
    cands = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidates_v2.tsv", sep="\t")
    ex = pd.read_csv(AUDIT_DIR / "off_target_exact_21mer_matches_v2.tsv", sep="\t")
    # Note: off_target_near_guide_matches_v2.tsv is the exact filename written by run_post_execution_audit.py
    near = pd.read_csv(AUDIT_DIR / "off_target_near_guide_matches_v2.tsv", sep="\t")
    cand_seqs = {r["candidate_id"]: str(r["sequence"]).upper() for _, r in cands.iterrows()}

    # 1. Verify both orientations (sense + antisense 21-mers) map back to each reported exact hit
    for _, r in ex.iterrows():
        cs = cand_seqs[r["candidate_id"]]
        rc_cs = str(Seq(cs).reverse_complement())
        km = str(r["matched_21mer"]).upper()
        assert km in cs or km in rc_cs, f"{km} not in either strand of {r['candidate_id']}"

    # 2. Independently scan all beneficial fungi FASTAs across BOTH candidate orientations
    all_kmers = set()
    cand_to_kmers = {}
    for cid, cs in cand_seqs.items():
        rc_cs = str(Seq(cs).reverse_complement())
        kms = {cs[i : i + 21] for i in range(len(cs) - 20)} | {rc_cs[i : i + 21] for i in range(len(rc_cs) - 20)}
        cand_to_kmers[cid] = kms
        all_kmers.update(kms)

    manifest = pd.read_csv(OUT_DIR / "08_off_target_screening" / "blast_databases" / "database_manifest.tsv", sep="\t")
    indep_hits = set()
    for _, mr in manifest[manifest["group"] == "beneficial_fungi"].iterrows():
        for rec in SeqIO.parse(str(mr["fasta_file"]), "fasta"):
            s = str(rec.seq).upper()
            for i in range(max(0, len(s) - 20)):
                sub = s[i : i + 21]
                if sub in all_kmers:
                    for cid, kms in cand_to_kmers.items():
                        if sub in kms:
                            indep_hits.add((cid, rec.id, sub))

    bf_ex = ex[ex["database_group"] == "beneficial_fungi"]
    reported_hits = set(zip(bf_ex["candidate_id"], bf_ex["subject_accession"], bf_ex["matched_21mer"]))
    assert indep_hits == reported_hits

    # 3. Verify guide BLAST near-matches recorded both orientations (plus & minus)
    assert {"plus", "minus"}.issubset(set(near["sstrand"].dropna().astype(str).unique()))


def test_12_biological_scoring_rubric_and_expression_exclusion():
    rubric = json.loads((AUDIT_DIR / "biological_target_scoring_rubric_v2.json").read_text(encoding="utf-8"))
    weights = {c["criterion_id"]: float(c["weight"]) for c in rubric["criteria"]}
    assert abs(sum(weights.values()) - 1.0) <= 1e-6

    bio_p = pd.read_csv(AUDIT_DIR / "BRD_Target_Biological_Prioritisation_v2.tsv", sep="\t").set_index("target")
    ev = pd.read_csv(AUDIT_DIR / "BRD_Biological_Scoring_Evidence_Matrix_v2.tsv", sep="\t")
    expr = pd.read_csv(AUDIT_DIR / "GSE130517_expression_distribution_audit_v2.tsv", sep="\t").set_index("target")
    rk = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidate_Ranking_v2.tsv", sep="\t")

    assert (bio_p["expression_excluded_from_composite"] == True).all()
    assert (rk["expression_excluded_from_composite"] == True).all()

    for tname, row in bio_p.iterrows():
        sub_ev = ev[ev["target"] == tname]
        assert set(sub_ev["criterion_id"]) == set(weights.keys())
        recomp = round(sum(float(r["weight"]) * float(r["criterion_score_0_1"]) for _, r in sub_ev.iterrows()), 4)
        assert abs(float(row["biological_target_score"]) - recomp) <= 1e-4

        cand_ev = rk[rk["target"] == tname]["target_evidence_score"].unique()
        assert len(cand_ev) == 1 and abs(float(cand_ev[0]) - float(row["biological_target_score"])) <= 1e-4

    # Verify CsFUS3 stays #1 biologically despite having 0 clean dsRNA candidates (#4 dsRNA targetability)
    assert int(bio_p.loc["CsFUS3", "biological_target_rank"]) == 1
    assert int(bio_p.loc["CsFUS3", "dsrna_targetability_rank"]) == 4

    # Verify GSE130517 single-read spike provenance documented explicitly
    assert str(expr.loc["CsSp1", "v1_peak_source_sample"]).startswith("3A_")
    assert str(expr.loc["CsSTE11", "v1_peak_source_sample"]).startswith("2B_")
    assert (expr["valid_3B_4B_5B_raw_reads_sum"] > 0).all()


def test_13_filtering_funnel_and_required_deliverables_present():
    funnel = pd.read_csv(AUDIT_DIR / "off_target_filtering_funnel_v2.tsv", sep="\t")
    rk = pd.read_csv(AUDIT_DIR / "BRD_dsRNA_Candidate_Ranking_v2.tsv", sep="\t")
    cumul = list(funnel["candidates_passing_cumulatively"])
    assert cumul[0] == 50
    assert all(cumul[i] >= cumul[i + 1] for i in range(len(cumul) - 1))
    hard_row = funnel[funnel["filter_step"] == "5_Conservation_ge85pct_Hard_Filter_Pass"].iloc[0]
    assert int(hard_row["candidates_passing_cumulatively"]) == int(rk["passes_computational_filters_v2"].sum())

    md_txt = (AUDIT_DIR / "BRD_Final_Research_Report_v2.md").read_text(encoding="utf-8")
    ex_txt = (AUDIT_DIR / "BRD_Execution_Summary_v2.txt").read_text(encoding="utf-8")
    for _, fr in funnel.iterrows():
        assert str(fr["filter_step"]) in md_txt and str(fr["filter_step"]) in ex_txt

    required = [
        "BRD_Final_Research_Report_v2.md",
        "BRD_Final_Research_Report_v2.html",
        "BRD_Final_Research_Report_v2.pdf",
        "BRD_Execution_Summary_v2.txt",
        "BRD_Target_Biological_Prioritisation_v2.tsv",
        "BRD_dsRNA_Candidate_Ranking_v2.tsv",
        "BRD_Off_Target_Analysis_v2.tsv",
        "BRD_dsRNA_Candidates_v2.fasta",
        "BRD_dsRNA_Candidates_v2.tsv",
        "BRD_Biological_Scoring_Evidence_Matrix_v2.tsv",
        "GSE130517_expression_distribution_audit_v2.tsv",
        "off_target_filtering_funnel_v2.tsv",
        "biological_target_scoring_rubric_v2.json",
        "off_target_classification_criteria.json",
    ]
    for fn in required:
        p = AUDIT_DIR / fn
        assert p.exists() and p.stat().st_size > 100, f"Missing or empty deliverable: {fn}"
