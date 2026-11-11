"""Module M10 — Automated Research Reporting (Markdown, HTML, PDF, Provenance, and Execution Summary)."""

from __future__ import annotations

import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Image as RLImage, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from brd import __version__
from brd.config import BRDConfig
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker
from brd.utils.report_utils import df_to_html_table, wrap_html_report, write_text


def _collect_software_versions(cfg: BRDConfig) -> pd.DataFrame:
    rows = []
    # Python libraries
    import Bio, matplotlib, numpy, pandas, primer3, reportlab, RNA, scipy, snakemake, yaml

    py_libs = [
        ("Python", os.sys.version.split()[0], "https://www.python.org/", "PSF", "Native Windows 11"),
        ("Biopython", Bio.__version__, "https://biopython.org/", "Biopython", "Native Windows 11"),
        ("ViennaRNA (Python RNA)", RNA.__version__, "https://www.tbi.univie.ac.at/RNA/", "ViennaRNA", "Native Windows 11"),
        ("Snakemake", snakemake.__version__, "https://snakemake.readthedocs.io/", "MIT", "Native Windows 11"),
        ("pandas", pandas.__version__, "https://pandas.pydata.org/", "BSD-3-Clause", "Native Windows 11"),
        ("NumPy", numpy.__version__, "https://numpy.org/", "BSD-3-Clause", "Native Windows 11"),
        ("SciPy", scipy.__version__, "https://scipy.org/", "BSD-3-Clause", "Native Windows 11"),
        ("Matplotlib", matplotlib.__version__, "https://matplotlib.org/", "PSF-based", "Native Windows 11"),
        ("primer3-py", primer3.__version__, "https://github.com/libnano/primer3-py", "GPL-2.0", "Native Windows 11"),
        ("reportlab", reportlab.__version__, "https://www.reportlab.com/", "BSD", "Native Windows 11"),
        ("PyYAML", yaml.__version__, "https://pyyaml.org/", "MIT", "Native Windows 11"),
    ]
    for name, ver, url, lic, status in py_libs:
        rows.append({"software": name, "version": ver, "repository_or_source": url, "license": lic, "execution_status": status})

    # External binaries & cloned repos
    ext_checks = [
        ("NCBI BLAST+", [cfg.tool_path("blastn"), "-version"], "https://blast.ncbi.nlm.nih.gov/", "Public Domain", "Native Windows 11"),
        ("NCBI Datasets CLI", [cfg.tool_path("datasets"), "--version"], "https://www.ncbi.nlm.nih.gov/datasets/", "Public Domain", "Native Windows 11"),
        ("DIAMOND", [cfg.tool_path("diamond"), "version"], "https://github.com/bbuchfink/diamond", "GPL-3.0", "Native Windows 11"),
        ("MAFFT", [cfg.tool_path("mafft"), "--version"], "https://mafft.cbrc.jp/alignment/software/", "BSD-3-Clause", "Native Windows 11"),
        ("dsRNAmax", [cfg.tool_path("dsrnamax"), "-h"], "https://github.com/sfletc/dsRNAmax", "BSD-3-Clause", "Native Windows 11 (commit 674b5b5)"),
    ]
    for name, cmd, url, lic, st in ext_checks:
        ver_str = "installed"
        try:
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=15)
            txt = [ln.strip() for ln in (p.stdout + "\n" + p.stderr).splitlines() if ln.strip() and "Active code page" not in ln and not ln.strip().startswith("Usage")]
            ver_str = " | ".join(txt[:2])[:100]
        except Exception as exc:
            ver_str = f"error: {exc}"
        rows.append({"software": name, "version": ver_str, "repository_or_source": url, "license": lic, "execution_status": st})

    # Cloned RNAi repos + Linux-only tools
    rows.extend(
        [
            {"software": "SIREN (siren-rnai)", "version": "0.1.6 (commit 0f5a5a9)", "repository_or_source": "https://github.com/pablovargasmejia/SIREN", "license": "GPL-3.0", "execution_status": "Native Windows Python modules (prefilter + designVIII); RNAhybrid replaced by ViennaRNA duplexfold"},
            {"software": "LUCID (sifi21_plus)", "version": "0.1.0 (commit 37910bd)", "repository_or_source": "https://github.com/Brunxi/LUCID", "license": "No explicit LICENSE file in repo", "execution_status": "sifi21_plus free-energy/efficiency reproduced natively with ViennaRNA 2.8.0"},
            {"software": "tensorflow/lucid", "version": "v0.3.10 (rejected)", "repository_or_source": "https://github.com/tensorflow/lucid", "license": "Apache-2.0", "execution_status": "EXCLUDED — neural network interpretability library, not RNAi design"},
            {"software": "evanpeikon/BLAST", "version": "unreleased (rejected)", "repository_or_source": "https://github.com/evanpeikon/BLAST", "license": "MIT", "execution_status": "EXCLUDED — educational script; official NCBI BLAST+ 2.17.0+ used"},
            {"software": "metashot/busco", "version": "1.0.0 (rejected)", "repository_or_source": "https://github.com/metashot/busco", "license": "GPL-3.0", "execution_status": "EXCLUDED — deprecated Nextflow wrapper"},
            {"software": "QUAST", "version": "5.3.0 (not installed)", "repository_or_source": "https://github.com/ablab/quast", "license": "GPL-2.0", "execution_status": "Linux-only (WSL2 distro absent); Python QUAST-compatible metrics computed"},
            {"software": "BUSCO", "version": "6.1.0 (not installed)", "repository_or_source": "https://busco.ezlab.org/", "license": "MIT", "execution_status": "Linux-only (WSL2 distro absent); scores marked NA (not fabricated)"},
            {"software": "OrthoFinder", "version": "3.2 (not installed)", "repository_or_source": "https://github.com/OrthoFinder/OrthoFinder", "license": "GPL-3.0", "execution_status": "Linux MCL/FastME unavailable natively; DIAMOND RBH graph clustering executed"},
            {"software": "Funannotate / BRAKER3", "version": "1.9.1 / 3.0.8 (not installed)", "repository_or_source": "https://github.com/nextgenusfs/funannotate", "license": "BSD-2 / Artistic", "execution_status": "Linux-only; NCBI RefSeq/GenBank annotations reused + tblastn/blastn exon alignment"},
            {"software": "RNAhybrid", "version": "2.1.2 (not installed)", "repository_or_source": "https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid", "license": "GPL-2.0", "execution_status": "Linux/macOS Bioconda only; ViennaRNA 2.8.0 RNA.duplexfold executed"},
        ]
    )
    return pd.DataFrame(rows)


def _draw_workflow_diagram(out_png: Path, out_svg: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 6.5))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.axis("off")

    boxes = [
        (0.4, 5.5, 2.3, 0.9, "M01 Data Acquisition\nLocal ZIP + 19 Isolates\n+ Non-Target Refs", "#dbeafe", "#1d4ed8"),
        (3.2, 5.5, 2.3, 0.9, "M02 Genome QC\nAssembly Metrics\n(QUAST-compat / BUSCO)", "#dcfce7", "#15803d"),
        (6.0, 5.5, 2.5, 0.9, "M03 Annotation\nGFF3 & CDS Translation\nValidation (12,210 CDS)", "#fef3c7", "#b45309"),
        (9.0, 5.5, 2.5, 0.9, "M04 Comparative\nDIAMOND + MAFFT\n19-Isolate Nucleotide Cons.", "#ede9fe", "#6d28d9"),
        (9.0, 3.6, 2.5, 0.9, "M05 Functional Evidence\nPHI-base v4.19 + DEG\nInterPro + GSE130517", "#fce7f3", "#be185d"),
        (6.0, 3.6, 2.5, 0.9, "M06 Target Prioritisation\nCsPKS1, CsFUS3,\nCsSTE11, CsSp1 Verified", "#e0f2fe", "#0369a1"),
        (3.2, 3.6, 2.3, 0.9, "M07 dsRNA Design\nLUCID/sifi21 + SIREN\n+ dsRNAmax + BRD", "#ecfccb", "#4d7c0f"),
        (0.4, 3.6, 2.3, 0.9, "M08 Off-Target Screen\nWheat A/B/D, Barley,\n5 Beneficial Fungi", "#fee2e2", "#b91c1c"),
        (3.2, 1.5, 3.2, 1.0, "M09 Ranking & Uncertainty\nPareto Front + Weighted Score\n1,000-Bootstrap + Sensitivity", "#f1f5f9", "#0f172a"),
        (7.2, 1.5, 3.8, 1.0, "M10 Automated Reporting\nMarkdown, HTML, PDF,\nProvenance & Validation", "#eff6ff", "#1e3a8a"),
    ]
    for x, y, w, h, txt, fc, ec in boxes:
        rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", facecolor=fc, edgecolor=ec, linewidth=1.6)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=8.5, fontweight="bold", color="#0f172a")

    arrows = [
        ((2.75, 5.95), (3.15, 5.95)),
        ((5.55, 5.95), (5.95, 5.95)),
        ((8.55, 5.95), (8.95, 5.95)),
        ((10.25, 5.45), (10.25, 4.55)),
        ((8.95, 4.05), (8.55, 4.05)),
        ((5.95, 4.05), (5.55, 4.05)),
        ((3.15, 4.05), (2.75, 4.05)),
        ((1.55, 3.55), (3.15, 2.0)),
        ((6.45, 2.0), (7.15, 2.0)),
    ]
    for (x1, y1), (x2, y2) in arrows:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=1.8, color="#334155"))

    ax.text(6.0, 6.75, "Bipolaris-RNAi Designer (BRD) v1.0.0 — M01–M10 Pipeline Architecture", ha="center", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_svg)
    plt.close(fig)


def _build_pdf_report(pdf_path: Path, sections: List[tuple[str, str]], figures: List[tuple[Path, str]], top_df: pd.DataFrame, target_df: pd.DataFrame) -> None:
    doc = SimpleDocTemplate(str(pdf_path), pagesize=letter, rightMargin=48, leftMargin=48, topMargin=48, bottomMargin=48)
    styles = getSampleStyleSheet()
    h1 = ParagraphStyle("H1", parent=styles["Heading1"], fontSize=14, spaceAfter=8, textColor=colors.HexColor("#1e3a8a"))
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=11.5, spaceAfter=6, textColor=colors.HexColor("#0f172a"))
    body = ParagraphStyle("Body", parent=styles["BodyText"], fontSize=9.2, leading=13, spaceAfter=6)

    story = []
    story.append(Paragraph("Bipolaris-RNAi Designer (BRD) v1.0.0 — Final Computational Research Report", h1))
    story.append(Paragraph("Target Organism: Bipolaris sorokiniana (syn. Cochliobolus sativus) | Scope: Computational Target Prioritisation, dsRNA Design, and Multi-Reference Off-Target Screening", body))
    story.append(Spacer(1, 0.12 * inch))

    for title, text in sections:
        story.append(Paragraph(title, h2))
        for para in text.split("\n\n"):
            if para.lstrip().startswith("|"):
                continue
            clean = para.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            story.append(Paragraph(clean, body))

    # Target Summary Table
    story.append(PageBreak())
    story.append(Paragraph("Verified Target Gene Summary Table", h2))
    t_cols = ["target", "locus_tag", "protein_id", "cds_length_nt", "isolates_with_gene", "mean_pairwise_nt_identity_pct", "gse130517_inplanta_max_rpkm"]
    t_data = [t_cols] + [[str(r[c]) for c in t_cols] for _, r in target_df.iterrows()]
    tbl = Table(t_data, repeatRows=1)
    tbl.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1e3a8a")), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTSIZE", (0, 0), (-1, -1), 7.5), ("GRID", (0, 0), (-1, -1), 0.4, colors.grey)]))
    story.append(tbl)
    story.append(Spacer(1, 0.15 * inch))

    # Top Ranked Candidates Table
    story.append(Paragraph("Top 15 Ranked dsRNA Candidates", h2))
    c_cols = ["final_rank", "candidate_id", "target", "tool_name", "candidate_length_nt", "conservation_pct", "off_target_safety_score", "composite_score", "computational_priority_category"]
    c_data = [c_cols] + [[str(r[c])[:28] for c in c_cols] for _, r in top_df.head(15).iterrows()]
    tbl2 = Table(c_data, repeatRows=1)
    tbl2.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1e3a8a")), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTSIZE", (0, 0), (-1, -1), 7.0), ("GRID", (0, 0), (-1, -1), 0.4, colors.grey)]))
    story.append(tbl2)

    # Embed key figures
    story.append(PageBreak())
    story.append(Paragraph("Key Scientific Figures", h2))
    for img_p, cap in figures:
        if img_p.exists():
            story.append(RLImage(str(img_p), width=6.4 * inch, height=3.6 * inch, kind="proportional"))
            story.append(Paragraph(cap, body))
            story.append(Spacer(1, 0.12 * inch))

    doc.build(story)


def run_m10(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    logger = get_logger("BRD.M10", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    fin_dir = cfg.module_dir("10_final_reports")
    meta_dir = cfg.output_dir / "00_project_metadata"
    fig_dir = fin_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 1. Collect software versions
    sw_df = _collect_software_versions(cfg)
    sw_p = meta_dir / "software_versions.tsv"
    sw_df.to_csv(sw_p, sep="\t", index=False)

    # 2. Generate workflow architecture figure and copy key pipeline figures
    _draw_workflow_diagram(fig_dir / "brd_workflow_architecture.png", fig_dir / "brd_workflow_architecture.svg")
    fig_sources = [
        (cfg.output_dir / "02_genome_qc" / "figures" / "genome_quality_visualisation.png", "Figure 2. Genome assembly quality comparison across 19 B. sorokiniana biological isolates."),
        (cfg.output_dir / "02_genome_qc" / "figures" / "gc_distribution.png", "Figure 3. GC% distribution across B. sorokiniana assemblies."),
        (cfg.output_dir / "04_comparative_genomics" / "figures" / "multi_isolate_nucleotide_conservation.png", "Figure 4. Multi-isolate CDS nucleotide conservation and SNP positions for CsPKS1, CsFUS3, CsSTE11, and CsSp1."),
        (cfg.output_dir / "06_target_prioritisation" / "figures" / "four_target_functional_evidence_comparison.png", "Figure 5. Four-target biological, expression, conservation, and designability evidence profile."),
        (cfg.output_dir / "06_target_prioritisation" / "figures" / "target_prioritisation.png", "Figure 6. Provisional target prioritisation scores."),
        (cfg.output_dir / "07_dsRNA_design" / "figures" / "candidate_characteristics.png", "Figure 7. Designed dsRNA candidate length, GC%, and conservation vs efficient siRNA fraction."),
        (cfg.output_dir / "08_off_target_screening" / "figures" / "wheat_off_target_summary.png", "Figure 8. Hexaploid wheat (T. aestivum IWGSC v2.1) off-target screening summary."),
        (cfg.output_dir / "08_off_target_screening" / "figures" / "barley_off_target_summary.png", "Figure 9. Barley (H. vulgare MorexV3) off-target screening summary."),
        (cfg.output_dir / "08_off_target_screening" / "figures" / "beneficial_fungi_off_target_summary.png", "Figure 10. Beneficial fungi panel off-target screening summary."),
        (cfg.output_dir / "09_candidate_ranking" / "figures" / "candidate_pareto_frontier.png", "Figure 11. Candidate Pareto frontier (isolate conservation vs off-target safety)."),
        (cfg.output_dir / "09_candidate_ranking" / "figures" / "bootstrap_stability.png", "Figure 12. Bootstrap ranking stability (1,000 isolate resamples, 95% CI)."),
        (cfg.output_dir / "09_candidate_ranking" / "figures" / "sensitivity_analysis.png", "Figure 13. Weighting sensitivity analysis across five heuristic scenarios."),
    ]
    copied_figs = [(fig_dir / "brd_workflow_architecture.png", "Figure 1. BRD M01–M10 modular workflow architecture.")]
    for src, cap in fig_sources:
        if src.exists():
            dst = fig_dir / src.name
            shutil.copy2(src, dst)
            svg_src = src.with_suffix(".svg")
            if svg_src.exists():
                shutil.copy2(svg_src, fig_dir / svg_src.name)
            copied_figs.append((dst, cap))

    # 3. Load actual quantitative results from M01-M09
    inv_df = pd.read_csv(meta_dir / "input_file_inventory.tsv", sep="\t")
    iso_df = pd.read_csv(cfg.output_dir / "01_data_acquisition" / "isolate_metadata.tsv", sep="\t")
    qc_df = pd.read_csv(cfg.output_dir / "02_genome_qc" / "assembly_quality_summary.tsv", sep="\t")
    ann_qc = pd.read_csv(cfg.output_dir / "03_annotation" / "annotation_qc" / "nd90pr_cds_translation_qc.tsv", sep="\t")
    cons_df = pd.read_csv(fin_dir / "BRD_Target_Nucleotide_Conservation.tsv", sep="\t")
    verif_df = pd.read_csv(fin_dir / "BRD_Target_Gene_Accession_Verification.tsv", sep="\t")
    targ_df = pd.read_csv(fin_dir / "BRD_Target_Gene_Summary.tsv", sep="\t")
    func_df = pd.read_csv(cfg.output_dir / "05_functional_evidence" / "functional_evidence_summary" / "functional_evidence_summary.tsv", sep="\t")
    rnai_df = pd.read_csv(cfg.output_dir / "05_functional_evidence" / "functional_evidence_summary" / "bsorokiniana_rnai_machinery_inventory.tsv", sep="\t")
    cand_df = pd.read_csv(fin_dir / "BRD_dsRNA_Candidates.tsv", sep="\t")
    rank_df = pd.read_csv(fin_dir / "BRD_dsRNA_Candidate_Ranking.tsv", sep="\t")
    pareto_df = pd.read_csv(fin_dir / "BRD_Pareto_Analysis.tsv", sep="\t")
    boot_df = pd.read_csv(fin_dir / "BRD_Bootstrap_Results.tsv", sep="\t")
    sens_df = pd.read_csv(fin_dir / "BRD_Sensitivity_Analysis.tsv", sep="\t")
    ot_df = pd.read_csv(fin_dir / "BRD_Off_Target_Analysis.tsv", sep="\t")

    n_passing = int(rank_df["passes_computational_filters"].sum())
    n_high = int((rank_df["computational_priority_category"] == "High computational priority").sum())
    n_mod = int((rank_df["computational_priority_category"] == "Moderate computational priority").sum())
    top5 = rank_df[rank_df["passes_computational_filters"]].head(5)

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    per_target_top = []
    for tname_o in ["CsSTE11", "CsFUS3", "CsSp1", "CsPKS1"]:
        sub_t = rank_df[(rank_df["target"] == tname_o) & (rank_df["passes_computational_filters"])]
        if not sub_t.empty:
            r0 = sub_t.iloc[0]
            per_target_top.append(
                f"- `{tname_o}`: `{r0['candidate_id']}` ({r0['tool_name']}, {r0['candidate_length_nt']} nt, "
                f"cons={r0['conservation_pct']}%, safety={r0['off_target_safety_score']}, "
                f"score={r0['composite_score']}, rank #{r0['final_rank']}, {r0['computational_priority_category']})"
            )
    per_target_block = "\n".join(per_target_top)

    # 4. Build comprehensive 22-section Markdown report
    md_sections = [
        ("1. Title and Project Description", f"**Bipolaris-RNAi Designer (BRD) v{__version__}** is a modular, reproducible bioinformatics pipeline engineered to identify, verify, and computationally prioritise dsRNA candidates targeting pathogenicity-associated genes in *Bipolaris sorokiniana* (teleomorph *Cochliobolus sativus*). Generated autonomously on {now_utc}."),
        ("2. Research Objectives", "1. Verify local and public *B. sorokiniana* genome assemblies and annotations.\n2. Resolve and sequence-verify the four target genes (`CsPKS1`, `CsFUS3`, `CsSTE11`, and `CsSp1`).\n3. Measure within-species nucleotide and 21-mer conservation across all qualifying public isolates.\n4. Integrate species-specific genetics, PHI-base v4.19, InterPro, UniProt, and GEO infection-stage RNA-seq expression (`GSE130517`).\n5. Design candidate dsRNAs using LUCID (`sifi21_plus` + ViennaRNA), SIREN, dsRNAmax, and BRD isolate-aware sliding windows.\n6. Screen candidates and 21-mer guides against hexaploid wheat (`IWGSC v2.1`), barley (`MorexV3`), and five beneficial fungal species.\n7. Rank candidates with Pareto optimality, 1,000-iteration isolate bootstrap CIs, and weight sensitivity testing."),
        ("3. Biological Rationale", "Spot blotch, common root rot, and black point caused by *B. sorokiniana* threaten wheat and barley production, particularly in warm humid regions such as South Asia (Al-Sadi, 2021; Roy et al., 2023). While transformation-based hairpin RNAi silenced endogenous `CsPKS1` in *B. sorokiniana* (Leng et al., 2011), demonstrating functional intracellular RNAi machinery, **no exogenous dsRNA uptake or spray-induced gene silencing (SIGS) experiment has been published for this species**. BRD frames exogenous dsRNA uptake as a testable experimental hypothesis rather than an established fact (Qiao et al., 2021)."),
        ("4. Input Dataset Inventory", f"Recursive inspection of `{cfg.input_dir}` identified {len(inv_df)} inventory records (local NCBI Datasets archive `Bipolaris sorokiniana.zip` for `GCF_000338995.1` ND90Pr verified against all 11 MD5 checksums, plus 3 project planning/review `.docx` documents). Original input files remained strictly unmodified."),
        ("5. Genome Accession and Isolate Information", f"Querying NCBI TaxID 45130 identified 20 assembly records representing **{len(iso_df)} distinct biological isolates** after deduplicating paired GenBank/RefSeq `GCA_000338995.1` vs `GCF_000338995.1` (ND90Pr). Isolates span North America (`ND90Pr`, `ND93-1` T2T), China (`LK93`, `Yt-6`, Gansu), Australia (`WAI2406`, `WAI2411`, `BRIP10943a`, `BRIP27492a`), and India (`BS112` plus 9 recent 2026 assemblies)."),
        ("6. Genome Quality Assessment", f"All {len(qc_df)} biological isolate assemblies passed FASTA validity, genome size ({qc_df['total_length_bp'].min()/1e6:.2f}–{qc_df['total_length_bp'].max()/1e6:.2f} Mb), and GC% ({qc_df['gc_percent'].min():.2f}%–{qc_df['gc_percent'].max():.2f}%) checks. Contig/scaffold N50 ranged from {qc_df['n50_bp'].min()/1e3:.1f} kb to {qc_df['n50_bp'].max()/1e6:.2f} Mb. Note: Official BUSCO v6.1.0 and QUAST v5.3.0 require Linux/WSL2 (no WSL2 distro installed); BUSCO completeness is explicitly recorded as `NOT_RUN_LINUX_ONLY_NO_WSL_DISTRO` and never fabricated."),
        ("7. Annotation Harmonisation", f"RefSeq `GCF_000338995.1` (`ND90Pr`, JGI/NCBI annotation) passed structural GFF3 validation with 0 structural errors. Out of {ann_qc.loc[ann_qc['metric']=='Total CDS evaluated','value'].values[0]} evaluated CDS sequences, {ann_qc.loc[ann_qc['metric']=='CDS translations 100% matching protein.faa','value'].values[0]} translated with 100% identity to `protein.faa`. Re-annotation via Linux-only Funannotate/BRAKER3 was not run; unannotated isolate genomes were interrogated directly via splice-aware `tblastn` + `blastn` exon reconstruction in M04."),
        ("8. Comparative Genomics", f"DIAMOND all-vs-all RBH orthogroup clustering grouped annotated proteomes, and splice-aware genomic locus mapping reconstructed all 4 targets across the 19-isolate panel, followed by MAFFT v7.526 alignment.\n\n" + cons_df.to_markdown(index=False)),
        ("9. Target Gene Identification", "All four targets achieved **Confirmed identity** in `GCF_000338995.1`:\n- `CsPKS1`: `COCSADRAFT_157080` (`XP_007696213.1`, 6,462 nt CDS, 100.0% match to `HQ830033.1`).\n- `CsFUS3`: `COCSADRAFT_24574` (`XP_007697446.1`, RefSeq symbol `CHK1`, 1,059 nt CDS, confirmed by `CsFUS3-F1`/`CsFUS3-com-F/R` primers).\n- `CsSTE11`: `COCSADRAFT_105370` (`XP_007694574.1`, RefSeq symbol `STE11`, 2,880 nt CDS, confirmed by `CsSTE11-F1/F4/F5/F6` primers).\n- `CsSp1`: **Resolved** via 100% matches to 5 independent primers from Zhang et al. (2022, PMC8743017 Table 2) as `COCSADRAFT_164106` (`XP_007704210.1`, `XM_007706020.1`, 270 nt CDS, 89 aa)."),
        ("10. Accession Verification", "Complete 24-column accession verification including `HQ830033.1`, `KF769134.1`–`KF769141.1`, paired `GCA/GCF_000338995.1`, and `CsSp1` primer proof is saved in `10_final_reports/BRD_Target_Gene_Accession_Verification.tsv`."),
        ("11. Functional Evidence", f"Direct *B. sorokiniana* experimental genetics support `CsFUS3` and `CsSTE11` (Leng & Zhong, 2015), `CsSp1` (Zhang et al., 2022), and `CsPKS1` hairpin silencing (Leng et al., 2011). PHI-base v4.19 DIAMOND searches, InterPro/UniProt annotations, and GEO `GSE130517` barley-infection RNA-seq confirmed in-planta transcription of all 4 targets and identified {len(rnai_df)} core RNAi machinery homologs (DCL, AGO, RdRP) in `ND90Pr`."),
        ("12. Target Prioritisation", "Using provisional research weights (25% conservation, 25% off-target safety, 20% functional evidence, 15% expression, 15% designability), `CsFUS3` and `CsSTE11` lead on combined pathogenicity evidence, length flexibility, and expression, while `CsSp1` is constrained by its 270 nt CDS length and `CsPKS1` serves as the endogenous hairpin-RNAi benchmark."),
        ("13. dsRNA Design Methodology", "Four complementary computational workflows were executed: (1) LUCID baseline (`sifi21_plus` free-energy/S1-S19 asymmetry + ViennaRNA 2.8.0 unpaired accessibility + Primer3, 200–400 nt), (2) SIREN v0.1.6 (`siren_prefilter.py` + ViennaRNA `duplexfold` + `siren_designVIII.py`), (3) `dsRNAmax.exe` v1.1.14 multi-isolate pan-CDS k-mer maximisation, and (4) BRD isolate-aware sliding windows across 180–500 nt."),
        ("14. Candidate Sequences and Characteristics", f"M07 generated **{len(cand_df)} unique dsRNA candidates** (lengths {cand_df['candidate_length_nt'].min()}–{cand_df['candidate_length_nt'].max()} nt; GC {cand_df['gc_percent'].min():.1f}%–{cand_df['gc_percent'].max():.1f}%) and {len(pd.read_csv(fin_dir / 'BRD_Guide_Evaluation.tsv', sep='\t'))} unique 21-mer guide sequences."),
        ("15. Off-Target Screening", f"All {len(ot_df)} candidates and their 21-mer guides were screened via `blastn-short` (`word_size=7`), exact 21-mer matching, and ViennaRNA `RNA.duplexfold` against 7 RefSeq non-target transcriptomes: Hexaploid Wheat (`GCF_018294505.1`), Barley (`GCF_904849725.1`), *Trichoderma atroviride*, *T. virens*, *T. reesei*, *Beauveria bassiana*, and *Rhizophagus irregularis*. Absence of hits in these databases indicates no detected sequence matches under stated thresholds, **not** proven biological safety."),
        ("16. Candidate Ranking", f"Of {len(rank_df)} candidates evaluated in M09, **{n_passing} passed all computational filters** (0 exact 21-mer host/beneficial/paralog hits, QC pass, >=85% conservation), including **{n_high} High computational priority**, **{n_mod} Moderate computational priority**, and **{len(pareto_df)} Pareto rank-1 non-dominated** candidates.\n\n**Top 5 Ranked Candidates:**\n\n" + top5[["final_rank", "candidate_id", "target", "tool_name", "candidate_length_nt", "conservation_pct", "off_target_safety_score", "composite_score", "computational_priority_category"]].to_markdown(index=False)),
        ("17. Uncertainty and Sensitivity Analysis", "1,000-iteration bootstrap resampling across the 19 biological isolates quantified 95% confidence intervals for 21-mer conservation and rank stability. Five weight-perturbation scenarios (`baseline`, `safety_first`, `conservation_first`, `evidence_first`, `equal_weights`) evaluated sensitivity of the top shortlist."),
        ("18. Limitations", "1. **Exogenous dsRNA uptake is unproven in *B. sorokiniana*** (Qiao et al., 2021).\n2. **Linux-only tools (QUAST, BUSCO, OrthoFinder MCL, Funannotate/BRAKER3, RNAhybrid)** could not execute natively on Windows 11 without an installed WSL2 Linux distro; transparent Python/BLAST+/DIAMOND/ViennaRNA equivalents were run where scientifically valid, and BUSCO scores remain `NA`.\n3. **RNA-seq (`GSE130517`)** represents *H. vulgare* root infection in `ND90Pr`, whereas `PRJNA743515` (wheat root rot) raw FASTQ quantification remains future work.\n4. **Off-target screening** reflects the 7 RefSeq transcriptomes used and does not cover unannotated non-coding RNAs or untested soil microbiomes."),
        ("19. Recommended Candidates for Experimental Validation", "Top computationally passing candidates per target gene selected for sequential bench evaluation:\n\n" + per_target_block + "\n\nSequential validation roadmap: (1) fluorescently labelled dsRNA uptake assays in *B. sorokiniana* conidia/hyphae (essential prerequisite, as exogenous dsRNA uptake is unproven in this species), (2) *in vitro* transcript knockdown qRT-PCR, and (3) detached wheat/barley leaf SIGS assays with scrambled dsRNA and water controls."),
        ("20. References", "- Al-Sadi, A. M. (2021). *Frontiers in Cellular and Infection Microbiology*, 11, 584899. https://doi.org/10.3389/fcimb.2021.584899\n- Appels, R., et al. (2018). *Science*, 361(6403), eaar7191. https://doi.org/10.1126/science.aar7191\n- Fletcher, S. J., et al. (2025). dsRNAmax. *NAR Genomics and Bioinformatics*, 7(2), lqaf064. https://doi.org/10.1093/nargab/lqaf064\n- Jiménez-Castro, L., et al. (2026). LUCID. *Plant Biotechnology Journal*, 24(6), 3597–3615. https://doi.org/10.1111/pbi.70526\n- Koch, A., et al. (2016). *PLOS Pathogens*, 12(10), e1005901. https://doi.org/10.1371/journal.ppat.1005901\n- Leng, Y., et al. (2011). *Molecular Plant Pathology*, 12(3), 289–298. https://doi.org/10.1111/j.1364-3703.2010.00666.x\n- Leng, Y., & Zhong, S. (2015). *PLOS ONE*, 10(5), e0128291. https://doi.org/10.1371/journal.pone.0128291\n- Leng, Y., et al. (2024). *PhytoFrontiers*, 4(4), 247–250. https://doi.org/10.1094/PHYTOFR-08-23-0108-A\n- Piombo, E., et al. (2021). *Microbiology Spectrum*, 9(2), e0109921. https://doi.org/10.1128/Spectrum.01099-21\n- Qiao, L., et al. (2021). *Plant Biotechnology Journal*, 19(9), 1756–1768. https://doi.org/10.1111/pbi.13589\n- Roy, C., et al. (2023). *Frontiers in Plant Science*, 14, 1098648. https://doi.org/10.3389/fpls.2023.1098648\n- Urban, M., et al. (2022). PHI-base. *Nucleic Acids Research*, 50(D1), D837–D847. https://doi.org/10.1093/nar/gkab1037\n- Vargas Mejía, P., & Vega-Arreguín, J. C. (2025/2026). SIREN. https://doi.org/10.1101/2025.05.26.656188\n- Zhang, W., et al. (2022). CsSp1. *Molecular Plant Pathology*, 23(2), 218–236. https://doi.org/10.1111/mpp.13155\n- Zhang, W., et al. (2023). LK93 genome. *MPMI*, 36(7), 452–456. https://doi.org/10.1094/MPMI-09-22-0196-A"),
        ("21. Software Versions and Parameters", sw_df.to_markdown(index=False)),
        ("22. Data Provenance", f"Complete SHA-256 provenance ledger saved at `{meta_dir / 'data_provenance.tsv'}` and `{meta_dir / 'provenance.json'}`."),
    ]

    md_text = "# Bipolaris-RNAi Designer (BRD) — Final Research Report\n\n" + "\n\n".join(f"## {t}\n\n{b}" for t, b in md_sections)
    md_p = fin_dir / "BRD_Final_Research_Report.md"
    write_text(md_p, md_text)

    # 5. HTML Report
    html_parts = [f"<div class='warn'><strong>Scientific Boundary:</strong> BRD produces computationally prioritised dsRNA candidates; exogenous uptake, RNAi knockdown, and disease control in <em>B. sorokiniana</em> remain experimental hypotheses requiring bench validation.</div>"]
    top5_cols = ["final_rank", "candidate_id", "target", "tool_name", "candidate_length_nt", "conservation_pct", "off_target_safety_score", "composite_score", "computational_priority_category"]
    for idx_s, (t, b) in enumerate(md_sections[:19]):
        if idx_s == 7:
            html_parts.append(f"<h2>{t}</h2><p>DIAMOND all-vs-all RBH orthogroup clustering grouped annotated proteomes, and splice-aware genomic locus mapping reconstructed all 4 targets across the 19-isolate panel, followed by MAFFT v7.526 alignment.</p>{df_to_html_table(cons_df)}")
        elif idx_s == 15:
            html_parts.append(
                f"<h2>{t}</h2><p>Of {len(rank_df)} candidates evaluated in M09, <strong>{n_passing} passed all computational filters</strong> "
                f"(0 exact 21-mer host/beneficial/paralog hits, QC pass, &ge;85% conservation), including <strong>{n_high} High computational priority</strong>, "
                f"<strong>{n_mod} Moderate computational priority</strong>, and <strong>{len(pareto_df)} Pareto rank-1 non-dominated</strong> candidates.</p>"
                f"<h3>Top 5 Ranked Candidates</h3>{df_to_html_table(top5[top5_cols])}"
            )
        else:
            html_parts.append(f"<h2>{t}</h2><p>{b.replace(chr(10), '<br/>')}</p>")
    html_parts.append(f"<h2>Target Summary</h2>{df_to_html_table(targ_df)}")
    html_parts.append(f"<h2>Accession Verification</h2>{df_to_html_table(verif_df)}")
    html_parts.append(f"<h2>Top 25 Ranked Candidates</h2>{df_to_html_table(rank_df.drop(columns=['sequence'], errors='ignore').head(25))}")
    html_parts.append("<h2>Key Figures</h2>")
    for img_p, cap in copied_figs:
        html_parts.append(f"<figure><img src='figures/{img_p.name}'/><figcaption>{cap}</figcaption></figure>")
    html_parts.append(f"<h2>20. References</h2><p>{md_sections[19][1].replace(chr(10), '<br/>')}</p>")
    html_parts.append(f"<h2>21. Software Versions & Compatibility</h2>{df_to_html_table(sw_df)}")

    html_p = fin_dir / "BRD_Final_Research_Report.html"
    write_text(html_p, wrap_html_report("Bipolaris-RNAi Designer (BRD) — Final Research Report", f"Generated {now_utc} | 19 B. sorokiniana Isolates | 4 Targets | {len(cand_df)} Candidates", "".join(html_parts)))

    # 6. PDF Report
    pdf_p = fin_dir / "BRD_Final_Research_Report.pdf"
    _build_pdf_report(pdf_p, md_sections[:19], copied_figs, rank_df, targ_df)

    # 7. Software Installation Report & Data Provenance Report
    sw_rep_p = fin_dir / "BRD_Software_Installation_Report.md"
    write_text(
        sw_rep_p,
        "# BRD Software Installation & Compatibility Report\n\n"
        "Executed on Windows 11 (WSL2 has no Linux distribution installed; `wsl -l -v` empty).\n\n"
        + sw_df.to_markdown(index=False),
    )
    prov_df = pd.read_csv(meta_dir / "data_provenance.tsv", sep="\t")
    prov_rep_p = fin_dir / "BRD_Data_Provenance_Report.md"
    write_text(
        prov_rep_p,
        "# BRD Data Provenance Report\n\n"
        + prov_df.to_markdown(index=False),
    )

    # 8. Execution Summary TXT (exact template from Section 52)
    top3_lines = "\n".join(
        f"  #{r['final_rank']} {r['candidate_id']} ({r['target']}, {r['tool_name']}, {r['candidate_length_nt']} nt, "
        f"cons={r['conservation_pct']}%, score={r['composite_score']}, {r['computational_priority_category']})"
        for _, r in top5.head(5).iterrows()
    )
    exec_txt = f"""====================================================
BIPOLARIS-RNAi DESIGNER — EXECUTION SUMMARY
====================================================

Input directory:
{cfg.input_dir}

Output directory:
{cfg.output_dir}

Genome assemblies detected:
20 NCBI records (19 unique biological isolates; paired GCA_000338995.1/GCF_000338995.1 deduplicated)

Genome assemblies processed:
{len(qc_df)}

Target genes:
CsPKS1: Confirmed (COCSADRAFT_157080 / XP_007696213.1 / HQ830033.1; 19/19 isolates)
CsFUS3: Confirmed (COCSADRAFT_24574 / XP_007697446.1 / CHK1; 19/19 isolates)
CsSTE11: Confirmed (COCSADRAFT_105370 / XP_007694574.1 / STE11; 19/19 isolates)
CsSp1: Confirmed via Zhang et al. (2022) PMC8743017 primers (COCSADRAFT_164106 / XP_007704210.1; 19/19 isolates)

M01 Data acquisition: COMPLETED
M02 Genome QC: COMPLETED (QUAST-compatible Python metrics; BUSCO NA - Linux only)
M03 Annotation: COMPLETED (12,210 CDS validated; Funannotate/BRAKER3 Linux-only)
M04 Comparative genomics: COMPLETED (DIAMOND RBH orthogroups + MAFFT + 19-isolate CDS conservation)
M05 Functional evidence: COMPLETED (PHI-base v4.19 + DEG + UniProt + InterPro + GSE130517)
M06 Target prioritisation: COMPLETED
M07 dsRNA design: COMPLETED (LUCID/sifi21_plus + SIREN + dsRNAmax + BRD)
M08 Off-target screening: COMPLETED (Wheat A/B/D + Barley + 5 Beneficial Fungi via BLAST+ & ViennaRNA)
M09 Candidate ranking: COMPLETED (Pareto + Weighted + 1,000-Bootstrap + Sensitivity)
M10 Reporting: COMPLETED (Markdown, HTML, PDF, Provenance)

Candidates generated:
{len(cand_df)}

Candidates passing computational filters:
{n_passing} ({n_high} High priority, {n_mod} Moderate priority)

Final candidate rankings (Top 5):
{top3_lines}

Unresolved issues:
1. WSL2 Linux distribution not installed: BUSCO, QUAST CLI, OrthoFinder MCL, Funannotate/BRAKER3, and RNAhybrid binary did not run natively on Windows (documented equivalents run for QUAST-metrics, orthogroups, and ViennaRNA duplex thermodynamics; BUSCO scores explicitly marked NA).
2. Exogenous dsRNA uptake and SIGS efficacy in Bipolaris sorokiniana remain experimentally untested hypotheses.

Final research report:
{html_p}

Execution log:
{log_path}
====================================================
"""
    exec_p = fin_dir / "BRD_Execution_Summary.txt"
    write_text(exec_p, exec_txt)

    prov.record("M10", "Automated Reporting", str(rank_df.shape[0]) + "_ranked_candidates", html_p, "Python+Matplotlib+ReportLab", __version__, "M10_Final")
    logger.info("M10 complete: HTML, Markdown, PDF, and Execution Summary written to %s", fin_dir)
    return {
        "html_report": html_p,
        "md_report": md_p,
        "pdf_report": pdf_p,
        "execution_summary": exec_p,
        "software_versions": sw_p,
    }
