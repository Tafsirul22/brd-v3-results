"""Module M05 — Functional Evidence Integration (NCBI, PHI-base, DEG, InterPro, UniProt, Expression, RNAi Machinery)."""

from __future__ import annotations

import gzip
import io
import json
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd
import requests
from Bio import SeqIO

from brd.config import BRDConfig
from brd.integrations.diamond import build_diamond_db, run_diamond_blastp
from brd.integrations.ncbi import NCBIClient
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker


def _fetch_uniprot_and_interpro(locus_tags: Dict[str, str], logger) -> tuple[pd.DataFrame, pd.DataFrame]:
    uniprot_rows = []
    interpro_rows = []
    for tname, locus in locus_tags.items():
        url = f"https://rest.uniprot.org/uniprotkb/search?query={locus}&format=json&size=5"
        try:
            r = requests.get(url, timeout=20)
            r.raise_for_status()
            results = r.json().get("results", [])
        except Exception as exc:
            logger.warning("UniProt query failed for %s (%s): %s", tname, locus, exc)
            results = []

        if not results:
            uniprot_rows.append({"target": tname, "locus_tag": locus, "uniprot_accession": "Not found", "protein_name": "", "organism": "", "length": 0, "go_terms": ""})
            continue

        # Prefer B. sorokiniana / C. sativus match
        entry = results[0]
        for cand in results:
            org = cand.get("organism", {}).get("scientificName", "")
            if "bipolaris" in org.lower() or "cochliobolus" in org.lower():
                entry = cand
                break

        uacc = entry.get("primaryAccession", "")
        pname = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        if not pname:
            sub = entry.get("proteinDescription", {}).get("submissionNames", [])
            if sub:
                pname = sub[0].get("fullName", {}).get("value", "")
        org_name = entry.get("organism", {}).get("scientificName", "")
        seq_len = entry.get("sequence", {}).get("length", 0)
        go_terms = []
        for xref in entry.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "GO":
                props = {p["key"]: p["value"] for p in xref.get("properties", [])}
                go_terms.append(f"{xref.get('id')}:{props.get('GoTerm', '')}")

        uniprot_rows.append(
            {
                "target": tname,
                "locus_tag": locus,
                "uniprot_accession": uacc,
                "protein_name": pname,
                "organism": org_name,
                "length": seq_len,
                "go_terms": "; ".join(go_terms[:10]),
            }
        )

        # Query InterPro REST API for this UniProt accession
        ipr_url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{uacc}?page_size=50"
        try:
            ir = requests.get(ipr_url, timeout=20)
            if ir.status_code == 200:
                for item in ir.json().get("results", []):
                    meta = item.get("metadata", {})
                    locs = []
                    for prot in item.get("proteins", []):
                        for epl in prot.get("entry_protein_locations", []):
                            for frag in epl.get("fragments", []):
                                locs.append(f"{frag.get('start')}-{frag.get('end')}")
                    interpro_rows.append(
                        {
                            "target": tname,
                            "uniprot_accession": uacc,
                            "interpro_accession": meta.get("accession", ""),
                            "name": meta.get("name", ""),
                            "type": meta.get("type", ""),
                            "source_database": meta.get("source_database", ""),
                            "locations": ";".join(locs),
                        }
                    )
        except Exception as exc:
            logger.warning("InterPro query failed for %s (%s): %s", tname, uacc, exc)

    return pd.DataFrame(uniprot_rows), pd.DataFrame(interpro_rows)


def _fetch_gse130517_sample_map() -> Dict[str, str]:
    """Fetch GEO series matrix metadata for GSE130517 to map sample prefixes to biological conditions."""
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130517/matrix/GSE130517_series_matrix.txt.gz"
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        txt = gzip.decompress(r.content).decode("utf-8", errors="replace")
        titles = []
        chars = []
        for line in txt.splitlines():
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "treatment" in line.lower():
                chars = [x.strip('"') for x in line.split("\t")[1:]]
        return {t: (chars[i] if i < len(chars) else t) for i, t in enumerate(titles)}
    except Exception:
        return {}


def run_m05(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    cmd_log = cfg.output_dir / "logs" / "commands.jsonl"
    logger = get_logger("BRD.M05", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    m05_dir = cfg.module_dir("05_functional_evidence")
    ncbi_dir = m05_dir / "ncbi"
    phi_dir = m05_dir / "phi_base"
    ipr_dir = m05_dir / "interpro"
    uni_dir = m05_dir / "uniprot"
    expr_dir = m05_dir / "expression"
    sum_dir = m05_dir / "functional_evidence_summary"
    for d in (ncbi_dir, phi_dir, ipr_dir, uni_dir, expr_dir, sum_dir):
        d.mkdir(parents=True, exist_ok=True)

    ncbi = NCBIClient(email=cfg.ncbi_email, api_key=cfg.ncbi_api_key, log_path=cmd_log)
    diamond_bin = cfg.tool_path("diamond")

    # Load ND90Pr proteins and GFF3
    nd90_prot_faa = cfg.output_dir / "03_annotation" / "proteins" / "genome_proteins.fasta"
    nd90_gff = cfg.output_dir / "03_annotation" / "harmonised" / "harmonised_annotation.gff3"
    prot_seqs = {r.id: r for r in SeqIO.parse(str(nd90_prot_faa), "fasta")}

    targets_cfg = cfg.raw["targets"]
    locus_map = {t: info["nd90pr_locus_tag"] for t, info in targets_cfg.items()}
    pid_to_target = {info["nd90pr_protein_id"]: t for t, info in targets_cfg.items()}

    # Write 4-target protein FASTA
    target_prot_fa = m05_dir / "target_proteins.faa"
    SeqIO.write([prot_seqs[info["nd90pr_protein_id"]] for info in targets_cfg.values()], str(target_prot_fa), "fasta")

    # 1. NCBI Gene & Protein & Conserved Domain Evidence
    gene_ids = [str(info["nd90pr_gene_id"]) for info in targets_cfg.values()]
    prot_ids = [str(info["nd90pr_protein_id"]) for info in targets_cfg.values()]

    gene_sums = ncbi.esummary("gene", gene_ids)
    doc_list = gene_sums.get("DocumentSummarySet", {}).get("DocumentSummary", []) if isinstance(gene_sums, dict) else list(gene_sums)
    ncbi_gene_rows = []
    for tname, info in targets_cfg.items():
        gid = str(info["nd90pr_gene_id"])
        matched = {}
        for d in doc_list:
            if isinstance(d, dict):
                uid = str(getattr(d, "attributes", {}).get("uid", d.get("Id", "")))
                if uid == gid or info["nd90pr_locus_tag"] in str(d):
                    matched = d
                    break
        ncbi_gene_rows.append(
            {
                "target": tname,
                "gene_id": gid,
                "locus_tag": info["nd90pr_locus_tag"],
                "protein_id": info["nd90pr_protein_id"],
                "ncbi_name": matched.get("Name", "") if isinstance(matched, dict) else "",
                "ncbi_description": matched.get("Description", info["description"]) if isinstance(matched, dict) else info["description"],
                "organism": "Bipolaris sorokiniana ND90Pr",
            }
        )
    ncbi_gene_p = ncbi_dir / "ncbi_gene_evidence.tsv"
    pd.DataFrame(ncbi_gene_rows).to_csv(ncbi_gene_p, sep="\t", index=False)

    # Fetch GenPept records for conserved domains & protein metadata
    gp_text = ncbi.efetch_text("protein", prot_ids, rettype="gp", retmode="text")
    (ncbi_dir / "target_proteins.gp").write_text(gp_text, encoding="utf-8")
    gp_recs = list(SeqIO.parse(io.StringIO(gp_text), "genbank"))
    ncbi_prot_rows = []
    cdd_rows = []
    for rec in gp_recs:
        tname = pid_to_target.get(rec.id, rec.id)
        ncbi_prot_rows.append(
            {
                "target": tname,
                "protein_accession": rec.id,
                "length_aa": len(rec.seq),
                "description": rec.description,
                "organism": rec.annotations.get("organism", ""),
            }
        )
        for feat in rec.features:
            if feat.type in ("Region", "Site"):
                cdd_rows.append(
                    {
                        "target": tname,
                        "protein_accession": rec.id,
                        "feature_type": feat.type,
                        "region_name": ",".join(feat.qualifiers.get("region_name", [])),
                        "note": ",".join(feat.qualifiers.get("note", [])),
                        "db_xref": ",".join(feat.qualifiers.get("db_xref", [])),
                        "location": str(feat.location),
                    }
                )
    ncbi_prot_p = ncbi_dir / "ncbi_protein_evidence.tsv"
    cdd_p = ncbi_dir / "conserved_domains.tsv"
    pd.DataFrame(ncbi_prot_rows).to_csv(ncbi_prot_p, sep="\t", index=False)
    pd.DataFrame(cdd_rows).to_csv(cdd_p, sep="\t", index=False)

    # 2. PHI-base v4.19 DIAMOND search (targets + full ND90Pr proteome)
    phi_fas = phi_dir / "phi-base_current.fas"
    phi_csv = phi_dir / "phi-base_current.csv"
    phi_db = phi_dir / "phi_base_db"
    build_diamond_db(diamond_bin, phi_fas, phi_db, log_path=cmd_log)
    phi_hits_raw = run_diamond_blastp(
        diamond_bin, target_prot_fa, phi_db, phi_dir / "targets_vs_phibase_diamond.tsv", evalue=1e-5, max_target_seqs=10, log_path=cmd_log
    )
    # Also run full proteome vs PHI-base for LUCID CEP/CNAP support
    run_diamond_blastp(
        diamond_bin, nd90_prot_faa, phi_db, phi_dir / "nd90pr_proteome_vs_phibase_diamond.tsv", evalue=1e-5, max_target_seqs=3, log_path=cmd_log
    )

    # Parse PHI-base FASTA headers & CSV annotations
    phi_meta = {}
    if phi_csv.exists():
        try:
            df_phic = pd.read_csv(phi_csv, low_memory=False)
            if "PHI_MolConn_ID" in df_phic.columns or "PHI-base accession" in df_phic.columns:
                id_col = "PHI-base accession" if "PHI-base accession" in df_phic.columns else df_phic.columns[0]
                for _, pr in df_phic.iterrows():
                    phi_meta[str(pr[id_col])] = pr.to_dict()
        except Exception as exc:
            logger.warning("Could not parse phi-base_current.csv fully: %s", exc)

    phi_rows = []
    for _, hr in phi_hits_raw.iterrows():
        q = str(hr["qseqid"])
        tname = pid_to_target.get(q, q)
        s = str(hr["sseqid"])
        # PHI-base FASTA header format: >UniProt#PHI:123#Gene#Taxon#Species#Phenotype
        parts = s.split("#")
        phi_id = parts[1] if len(parts) > 1 else s
        phi_gene = parts[2] if len(parts) > 2 else ""
        phi_sp = parts[4].replace("_", " ") if len(parts) > 4 else ""
        phi_pheno = parts[5].replace("__", "; ").replace("_", " ") if len(parts) > 5 else ""
        csv_row = phi_meta.get(phi_id, {})
        phi_rows.append(
            {
                "target": tname,
                "query_protein": q,
                "phi_subject_header": s,
                "phi_id": phi_id,
                "phi_gene": phi_gene,
                "phi_species": phi_sp,
                "mutant_phenotype": phi_pheno or csv_row.get("Mutant Phenotype", ""),
                "pident": hr["pident"],
                "qcovhsp": hr["qcovhsp"],
                "evalue": hr["evalue"],
                "bitscore": hr["bitscore"],
                "evidence_class": "Direct B. sorokiniana" if "bipolaris" in phi_sp.lower() or "cochliobolus sativus" in phi_sp.lower() else "Homolog in related fungus (PHI-base)",
            }
        )
    phi_p = phi_dir / "phi_base_evidence.tsv"
    pd.DataFrame(phi_rows).to_csv(phi_p, sep="\t", index=False)

    # 3. LUCID Essential Protein Database (DEG) search
    deg_fa = Path(cfg.tool_path("lucid_repo")) / "data" / "references" / "Essential_protein_database.fasta"
    deg_hits = pd.DataFrame()
    if deg_fa.exists():
        deg_db = phi_dir / "lucid_deg_db"
        build_diamond_db(diamond_bin, deg_fa, deg_db, log_path=cmd_log)
        deg_hits = run_diamond_blastp(
            diamond_bin, target_prot_fa, deg_db, phi_dir / "targets_vs_deg_diamond.tsv", evalue=1e-5, max_target_seqs=5, log_path=cmd_log
        )
        run_diamond_blastp(
            diamond_bin, nd90_prot_faa, deg_db, phi_dir / "nd90pr_proteome_vs_deg_diamond.tsv", evalue=1e-5, max_target_seqs=1, log_path=cmd_log
        )

    # 4. UniProt & InterPro annotations
    uni_df, ipr_df = _fetch_uniprot_and_interpro(locus_map, logger)
    uni_p = uni_dir / "uniprot_annotations.tsv"
    ipr_p = ipr_dir / "interpro_annotations.tsv"
    uni_df.to_csv(uni_p, sep="\t", index=False)
    ipr_df.to_csv(ipr_p, sep="\t", index=False)

    # 5. Infection-Stage Expression (GSE130517) + PRJNA743515 SRA Metadata
    gse_fp = expr_dir / "GSE130517_Bipolaris_sorokiniana_blast2GOannotation_rawrc_rpkm_rei.txt.gz"
    gse_map_fp = expr_dir / "GSE130517_Bipolaris_sorokiniana_JGI_geneID.txt.gz"
    if not gse_map_fp.exists():
        map_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130517/suppl/GSE130517_Bipolaris_sorokiniana_JGI_geneID.txt.gz"
        rm = requests.get(map_url, timeout=60)
        rm.raise_for_status()
        gse_map_fp.write_bytes(rm.content)

    bipso_to_cocsa: Dict[str, str] = {}
    if gse_map_fp.exists():
        with gzip.open(gse_map_fp, "rt", encoding="utf-8", errors="replace") as fh:
            df_jgi_map = pd.read_csv(fh, sep="\t")
        for _, mr in df_jgi_map.iterrows():
            m_id = re.search(r"jgi\|Cocsa1\|(\d+)\|", str(mr.get("jgi_protein_id", "")))
            if m_id:
                bipso_to_cocsa[str(mr["gene_id"])] = f"COCSADRAFT_{m_id.group(1)}"

    sample_map = _fetch_gse130517_sample_map()
    expr_rows = []
    rnai_expr_map: Dict[str, dict] = {}
    if gse_fp.exists():
        with gzip.open(gse_fp, "rt", encoding="utf-8", errors="replace") as fh:
            df_gse = pd.read_csv(fh, sep="\t", low_memory=False)
        rpkm_cols = [c for c in df_gse.columns if c.endswith("_rpkm") and not c.endswith("_rpkm_rei")]
        bs_ctrl_cols = [c for c in rpkm_cols if c.startswith("Bs_")]
        inplanta_cols = [c for c in rpkm_cols if re.match(r"^[1-5][AB]_", c)]
        svbs_cols = [c for c in rpkm_cols if c.startswith("SvBs_")]

        for _, gr in df_gse.iterrows():
            g_id = str(gr["gene"])
            desc = str(gr.get("seq_desc", ""))
            loc_key = bipso_to_cocsa.get(g_id, "")
            if not loc_key:
                m_loc = re.search(r"(COCSADRAFT_\d+)", desc)
                loc_key = m_loc.group(1) if m_loc else ""
            bs_mean = float(gr[bs_ctrl_cols].mean()) if bs_ctrl_cols else 0.0
            in_mean = float(gr[inplanta_cols].mean()) if inplanta_cols else 0.0
            in_max = float(gr[inplanta_cols].max()) if inplanta_cols else 0.0
            svbs_mean = float(gr[svbs_cols].mean()) if svbs_cols else 0.0
            if loc_key:
                rnai_expr_map[loc_key] = {"bs_ctrl_rpkm": round(bs_mean, 2), "inplanta_mean_rpkm": round(in_mean, 2), "inplanta_max_rpkm": round(in_max, 2)}

            for tname, locus in locus_map.items():
                if loc_key == locus or locus in desc:
                    raw_cols = [c.replace("_rpkm", "") for c in inplanta_cols]
                    in_raw_sum = int(gr[raw_cols].sum()) if raw_cols else 0
                    bs_raw_sum = int(gr[[c.replace("_rpkm", "") for c in bs_ctrl_cols]].sum()) if bs_ctrl_cols else 0
                    expr_rows.append(
                        {
                            "target": tname,
                            "locus_tag": locus,
                            "gse_gene_id": g_id,
                            "dataset": "GSE130517 (Hordeum vulgare-B. sorokiniana ND90Pr root interaction)",
                            "bs_axenic_mean_rpkm": round(bs_mean, 2),
                            "inplanta_mean_rpkm": round(in_mean, 2),
                            "inplanta_max_rpkm": round(in_max, 2),
                            "svbs_co_inoc_mean_rpkm": round(svbs_mean, 2),
                            "inplanta_raw_reads_sum": in_raw_sum,
                            "axenic_raw_reads_sum": bs_raw_sum,
                            "infection_expression_detected": in_max > 1.0,
                            "notes": "Measured RPKM from GSE130517 supplementary matrix mapped via JGI protein ID.",
                        }
                    )

    # Also record PRJNA743515 SRA metadata (Zhang et al. 2022 wheat root rot)
    try:
        sra_res = ncbi.esearch("sra", "PRJNA743515[BioProject]", retmax=20)
        if sra_res.get("IdList"):
            runinfo = ncbi.efetch_text("sra", sra_res["IdList"], rettype="runinfo", retmode="text")
            (expr_dir / "PRJNA743515_sra_runinfo.csv").write_text(runinfo, encoding="utf-8")
    except Exception as exc:
        logger.warning("Could not fetch PRJNA743515 SRA runinfo: %s", exc)

    expr_df = pd.DataFrame(expr_rows)
    expr_p = expr_dir / "infection_expression.tsv"
    expr_df.to_csv(expr_p, sep="\t", index=False)

    # 6. Annotate B. sorokiniana ND90Pr Core RNAi Machinery (DCL, AGO, RdRP)
    rnai_rows = []
    rnai_patterns = {
        "Dicer-like (DCL)": ["IPR000999", "dicer"],
        "Argonaute (AGO/PAZ/Piwi)": ["IPR003100", "IPR003165", "argonaute", "piwi"],
        "RdRP (QDE-1 homolog)": ["IPR007855", "rna-dependent rna polymerase"],
    }
    with open(nd90_gff, "r", encoding="utf-8", errors="replace") as fh:
        seen_rnai = set()
        for line in fh:
            if "\tCDS\t" not in line:
                continue
            low = line.lower()
            for comp, pats in rnai_patterns.items():
                if any(p.lower() in low for p in pats):
                    m_name = re.search(r"Name=([^;]+)", line)
                    m_note = re.search(r"Note=([^;]+)", line)
                    m_jgi = re.search(r"JGIDB:Cocsa1_(\d+)", line)
                    m_ipr = re.findall(r"(IPR\d+)", line)
                    pid = m_name.group(1) if m_name else ""
                    if pid and pid not in seen_rnai:
                        seen_rnai.add(pid)
                        ltag = f"COCSADRAFT_{m_jgi.group(1)}" if m_jgi else ""
                        ex_i = rnai_expr_map.get(ltag, {})
                        rnai_rows.append(
                            {
                                "rnai_component": comp,
                                "protein_id": pid,
                                "locus_tag": ltag,
                                "protein_length_aa": len(prot_seqs[pid].seq) if pid in prot_seqs else 0,
                                "interpro_domains": ";".join(sorted(set(m_ipr))),
                                "gff_note": m_note.group(1) if m_note else "",
                                "gse130517_axenic_rpkm": ex_i.get("bs_ctrl_rpkm", "NA"),
                                "gse130517_inplanta_max_rpkm": ex_i.get("inplanta_max_rpkm", "NA"),
                                "evidence_status": "Computational domain/GFF annotation in ND90Pr (not functionally knocked out in B. sorokiniana)",
                            }
                        )
    rnai_df = pd.DataFrame(rnai_rows)
    rnai_df.to_csv(sum_dir / "bsorokiniana_rnai_machinery_inventory.tsv", sep="\t", index=False)

    # 7. Consolidated Functional Evidence Summary
    lit_evidence = {
        "CsPKS1": {
            "direct_exp": "Direct B. sorokiniana experimental evidence (hairpin RNAi silencing -> albino/reduced-melanisation phenotypes, Leng et al., 2011; NOT exogenous SIGS)",
            "bio_role": "Polyketide synthase for DHN-melanin biosynthesis; linked to aggressiveness",
            "pmid_doi": "PMID:21356000; DOI:10.1111/j.1364-3703.2010.00666.x",
            "evidence_score": 0.72,
        },
        "CsFUS3": {
            "direct_exp": "Direct B. sorokiniana experimental evidence (targeted knockout delta-Csfus3 -> loss of pathogenicity on intact barley leaves, impaired conidiation & appressoria, oxidative stress hypersensitivity, Leng & Zhong, 2015)",
            "bio_role": "MAPK signalling component (PMK1/CHK1/FUS3 orthologue)",
            "pmid_doi": "PMID:26011429; DOI:10.1371/journal.pone.0128291",
            "evidence_score": 0.92,
        },
        "CsSTE11": {
            "direct_exp": "Direct B. sorokiniana experimental evidence (targeted knockout delta-Csste11 -> loss of pathogenicity on intact barley leaves, impaired conidiation & appressoria, oxidative stress hypersensitivity, Leng & Zhong, 2015)",
            "bio_role": "MAPKKK upstream signalling component (STE11 orthologue)",
            "pmid_doi": "PMID:26011429; DOI:10.1371/journal.pone.0128291",
            "evidence_score": 0.90,
        },
        "CsSp1": {
            "direct_exp": "Direct B. sorokiniana experimental evidence (split-marker knockout delta-Cssp1 in LK93 -> attenuated wheat pathogenicity, defective infectious hyphae growth, reduced conidiation, Zhang et al., 2022)",
            "bio_role": "Secreted effector essential for wheat colonisation and modulating host immunity",
            "pmid_doi": "PMID:34741560; DOI:10.1111/mpp.13155",
            "evidence_score": 0.84,
        },
    }

    sum_rows = []
    for tname, info in targets_cfg.items():
        pid = info["nd90pr_protein_id"]
        phi_sub = phi_hits_raw[phi_hits_raw["qseqid"] == pid] if not phi_hits_raw.empty else pd.DataFrame()
        best_phi = phi_sub.sort_values("bitscore", ascending=False).iloc[0] if not phi_sub.empty else None
        deg_sub = deg_hits[deg_hits["qseqid"] == pid] if not deg_hits.empty else pd.DataFrame()
        best_deg = deg_sub.sort_values("bitscore", ascending=False).iloc[0] if not deg_sub.empty else None
        ex_sub = expr_df[expr_df["target"] == tname] if not expr_df.empty else pd.DataFrame()
        ex_row = ex_sub.iloc[0] if not ex_sub.empty else None
        ipr_sub = ipr_df[ipr_df["target"] == tname] if not ipr_df.empty else pd.DataFrame()
        uni_sub = uni_df[uni_df["target"] == tname] if not uni_df.empty else pd.DataFrame()
        le = lit_evidence[tname]

        sum_rows.append(
            {
                "target": tname,
                "locus_tag": info["nd90pr_locus_tag"],
                "protein_id": pid,
                "uniprot_accession": uni_sub["uniprot_accession"].iloc[0] if not uni_sub.empty else "NA",
                "interpro_domains": "; ".join(f"{r['interpro_accession']}:{r['name']}" for _, r in ipr_sub.iterrows()) if not ipr_sub.empty else "None annotated in InterPro",
                "direct_species_evidence": le["direct_exp"],
                "direct_literature_ref": le["pmid_doi"],
                "phi_base_best_hit": f"{best_phi['sseqid']} (pident={best_phi['pident']}%, e={best_phi['evalue']})" if best_phi is not None else "No PHI-base hit (<1e-5)",
                "deg_essential_homolog_hit": f"{best_deg['sseqid']} (pident={best_deg['pident']}%, e={best_deg['evalue']})" if best_deg is not None else "No DEG hit (<1e-5)",
                "gse130517_axenic_rpkm": ex_row["bs_axenic_mean_rpkm"] if ex_row is not None else "NA",
                "gse130517_inplanta_mean_rpkm": ex_row["inplanta_mean_rpkm"] if ex_row is not None else "NA",
                "gse130517_inplanta_max_rpkm": ex_row["inplanta_max_rpkm"] if ex_row is not None else "NA",
                "functional_evidence_score_0_1": le["evidence_score"],
                "biological_role": le["bio_role"],
            }
        )

    sum_df = pd.DataFrame(sum_rows)
    sum_p = sum_dir / "functional_evidence_summary.tsv"
    sum_df.to_csv(sum_p, sep="\t", index=False)

    prov.record("M05", "Functional Evidence", str(target_prot_fa), sum_p, "DIAMOND+NCBI+UniProt+InterPro+GEO", "v2.2.8", "PHI-base_v4.19+GSE130517")
    logger.info("M05 complete: %d targets integrated with PHI-base, DEG, InterPro, UniProt, GSE130517, and %d RNAi machinery genes.", len(sum_df), len(rnai_df))
    return {
        "ncbi_gene_evidence": ncbi_gene_p,
        "ncbi_protein_evidence": ncbi_prot_p,
        "conserved_domains": cdd_p,
        "phi_base_evidence": phi_p,
        "interpro_annotations": ipr_p,
        "uniprot_annotations": uni_p,
        "infection_expression": expr_p,
        "functional_evidence_summary": sum_p,
    }
