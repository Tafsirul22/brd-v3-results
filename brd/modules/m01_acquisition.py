"""Module M01 — Data Acquisition for Bipolaris-RNAi Designer (BRD)."""

from __future__ import annotations

import csv
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import pandas as pd
import requests
from Bio import SeqIO

from brd.config import BRDConfig, load_yaml
from brd.integrations.ncbi import NCBIClient
from brd.logger import get_logger
from brd.provenance import ProvenanceTracker
from brd.utils.checksums import compute_md5, compute_sha256


def _verify_local_zip(zip_path: Path, extract_dir: Path, logger) -> List[Dict[str, str]]:
    """Extract local NCBI Datasets package without modifying original zip and verify MD5."""
    extract_dir.mkdir(parents=True, exist_ok=True)
    marker = extract_dir / ".extracted_ok"
    if not marker.exists():
        logger.info("Extracting local genome package %s -> %s", zip_path, extract_dir)
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(extract_dir)
        marker.write_text("ok\n", encoding="utf-8")
    else:
        logger.info("Local genome package already extracted at %s", extract_dir)

    md5_file = extract_dir / "md5sum.txt"
    records = []
    if md5_file.exists():
        for line in md5_file.read_text(encoding="utf-8").splitlines():
            parts = line.strip().split()
            if len(parts) >= 2:
                expected_md5, rel_path = parts[0], parts[1]
                fpath = extract_dir / rel_path
                actual_md5 = compute_md5(fpath) if fpath.exists() else "MISSING"
                records.append(
                    {
                        "file": rel_path,
                        "expected_md5": expected_md5,
                        "actual_md5": actual_md5,
                        "md5_match": str(expected_md5 == actual_md5),
                    }
                )
    return records


def run_m01(cfg: BRDConfig) -> Dict[str, Path]:
    log_path = cfg.output_dir / "logs" / "brd_execution.log"
    cmd_log = cfg.output_dir / "logs" / "commands.jsonl"
    logger = get_logger("BRD.M01", log_path)
    prov = ProvenanceTracker(cfg.output_dir)

    m01_dir = cfg.module_dir("01_data_acquisition")
    ext_dir = m01_dir / "external_references"
    ext_dir.mkdir(parents=True, exist_ok=True)
    local_dir = m01_dir / "local_nd90pr"
    genomes_dir = m01_dir / "genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)

    # Sync configuration snapshot to 00_project_metadata
    meta_cfg = cfg.output_dir / "00_project_metadata" / "project_configuration.yaml"
    cfg_src = cfg.output_dir / "config" / "project_configuration.yaml"
    if cfg_src.exists():
        meta_cfg.write_text(cfg_src.read_text(encoding="utf-8"), encoding="utf-8")

    # 1. Verify and extract local zip
    local_zip = Path(cfg.raw["project"]["local_zip"])
    md5_checks = _verify_local_zip(local_zip, local_dir, logger)
    md5_ok = all(r["md5_match"] == "True" for r in md5_checks)
    logger.info("Local zip MD5 verification: %s (%d files checked)", md5_ok, len(md5_checks))
    pd.DataFrame(md5_checks).to_csv(m01_dir / "local_zip_md5_verification.tsv", sep="\t", index=False)

    # 2. Query NCBI Datasets for all Bipolaris sorokiniana assemblies (taxid 45130)
    ncbi = NCBIClient(
        email=cfg.ncbi_email,
        api_key=cfg.ncbi_api_key,
        datasets_bin=cfg.tool_path("datasets"),
        log_path=cmd_log,
    )
    reports_json = m01_dir / "ncbi_taxon45130_reports.json"
    if not reports_json.exists():
        reports = ncbi.fetch_taxon_assemblies(cfg.raw["project"]["taxonomy_id"])
        reports_json.write_text(json.dumps(reports, indent=2), encoding="utf-8")
    else:
        reports = json.loads(reports_json.read_text(encoding="utf-8"))

    now_str = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    genome_rows = []
    isolate_rows = []
    accession_rows = []
    download_rows = []

    seen_biological = set()
    for rep in reports:
        acc = rep.get("accession", "")
        paired = rep.get("paired_accession", "")
        org = rep.get("organism", {}).get("organism_name", "")
        inf = rep.get("organism", {}).get("infraspecific_names", {})
        strain = inf.get("strain", "") or inf.get("isolate", "") or acc
        info = rep.get("assembly_info", {})
        stats = rep.get("assembly_stats", {})
        annot = rep.get("annotation_info", {})
        bs = info.get("biosample", {})
        attrs = {a["name"]: a["value"] for a in bs.get("attributes", [])}
        sra_ids = [s.get("value", "") for s in bs.get("sample_ids", []) if s.get("db") == "SRA"]
        geo = attrs.get("geo_loc_name", "Unknown")
        host = attrs.get("host", "Unknown")
        level = info.get("assembly_level", "")
        bp = info.get("bioproject_accession", "")
        bs_acc = bs.get("accession", "")
        tech = info.get("sequencing_tech", "")
        method = info.get("assembly_method", "")
        gsize = int(stats.get("total_sequence_length", 0) or 0)
        annot_status = annot.get("status", "None")
        annot_provider = annot.get("provider", "None")
        gene_count = annot.get("stats", {}).get("gene_counts", {}).get("total", 0)

        # Deduplicate paired GenBank/RefSeq GCA_000338995.1 vs GCF_000338995.1
        is_primary_pair = acc == "GCF_000338995.1"
        is_redundant_pair = acc == "GCA_000338995.1"
        bio_key = f"{strain}".strip().lower()

        include_flag = True
        excl_reason = ""
        if is_redundant_pair:
            include_flag = False
            excl_reason = "Paired GenBank representation of GCF_000338995.1 (ND90Pr); excluded to prevent duplicate isolate counting."
        elif "bipolaris sorokiniana" not in org.lower() and "cochliobolus sativus" not in org.lower():
            include_flag = False
            excl_reason = f"Non-target taxonomy: {org}"

        # Resolve local or downloaded path
        genome_fasta = ""
        gff_path = ""
        cds_path = ""
        prot_path = ""
        rna_path = ""

        if acc == "GCF_000338995.1":
            pkg_dir = local_dir / "ncbi_dataset" / "data" / "GCF_000338995.1"
            source_label = "Local NCBI Datasets ZIP (RefSeq)"
        elif include_flag:
            pkg_dir = genomes_dir / acc
            inc = ["genome"]
            if (annot_status and annot_status != "None") or int(gene_count or 0) > 0:
                inc.extend(["gff3", "cds", "protein", "rna"])
            try:
                ncbi.download_genome_package(acc, inc, pkg_dir, overwrite=(acc == "GCA_013416765.1"))
                source_label = "NCBI Datasets CLI download"
            except Exception as exc:
                logger.warning("Failed downloading %s: %s", acc, exc)
                include_flag = False
                excl_reason = f"Download failed: {exc}"
                source_label = "Download failed"
            pkg_dir = pkg_dir / "ncbi_dataset" / "data" / acc
        else:
            pkg_dir = Path("")
            source_label = "Excluded paired assembly"

        if include_flag and pkg_dir.exists():
            fnas = [p for p in sorted(pkg_dir.glob("*_genomic.fna")) if not p.name.startswith("cds_")]
            if fnas:
                genome_fasta = str(fnas[0])
            gffs = sorted(pkg_dir.glob("*.gff*"))
            if gffs:
                gff_path = str(gffs[0])
            cdss = sorted(pkg_dir.glob("cds_from_genomic.fna"))
            if cdss:
                cds_path = str(cdss[0])
            prots = sorted(pkg_dir.glob("protein.faa"))
            if prots:
                prot_path = str(prots[0])
            rnas = sorted(pkg_dir.glob("rna.fna"))
            if rnas:
                rna_path = str(rnas[0])

            if genome_fasta:
                sha = compute_sha256(Path(genome_fasta))
                download_rows.append(
                    {
                        "dataset_id": acc,
                        "category": "Pathogen genome",
                        "organism": org,
                        "strain_isolate": strain,
                        "source": source_label,
                        "file_path": genome_fasta,
                        "size_bytes": Path(genome_fasta).stat().st_size,
                        "sha256": sha,
                        "download_date": now_str,
                    }
                )

        genome_rows.append(
            {
                "assembly_accession": acc,
                "paired_accession": paired,
                "organism": org,
                "isolate_id": strain.replace(" ", "_").replace("(", "").replace(")", ""),
                "strain": strain,
                "biosample": bs_acc,
                "bioproject": bp,
                "sra_accession": ";".join(sra_ids) if sra_ids else "None_linked",
                "geographic_origin": geo,
                "host": host,
                "sequencing_platform": tech,
                "assembly_method": method,
                "assembly_level": level,
                "genome_size_bp": gsize,
                "contig_n50": stats.get("contig_n50", ""),
                "annotation_source": annot_provider,
                "annotation_status": annot_status,
                "gene_count": gene_count,
                "release_date": info.get("release_date", ""),
                "download_date": now_str,
                "include_in_panel": include_flag,
                "exclusion_reason": excl_reason,
                "genome_fasta": genome_fasta,
                "gff_path": gff_path,
                "cds_path": cds_path,
                "protein_path": prot_path,
                "rna_path": rna_path,
            }
        )
        if include_flag and bio_key not in seen_biological:
            seen_biological.add(bio_key)
            isolate_rows.append(
                {
                    "isolate_id": strain.replace(" ", "_").replace("(", "").replace(")", ""),
                    "strain": strain,
                    "assembly_accession": acc,
                    "biosample": bs_acc,
                    "bioproject": bp,
                    "geographic_origin": geo,
                    "host": host,
                    "assembly_level": level,
                    "has_ncbi_annotation": bool(gff_path),
                    "genome_fasta": genome_fasta,
                    "gff_path": gff_path,
                    "protein_path": prot_path,
                }
            )

    # 3. Retrieve supplied GenBank accessions (CsPKS1 + target mRNAs/proteins)
    pks1_accs = ["HQ830033.1"] + [f"KF7691{i}.1" for i in range(34, 42)]
    pks1_gb_path = ext_dir / "CsPKS1_supplied_records.gb"
    pks1_fa_path = ext_dir / "CsPKS1_supplied_records.fasta"
    if not pks1_gb_path.exists():
        gb_text = ncbi.efetch_text("nuccore", pks1_accs, rettype="gb", retmode="text")
        pks1_gb_path.write_text(gb_text, encoding="utf-8")
    gb_recs = list(SeqIO.parse(str(pks1_gb_path), "genbank"))
    SeqIO.write(gb_recs, str(pks1_fa_path), "fasta")

    for r in gb_recs:
        strain = ""
        isolate = ""
        cds_len = 0
        prot_id = ""
        for f in r.features:
            if f.type == "source":
                strain = ",".join(f.qualifiers.get("strain", []))
                isolate = ",".join(f.qualifiers.get("isolate", []))
            if f.type == "CDS":
                cds_len = len(f.extract(r.seq))
                prot_id = ",".join(f.qualifiers.get("protein_id", []))
        accession_rows.append(
            {
                "target": "CsPKS1",
                "accession": r.id,
                "database": "NCBI nuccore",
                "organism": r.annotations.get("organism", ""),
                "strain_or_isolate": strain or isolate,
                "length_bp": len(r.seq),
                "cds_length_bp": cds_len,
                "protein_id": prot_id,
                "description": r.description,
                "status": "Verified",
            }
        )

    # Record assembly accessions in accession_manifest
    for gr in genome_rows:
        accession_rows.append(
            {
                "target": "Genome_Panel",
                "accession": gr["assembly_accession"],
                "database": "NCBI Assembly",
                "organism": gr["organism"],
                "strain_or_isolate": gr["strain"],
                "length_bp": gr["genome_size_bp"],
                "cds_length_bp": gr["gene_count"],
                "protein_id": "",
                "description": f"{gr['assembly_level']} assembly ({gr['geographic_origin']})",
                "status": "Included" if gr["include_in_panel"] else f"Excluded ({gr['exclusion_reason']})",
            }
        )

    # 4. Download non-target reference transcriptomes (Wheat, Barley, Beneficial fungi)
    refs_cfg = load_yaml(cfg.output_dir / "config" / "references.yaml")
    nontarget_dir = cfg.output_dir / "08_off_target_screening" / "blast_databases" / "sources"
    nontarget_dir.mkdir(parents=True, exist_ok=True)
    for grp_name, entries in refs_cfg.get("non_target_panel", {}).items():
        for entry in entries:
            acc = entry["accession"]
            inc = entry["include"]
            dest = nontarget_dir / entry["id"]
            try:
                ncbi.download_genome_package(acc, inc, dest)
                data_sub = dest / "ncbi_dataset" / "data" / acc
                # Prefer rna.fna, fallback to cds_from_genomic.fna
                rna_f = data_sub / "rna.fna"
                cds_f = data_sub / "cds_from_genomic.fna"
                src_f = rna_f if rna_f.exists() else cds_f
                if src_f.exists():
                    download_rows.append(
                        {
                            "dataset_id": entry["id"],
                            "category": f"Non-target ({grp_name})",
                            "organism": entry["organism"],
                            "strain_isolate": acc,
                            "source": "NCBI Datasets CLI",
                            "file_path": str(src_f),
                            "size_bytes": src_f.stat().st_size,
                            "sha256": compute_sha256(src_f),
                            "download_date": now_str,
                        }
                    )
            except Exception as exc:
                logger.error("Failed downloading non-target %s (%s): %s", entry["id"], acc, exc)

    # 5. Download PHI-base current FASTA + CSV and GSE130517 expression table
    phi_dir = cfg.output_dir / "05_functional_evidence" / "phi_base"
    phi_dir.mkdir(parents=True, exist_ok=True)
    for fname in ["phi-base_current.fas", "phi-base_current.csv"]:
        fp = phi_dir / fname
        if not fp.exists():
            url = f"https://raw.githubusercontent.com/PHI-base/data/master/releases/{fname}"
            logger.info("Downloading %s", url)
            resp = requests.get(url, timeout=180)
            resp.raise_for_status()
            fp.write_bytes(resp.content)
        download_rows.append(
            {
                "dataset_id": fname,
                "category": "Functional DB (PHI-base)",
                "organism": "Multi-species pathogens",
                "strain_isolate": "v4.19",
                "source": "GitHub PHI-base/data",
                "file_path": str(fp),
                "size_bytes": fp.stat().st_size,
                "sha256": compute_sha256(fp),
                "download_date": now_str,
            }
        )

    expr_dir = cfg.output_dir / "05_functional_evidence" / "expression"
    expr_dir.mkdir(parents=True, exist_ok=True)
    gse_fname = "GSE130517_Bipolaris_sorokiniana_blast2GOannotation_rawrc_rpkm_rei.txt.gz"
    gse_fp = expr_dir / gse_fname
    if not gse_fp.exists():
        gse_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130517/suppl/{gse_fname}"
        logger.info("Downloading %s", gse_url)
        resp = requests.get(gse_url, timeout=120)
        resp.raise_for_status()
        gse_fp.write_bytes(resp.content)
    download_rows.append(
        {
            "dataset_id": "GSE130517",
            "category": "Expression matrix (GEO)",
            "organism": "Bipolaris sorokiniana ND90Pr",
            "strain_isolate": "ND90Pr on Hordeum vulgare",
            "source": "NCBI GEO FTP",
            "file_path": str(gse_fp),
            "size_bytes": gse_fp.stat().st_size,
            "sha256": compute_sha256(gse_fp),
            "download_date": now_str,
        }
    )

    # Write M01 TSV outputs
    out_genome = m01_dir / "genome_metadata.tsv"
    out_isolate = m01_dir / "isolate_metadata.tsv"
    out_acc = m01_dir / "accession_manifest.tsv"
    out_dl = m01_dir / "download_manifest.tsv"
    pd.DataFrame(genome_rows).to_csv(out_genome, sep="\t", index=False)
    pd.DataFrame(isolate_rows).to_csv(out_isolate, sep="\t", index=False)
    pd.DataFrame(accession_rows).to_csv(out_acc, sep="\t", index=False)
    pd.DataFrame(download_rows).to_csv(out_dl, sep="\t", index=False)

    prov.record("M01", "Data Acquisition", str(local_zip), out_genome, "datasets+Entrez", "18.41.0", "taxid:45130", f"{len(isolate_rows)} biological isolates")
    prov.record("M01", "Download Manifest", "NCBI+PHI-base+GEO", out_dl, "datasets+requests", "18.41.0", "multi-source", f"{len(download_rows)} files checksummed")
    logger.info("M01 complete: %d biological isolates, %d downloaded datasets.", len(isolate_rows), len(download_rows))

    return {
        "genome_metadata": out_genome,
        "isolate_metadata": out_isolate,
        "accession_manifest": out_acc,
        "download_manifest": out_dl,
    }
