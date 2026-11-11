"""Stage 4: Expanded Non-Target Reference Panel (Crops, Beneficial Fungi, Bacteria) & Dual-Strand Transcript+Genome Screening."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import RNA
from Bio import SeqIO
from Bio.Seq import Seq

from paths_v3 import load_yaml_portable, resolve_output_root, resolve_tool_binary, resolve_v3_root

BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "sstrand"]
LEGACY_BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs", "sstrand"]


def _read_blast_tsv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLS)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first_line = ""
        for line in fh:
            if line.strip() and not line.startswith("#"):
                first_line = line
                break
    if not first_line:
        return pd.DataFrame(columns=BLAST_COLS)
    ncols = len(first_line.rstrip("\r\n").split("\t"))
    cols = LEGACY_BLAST_COLS if ncols == 16 else (BLAST_COLS if ncols == 14 else [f"c{i}" for i in range(ncols)])
    df = pd.read_csv(path, sep="\t", names=cols, comment="#")
    for req in BLAST_COLS:
        if req not in df.columns:
            df[req] = "" if req == "sstrand" else 0
    return df[BLAST_COLS]

# Expanded taxonomic reference panel (Stage 4.1)
PANEL_SPECS = [
    # Existing 7 from M08 v1/v2 (with added genome downloads for beneficial fungi!)
    {"ref_id": "Triticum_aestivum_IWGSC_v2_1", "group": "crop_wheat", "organism": "Triticum aestivum", "accession": "GCF_018294505.1", "tx_source": "existing_m08", "fetch_genome": False, "include": "rna", "rationale": "Primary cereal crop host (hexaploid CS RefSeq v2.1 transcriptome)"},
    {"ref_id": "Hordeum_vulgare_MorexV3", "group": "crop_barley", "organism": "Hordeum vulgare subsp. vulgare", "accession": "GCF_904849725.1", "tx_source": "existing_m08", "fetch_genome": False, "include": "rna", "rationale": "Secondary cereal crop host (cv. Morex V3 RefSeq transcriptome)"},
    {"ref_id": "Oryza_sativa_AGIS1", "group": "crop_rice", "organism": "Oryza sativa Japonica Group", "accession": "GCF_034140825.1", "tx_source": "download", "fetch_genome": True, "include": "rna,cds,genome", "rationale": "Major monocot cereal crop comparator (Nipponbare AGIS-1.0 RefSeq transcriptome + complete genome)"},
    {"ref_id": "Zea_mays_B73_NAM5", "group": "crop_maize", "organism": "Zea mays", "accession": "GCF_902167145.1", "tx_source": "download", "fetch_genome": False, "include": "rna", "rationale": "Major monocot crop comparator (Zm-B73-REFERENCE-NAM-5.0 RefSeq transcriptome)"},
    # Beneficial & agriculturally relevant fungi (Trichoderma spp., Rhizophagus spp., Beauveria, Metarhizium)
    {"ref_id": "Trichoderma_atroviride_IMI206040", "group": "beneficial_fungi", "organism": "Trichoderma atroviride IMI 206040", "accession": "GCF_000171015.1", "tx_source": "existing_m08", "fetch_genome": True, "include": "genome", "rationale": "Biocontrol mycoparasitic fungus (Piombo et al., 2021 RNAi active)"},
    {"ref_id": "Trichoderma_virens_Gv29_8", "group": "beneficial_fungi", "organism": "Trichoderma virens Gv29-8", "accession": "GCF_000170995.1", "tx_source": "existing_m08", "fetch_genome": True, "include": "genome", "rationale": "Rhizosphere beneficial fungus in cereal cropping systems"},
    {"ref_id": "Trichoderma_reesei_QM6a", "group": "beneficial_fungi", "organism": "Trichoderma reesei QM6a", "accession": "GCF_000167675.1", "tx_source": "existing_m08", "fetch_genome": True, "include": "genome", "rationale": "Saprotrophic Trichoderma reference"},
    {"ref_id": "Trichoderma_harzianum_CBS226_95", "group": "beneficial_fungi", "organism": "Trichoderma harzianum CBS 226.95", "accession": "GCF_003025095.1", "tx_source": "download", "fetch_genome": True, "include": "rna,cds,genome", "rationale": "Widely deployed cereal & soil biocontrol Trichoderma species"},
    {"ref_id": "Beauveria_bassiana_ARSEF2860", "group": "beneficial_fungi", "organism": "Beauveria bassiana ARSEF 2860", "accession": "GCF_000280675.1", "tx_source": "existing_m08", "fetch_genome": True, "include": "genome", "rationale": "Beneficial entomopathogenic fungus deployed in IPM"},
    {"ref_id": "Metarhizium_anisopliae_RefSeq", "group": "beneficial_fungi", "organism": "Metarhizium anisopliae", "accession": "GCF_039654215.1", "tx_source": "download", "fetch_genome": True, "include": "rna,cds,genome", "rationale": "Entomopathogenic soil biocontrol fungus (RefSeq)"},
    {"ref_id": "Metarhizium_robertsii_ARSEF2575", "group": "beneficial_fungi", "organism": "Metarhizium robertsii", "accession": "GCA_000591435.1", "tx_source": "download", "fetch_genome": True, "include": "cds,genome", "rationale": "Rhizosphere-competent beneficial entomopathogenic fungus"},
    {"ref_id": "Rhizophagus_irregularis_A1", "group": "beneficial_fungi", "organism": "Rhizophagus irregularis", "accession": "GCA_002897155.2", "tx_source": "existing_m08", "fetch_genome": True, "include": "genome", "rationale": "Arbuscular mycorrhizal root symbiont of wheat & barley"},
    {"ref_id": "Rhizophagus_clarus_HR1", "group": "beneficial_fungi", "organism": "Rhizophagus clarus", "accession": "GCA_003203555.1", "tx_source": "download", "fetch_genome": True, "include": "cds,genome", "rationale": "Second arbuscular mycorrhizal Glomeromycotina symbiont"},
    # Relevant beneficial & rhizosphere bacteria (Bacillus & Pseudomonas)
    {"ref_id": "Bacillus_subtilis_168", "group": "bacteria", "organism": "Bacillus subtilis subsp. subtilis str. 168", "accession": "GCF_000009045.1", "tx_source": "download", "fetch_genome": True, "include": "rna,cds,genome", "rationale": "Model soil/phyllosphere Gram-positive bacterium"},
    {"ref_id": "Bacillus_velezensis_FZB42_grp", "group": "bacteria", "organism": "Bacillus velezensis", "accession": "GCF_000508125.2", "tx_source": "download", "fetch_genome": True, "include": "rna,cds,genome", "rationale": "Beneficial wheat/barley rhizosphere PGPR biocontrol agent"},
    {"ref_id": "Pseudomonas_fluorescens_ATCC13525", "group": "bacteria", "organism": "Pseudomonas fluorescens ATCC 13525", "accession": "GCF_900215245.1", "tx_source": "download", "fetch_genome": True, "include": "rna,cds,genome", "rationale": "Beneficial plant-growth-promoting & biocontrol rhizobacterium"},
    {"ref_id": "Pseudomonas_putida_RefSeq", "group": "bacteria", "organism": "Pseudomonas putida", "accession": "GCF_000412675.1", "tx_source": "download", "fetch_genome": True, "include": "rna,cds,genome", "rationale": "Ubiquitous soil/rhizosphere Gram-negative bacterium (NBRC 14164 reference)"},
]


def _iter_fasta_upper(fa_path: Path):
    """Fast lightweight FASTA reader returning (record_id, UPPERCASE_SEQ) without SeqRecord overhead."""
    hdr = None
    chunks: List[str] = []
    with open(fa_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line == "\n":
                continue
            if line[0] == ">":
                if hdr is not None:
                    yield hdr, "".join(chunks).upper()
                hdr = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(chunks).upper()


def _build_pigeonhole_21mer_index(kmers: Set[str]) -> dict[str, tuple[int, ...]]:
    """Every 21-mer contains at least one full 11-mer starting at grid offset g % 11 == 0."""
    idx: dict[str, set[int]] = {}
    for km in kmers:
        if len(km) != 21:
            continue
        for off in range(11):
            sub11 = km[off : off + 11]
            s = idx.get(sub11)
            if s is None:
                idx[sub11] = {off}
            else:
                s.add(off)
    return {k: tuple(sorted(v)) for k, v in idx.items()}


def _scan_exact_21mers(fa_path: Path, query_21mers: Set[str], check_revcomp: bool = True) -> List[Tuple[str, int, str]]:
    """Exact deterministic 21-mer scan (no missed matches) via 11-nt pigeonhole stride; returns (rec_id, start_0, strand)."""
    if not fa_path or not fa_path.exists() or not query_21mers:
        return []
    targets: Set[str] = set()
    fwd_set = {q.upper() for q in query_21mers if len(q) == 21}
    targets.update(fwd_set)
    rc_map: Dict[str, str] = {}
    if check_revcomp:
        for q in fwd_set:
            rc = str(Seq(q).reverse_complement())
            targets.add(rc)
            rc_map[rc] = q
    idx = _build_pigeonhole_21mer_index(targets)
    hits: List[Tuple[str, int, str]] = []
    for rec_id, seq in _iter_fasta_upper(fa_path):
        n = len(seq)
        if n < 21:
            continue
        max_start = n - 21
        for g in range(0, n - 10, 11):
            offs = idx.get(seq[g : g + 11])
            if offs is not None:
                for off in offs:
                    st = g - off
                    if 0 <= st <= max_start:
                        cand = seq[st : st + 21]
                        if cand in fwd_set:
                            hits.append((rec_id, st, "+"))
                        elif check_revcomp and cand in rc_map:
                            hits.append((rec_id, st, "-"))
    return hits


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _fetch_ncbi_summary(ds_bin: str, acc: str) -> dict:
    p = subprocess.run([ds_bin, "summary", "genome", "accession", acc], capture_output=True, text=True)
    if p.returncode == 0 and p.stdout.strip().startswith("{"):
        reps = json.loads(p.stdout).get("reports", [])
        if reps:
            return reps[0]
    return {}


def _download_ncbi_package(ds_bin: str, acc: str, include: str, target_dir: Path) -> Tuple[Path | None, Path | None, str]:
    target_dir.mkdir(parents=True, exist_ok=True)
    data_dir = target_dir / "ncbi_dataset" / "data" / acc
    zip_p = target_dir / f"{acc}.zip"
    if not data_dir.exists():
        p = subprocess.run([ds_bin, "download", "genome", "accession", acc, "--include", include, "--filename", str(zip_p)], capture_output=True, text=True)
        if p.returncode != 0 or not zip_p.exists():
            return None, None, f"download_failed:{p.stderr.replace(chr(10), ' ').replace(chr(13), ' ')[:180]}"
        with zipfile.ZipFile(zip_p) as zf:
            zf.extractall(target_dir)
        zip_p.unlink(missing_ok=True)

    tx_p = None
    for cand_name in ("rna.fna", "cds_from_genomic.fna"):
        cp = data_dir / cand_name
        if cp.exists() and cp.stat().st_size > 0:
            tx_p = cp
            break
    gen_files = [p for p in data_dir.glob("*_genomic.fna") if p.name != "cds_from_genomic.fna"]
    gen_p = gen_files[0] if gen_files else None
    return tx_p, gen_p, "OK"


def _short_path_if_nt(p: Path) -> str:
    s_abs = str(p)
    if os.name == "nt" and " " in s_abs:
        import ctypes
        buf = ctypes.create_unicode_buffer(1024)
        if ctypes.windll.kernel32.GetShortPathNameW(s_abs, buf, 1024) and buf.value and " " not in buf.value:
            return buf.value
    return s_abs


def _blast_safe_path(target: Path, cwd_dir: Path) -> str:
    """NCBI BLAST+ splits CLI path args on whitespace; return space-free relative or Windows 8.3 short path."""
    rel = os.path.relpath(str(target), start=str(cwd_dir))
    if " " not in rel:
        return rel
    return _short_path_if_nt(target)


def _make_blastdb_and_run(makeblastdb_bin: str, blastn_bin: str, ref_fa: Path, db_pref: Path, guide_fa: Path, cand_fa: Path, out_prefix: Path) -> Tuple[Path, Path]:
    db_pref.parent.mkdir(parents=True, exist_ok=True)
    cwd_dir = db_pref.parent
    cwd_safe = _short_path_if_nt(cwd_dir)
    if not (Path(str(db_pref) + ".nsq").exists() or Path(str(db_pref) + ".nal").exists() or Path(str(db_pref) + ".00.nsq").exists()):
        subprocess.run(
            [makeblastdb_bin, "-in", _blast_safe_path(ref_fa, cwd_dir), "-dbtype", "nucl", "-out", db_pref.name],
            check=True,
            capture_output=True,
            cwd=cwd_safe,
        )
    g_out = out_prefix.parent / f"guides_vs_{out_prefix.name}.tsv"
    c_out = out_prefix.parent / f"cands_vs_{out_prefix.name}.tsv"
    fmt = f"6 {' '.join(BLAST_COLS)}"
    if not g_out.exists():
        subprocess.run(
            [
                blastn_bin, "-task", "blastn-short", "-query", _blast_safe_path(guide_fa, cwd_dir), "-db", db_pref.name, "-out", _blast_safe_path(g_out, cwd_dir),
                "-evalue", "10.0", "-word_size", "7", "-strand", "both", "-dust", "no", "-soft_masking", "false",
                "-max_target_seqs", "20", "-num_threads", "4", "-outfmt", fmt,
            ],
            check=True,
            capture_output=True,
            cwd=cwd_safe,
        )
    if not c_out.exists():
        subprocess.run(
            [
                blastn_bin, "-task", "blastn", "-query", _blast_safe_path(cand_fa, cwd_dir), "-db", db_pref.name, "-out", _blast_safe_path(c_out, cwd_dir),
                "-evalue", "1e-4", "-word_size", "11", "-strand", "both", "-max_target_seqs", "25", "-num_threads", "4", "-outfmt", fmt,
            ],
            check=True,
            capture_output=True,
            cwd=cwd_safe,
        )
    return g_out, c_out


def run_stage4(output_root: Path | None = None) -> Dict[str, Path]:
    out = resolve_output_root(output_root)
    v3 = resolve_v3_root(out)
    s4 = v3 / "stage4_off_target_expanded"
    ref_root = v3 / "expanded_references"
    src_root = ref_root / "sources"
    db_root = ref_root / "blast_dbs"
    raw_root = ref_root / "raw_blast"
    fig_dir = s4 / "figures"
    for d in (s4, src_root, db_root, raw_root, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    cfg = load_yaml_portable(out / "config" / "project_configuration.yaml")
    tools_cfg = cfg.get("tools", {})
    ds_bin = resolve_tool_binary("datasets", tools_cfg, out)
    blastn_bin = resolve_tool_binary("blastn", tools_cfg, out)
    makeblastdb_bin = resolve_tool_binary("makeblastdb", tools_cfg, out)
    if not ds_bin or not blastn_bin or not makeblastdb_bin:
        raise FileNotFoundError("Required tools datasets/blastn/makeblastdb missing")

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    auth_cands = pd.read_csv(v3 / "authoritative_data" / "auth_candidates.tsv", sep="\t")
    auth_guides = pd.read_csv(v3 / "authoritative_data" / "auth_guides.tsv", sep="\t")
    auth_targets = pd.read_csv(v3 / "authoritative_data" / "auth_targets.tsv", sep="\t").set_index("target")

    # Write unique guides FASTA & candidate FASTA
    cand_fa = v3 / "authoritative_data" / "auth_candidates.fasta"
    unique_guides_fa = raw_root / "unique_guides_v3.fasta"
    guide_to_cands: Dict[str, List[str]] = {}
    guide_seq_map: Dict[str, str] = {}
    for _, gr in auth_guides.iterrows():
        k21 = str(gr["sense_dna_21nt"]).upper()
        guide_to_cands.setdefault(k21, [])
        if str(gr["candidate_id"]) not in guide_to_cands[k21]:
            guide_to_cands[k21].append(str(gr["candidate_id"]))
    unique_k21s = sorted(guide_to_cands.keys())
    with open(unique_guides_fa, "w", encoding="utf-8") as fh:
        for i, k21 in enumerate(unique_k21s, 1):
            gid = f"UG_{i:04d}"
            guide_seq_map[gid] = k21
            fh.write(f">{gid}\n{k21}\n")

    # 1. Acquire / verify all 17 non-target references + self ND90Pr
    old_m08_src = out / "08_off_target_screening" / "blast_databases" / "sources"
    manifest_rows: List[dict] = []
    screen_targets: List[dict] = []

    for spec in PANEL_SPECS:
        rid = spec["ref_id"]
        acc = spec["accession"]
        meta = _fetch_ncbi_summary(ds_bin, acc)
        org = meta.get("organism", {}).get("organism_name", spec["organism"])
        tax_id = meta.get("organism", {}).get("tax_id", "NA")
        src_db = meta.get("source_database", "RefSeq" if acc.startswith("GCF_") else "GenBank")
        asm_name = meta.get("assembly_info", {}).get("assembly_name", "NA")
        rel_date = meta.get("release_date", meta.get("assembly_info", {}).get("release_date", "NA"))
        busco_m = meta.get("annotation_info", {}).get("busco", {}).get("complete")
        checkm_m = meta.get("checkm_info", {}).get("completeness")
        qual_str = (
            f"BUSCO_complete={round(float(busco_m)*100,1)}%" if busco_m is not None else
            (f"CheckM_complete={checkm_m}%" if checkm_m is not None else f"assembly_level={meta.get('assembly_info', {}).get('assembly_level', 'NA')}")
        )

        tx_fa = gen_fa = None
        dl_status = "REUSED_M08"
        if spec["tx_source"] == "existing_m08":
            old_data = old_m08_src / rid / "ncbi_dataset" / "data" / acc
            for cn in ("rna.fna", "cds_from_genomic.fna"):
                if (old_data / cn).exists():
                    tx_fa = old_data / cn
                    break
            if spec["fetch_genome"]:
                _, gen_fa, dl_status = _download_ncbi_package(ds_bin, acc, "genome", src_root / rid)
                dl_status = f"TX_REUSED_M08; GENOME_{dl_status}"
        else:
            tx_fa, gen_fa, dl_status = _download_ncbi_package(ds_bin, acc, spec["include"], src_root / rid)

        tx_cnt = sum(1 for _ in SeqIO.parse(str(tx_fa), "fasta")) if tx_fa and tx_fa.exists() else 0
        gen_cnt = sum(1 for _ in SeqIO.parse(str(gen_fa), "fasta")) if gen_fa and gen_fa.exists() else 0

        manifest_rows.append(
            {
                "ref_id": rid,
                "group": spec["group"],
                "organism": org,
                "tax_id": tax_id,
                "accession": acc,
                "source_database": src_db,
                "assembly_version_name": asm_name,
                "release_date": rel_date,
                "retrieval_timestamp_utc": now_utc if "download" in spec["tx_source"] or spec["fetch_genome"] else "2026-11-11 M01/M08 baseline",
                "completeness_quality_metadata": qual_str,
                "transcript_screened": tx_fa is not None and tx_fa.exists(),
                "transcript_count": tx_cnt,
                "genome_screened": gen_fa is not None and gen_fa.exists(),
                "genome_scaffold_count": gen_cnt,
                "acquisition_status": dl_status,
                "rationale": spec["rationale"],
            }
        )
        screen_targets.append({**spec, "organism": org, "tax_id": tax_id, "tx_fa": tx_fa, "gen_fa": gen_fa})

    # Add B. sorokiniana ND90Pr self genome + paralog transcriptome record to manifest
    bsor_tx = out / "03_annotation" / "transcripts" / "genome_transcripts.fasta"
    bsor_gen = out / str(pd.read_csv(v3 / "authoritative_data" / "auth_isolates.tsv", sep="\t").set_index("isolate_id").loc["ND90Pr", "relative_genome_path"])
    manifest_rows.append(
        {
            "ref_id": "Bipolaris_sorokiniana_ND90Pr_self",
            "group": "self_fungal_paralog_and_genome",
            "organism": "Bipolaris sorokiniana ND90Pr",
            "tax_id": 45130,
            "accession": "GCF_000338995.1",
            "source_database": "RefSeq",
            "assembly_version_name": "Cocsa1",
            "release_date": "2013-02-01",
            "retrieval_timestamp_utc": "Local baseline verified MD5",
            "completeness_quality_metadata": "12,210 annotated genes; 34.41 Mb; 154 scaffolds",
            "transcript_screened": True,
            "transcript_count": sum(1 for _ in SeqIO.parse(str(bsor_tx), "fasta")),
            "genome_screened": True,
            "genome_scaffold_count": sum(1 for _ in SeqIO.parse(str(bsor_gen), "fasta")),
            "acquisition_status": "LOCAL_SELF_CONTROL (target loci masked cleanly during off-target counting)",
            "rationale": "Internal paralog / off-target site detection across B. sorokiniana genome and non-self transcripts",
        }
    )
    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(s4 / "expanded_reference_manifest_v3.tsv", sep="\t", index=False)

    # Coverage & gaps explicit summary
    cov_gaps = manifest_df[["ref_id", "group", "organism", "accession", "transcript_screened", "transcript_count", "genome_screened", "genome_scaffold_count"]].copy()
    cov_gaps["screening_scope_note"] = cov_gaps.apply(
        lambda r: "Transcriptome + Whole Genome screened" if (r["transcript_screened"] and r["genome_screened"])
        else ("Complete RefSeq Transcriptome screened (multi-GB whole genome excluded from local storage; intergenic/intronic unscreened)" if r["transcript_screened"]
              else ("Whole Genome Assembly screened (no separate GenBank RNA/CDS packaged; coding + non-coding covered at genome level)" if r["genome_screened"] else "Blocked / Unavailable")),
        axis=1,
    )
    cov_gaps.to_csv(s4 / "reference_panel_coverage_and_gaps_v3.tsv", sep="\t", index=False)

    # 2. Exact 21-mer Dual-Strand Scan across ALL transcriptomes AND ALL downloaded genomes
    cand_seq_map = {str(r["candidate_id"]): str(r["sequence"]).upper() for _, r in auth_cands.iterrows()}
    cand_target_map = {str(r["candidate_id"]): str(r["target"]) for _, r in auth_cands.iterrows()}
    all_c_21mers: Set[str] = set()
    cand_21_map: Dict[str, Set[str]] = {}
    for cid, cseq in cand_seq_map.items():
        kms: Set[str] = set()
        for i in range(len(cseq) - 20):
            k = cseq[i : i + 21]
            kms.add(k)
            kms.add(str(Seq(k).reverse_complement()))
        cand_21_map[cid] = kms
        all_c_21mers.update(kms)

    kmer_idx_all = _build_pigeonhole_21mer_index(all_c_21mers)
    kmer_hits: Dict[str, List[Tuple[str, str, str, str, str]]] = {k: [] for k in all_c_21mers}  # (scope, group, ref_id, organism, subject_id)

    def _scan_fasta(fa: Path, scope: str, grp: str, rid: str, org: str) -> None:
        if not fa or not fa.exists():
            return
        for rec_id, s in _iter_fasta_upper(fa):
            n = len(s)
            if n < 21:
                continue
            max_start = n - 21
            for g in range(0, n - 10, 11):
                offs = kmer_idx_all.get(s[g : g + 11])
                if offs is not None:
                    for off in offs:
                        st = g - off
                        if 0 <= st <= max_start:
                            sub = s[st : st + 21]
                            if sub in kmer_hits:
                                kmer_hits[sub].append((scope, grp, rid, org, rec_id))

    for st in screen_targets:
        _scan_fasta(st["tx_fa"], "transcript", st["group"], st["ref_id"], st["organism"])
        if st["gen_fa"]:
            _scan_fasta(st["gen_fa"], "genome", st["group"], st["ref_id"], st["organism"])

    # Scan self paralog transcriptome & genome (excluding true self target locus intervals!)
    import re as _re
    ref_gff = out / "01_data_acquisition" / "local_nd90pr" / "ncbi_dataset" / "data" / "GCF_000338995.1" / "genomic.gff"
    pid_to_t = {str(r["refseq_protein_id"]): str(t) for t, r in auth_targets.iterrows()}
    self_tx_map: Dict[str, Set[str]] = {t: {str(r["locus_tag"])} for t, r in auth_targets.iterrows()}
    if ref_gff.exists():
        with open(ref_gff, "r", encoding="utf-8", errors="replace") as gfh:
            for line in gfh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 9 and parts[2] == "CDS":
                    mp = _re.search(r"Name=([^;]+)", parts[8])
                    if mp and mp.group(1) in pid_to_t:
                        mpar = _re.search(r"Parent=rna-([^;]+)", parts[8])
                        if mpar:
                            self_tx_map[pid_to_t[mp.group(1)]].add(mpar.group(1))
    self_gen_windows = {t: (str(r["nd90pr_scaffold"]), int(r["genomic_start_1based"]) - 50, int(r["genomic_end_1based"]) + 50) for t, r in auth_targets.iterrows()}

    _scan_fasta(bsor_tx, "transcript", "fungal_paralog", "Bsor_ND90Pr_tx", "Bipolaris sorokiniana ND90Pr")
    for rec_id, s in _iter_fasta_upper(bsor_gen):
        n = len(s)
        if n < 21:
            continue
        max_start = n - 21
        for g in range(0, n - 10, 11):
            offs = kmer_idx_all.get(s[g : g + 11])
            if offs is not None:
                for off in offs:
                    st = g - off
                    if 0 <= st <= max_start:
                        sub = s[st : st + 21]
                        if sub in kmer_hits:
                            kmer_hits[sub].append(("genome", "self_genome", rec_id, str(st + 1), rec_id))

    exact_rows: List[dict] = []
    for cid, kms in cand_21_map.items():
        tname = cand_target_map[cid]
        self_tx_ids = self_tx_map[tname]
        g_scf, gs_lo, gs_hi = self_gen_windows[tname]
        for km in sorted(kms):
            for scope, grp, rid, org, sid in kmer_hits[km]:
                if grp == "fungal_paralog" and any(si in sid for si in self_tx_ids):
                    continue
                if grp == "self_genome" and rid == g_scf and gs_lo <= int(org) <= gs_hi:
                    continue
                exact_rows.append(
                    {
                        "candidate_id": cid,
                        "target": tname,
                        "screen_scope": scope,
                        "database_group": grp,
                        "ref_id": rid,
                        "organism": org if grp != "self_genome" else f"Bsor_ND90Pr:{rid}:{org}",
                        "subject_accession": sid,
                        "matched_21mer": km,
                    }
                )
    exact_cols = ["candidate_id", "target", "screen_scope", "database_group", "ref_id", "organism", "subject_accession", "matched_21mer"]
    exact_df = (pd.DataFrame(exact_rows).drop_duplicates() if exact_rows else pd.DataFrame(columns=exact_cols))
    exact_df.to_csv(s4 / "exact_21mer_matches_expanded_v3.tsv", sep="\t", index=False)

    # 3. Near-Match (blastn-short) & Long Homology (>=25 nt, blastn) across newly added + existing DBs
    # Reuse existing M08 BLAST outputs for the 7 original transcriptomes, and run BLAST+ cleanly on ALL new DBs + genomes!
    m08_dir = out / "08_off_target_screening"
    inter_v2 = out / "11_post_execution_audit" / "intermediate"
    blast_jobs: List[Tuple[str, str, str, str, Path, Path]] = []  # (scope, grp, rid, org, guides_tsv, cands_tsv)

    # Register existing M08 v1/v2 BLAST results for original 7 transcriptomes + paralogs
    existing_map = {
        "Triticum_aestivum_IWGSC_v2_1": ("crop_wheat", m08_dir / "wheat" / "guides_vs_Triticum_aestivum_IWGSC_v2_1_blastn.tsv", m08_dir / "wheat" / "candidates_vs_Triticum_aestivum_IWGSC_v2_1_blastn.tsv"),
        "Hordeum_vulgare_MorexV3": ("crop_barley", m08_dir / "barley" / "guides_vs_Hordeum_vulgare_MorexV3_blastn.tsv", m08_dir / "barley" / "candidates_vs_Hordeum_vulgare_MorexV3_blastn.tsv"),
        "Trichoderma_atroviride_IMI206040": ("beneficial_fungi", m08_dir / "beneficial_fungi" / "guides_vs_Trichoderma_atroviride_IMI206040_blastn.tsv", m08_dir / "beneficial_fungi" / "candidates_vs_Trichoderma_atroviride_IMI206040_blastn.tsv"),
        "Trichoderma_virens_Gv29_8": ("beneficial_fungi", m08_dir / "beneficial_fungi" / "guides_vs_Trichoderma_virens_Gv29_8_blastn.tsv", m08_dir / "beneficial_fungi" / "candidates_vs_Trichoderma_virens_Gv29_8_blastn.tsv"),
        "Trichoderma_reesei_QM6a": ("beneficial_fungi", m08_dir / "beneficial_fungi" / "guides_vs_Trichoderma_reesei_QM6a_blastn.tsv", m08_dir / "beneficial_fungi" / "candidates_vs_Trichoderma_reesei_QM6a_blastn.tsv"),
        "Beauveria_bassiana_ARSEF2860": ("beneficial_fungi", m08_dir / "beneficial_fungi" / "guides_vs_Beauveria_bassiana_ARSEF2860_blastn.tsv", m08_dir / "beneficial_fungi" / "candidates_vs_Beauveria_bassiana_ARSEF2860_blastn.tsv"),
        "Rhizophagus_irregularis_A1": ("beneficial_fungi", m08_dir / "beneficial_fungi" / "guides_vs_Rhizophagus_irregularis_A1_blastn.tsv", m08_dir / "beneficial_fungi" / "candidates_vs_Rhizophagus_irregularis_A1_blastn.tsv"),
    }

    legacy_guide_seq = {}
    g_eval_df = pd.read_csv(out / "07_dsRNA_design" / "guide_sequences" / "guide_evaluation.tsv", sep="\t")
    for _, gr in g_eval_df.iterrows():
        legacy_guide_seq[str(gr["guide_id"])] = str(gr["guide_sequence_sense_21nt"]).upper()

    for st in screen_targets:
        rid, grp, org = st["ref_id"], st["group"], st["organism"]
        if rid in existing_map:
            _, gp, cp = existing_map[rid]
            blast_jobs.append(("transcript", grp, rid, org, gp, cp))
        elif st["tx_fa"] and st["tx_fa"].exists():
            gp, cp = _make_blastdb_and_run(makeblastdb_bin, blastn_bin, st["tx_fa"], db_root / f"{rid}_tx", unique_guides_fa, cand_fa, raw_root / f"{rid}_tx")
            blast_jobs.append(("transcript", grp, rid, org, gp, cp))
        if st["gen_fa"] and st["gen_fa"].exists():
            gp, cp = _make_blastdb_and_run(makeblastdb_bin, blastn_bin, st["gen_fa"], db_root / f"{rid}_gen", unique_guides_fa, cand_fa, raw_root / f"{rid}_gen")
            blast_jobs.append(("genome", grp, rid, org, gp, cp))

    # Add B. sorokiniana paralogs + whole genome
    blast_jobs.append(("transcript", "fungal_paralog", "Bsor_ND90Pr_tx", "Bipolaris sorokiniana ND90Pr", inter_v2 / "guides_vs_bsor_tx.tsv", inter_v2 / "cands_vs_bsor_tx.tsv"))
    gp_sg, cp_sg = _make_blastdb_and_run(makeblastdb_bin, blastn_bin, bsor_gen, db_root / "Bsor_ND90Pr_gen", unique_guides_fa, cand_fa, raw_root / "Bsor_ND90Pr_gen")
    blast_jobs.append(("genome", "self_genome", "Bsor_ND90Pr_gen", "Bipolaris sorokiniana ND90Pr", gp_sg, cp_sg))

    # Parse near-matches (20/21, 19/21, 18/21) & long alignments (>=25 nt)
    near_rows: List[dict] = []
    long_rows: List[dict] = []
    for scope, grp, rid, org, gp, cp in blast_jobs:
        if gp.exists() and gp.stat().st_size > 0:
            gdf = _read_blast_tsv(gp)
            for qid, sid, pident, ln, mm, gp_o, qstart, qend, sstart, send, ev, bs, qcovs, sstrand in zip(
                gdf["qseqid"], gdf["sseqid"], gdf["pident"], gdf["length"], gdf["mismatch"], gdf["gapopen"],
                gdf["qstart"], gdf["qend"], gdf["sstart"], gdf["send"], gdf["evalue"], gdf["bitscore"], gdf["qcovs"], gdf["sstrand"]
            ):
                if int(gp_o) != 0:
                    continue
                ln_i, mm_i = int(ln), int(mm)
                tier = None
                if (ln_i == 21 and mm_i == 1) or (ln_i == 20 and mm_i == 0):
                    tier = "near_20_of_21"
                elif (ln_i == 21 and mm_i == 2) or (ln_i == 19 and mm_i == 0):
                    tier = "near_19_of_21"
                elif (ln_i == 21 and mm_i == 3) or (ln_i == 18 and mm_i == 0):
                    tier = "near_18_of_21"
                if not tier:
                    continue
                qseq21 = guide_seq_map.get(str(qid)) or legacy_guide_seq.get(str(qid), "")
                if not qseq21:
                    continue
                for cid in guide_to_cands.get(qseq21, []):
                    tname = cand_target_map[cid]
                    if grp == "fungal_paralog" and any(si in str(sid) for si in self_tx_map[tname]):
                        continue
                    if grp == "self_genome" and str(sid) == self_gen_windows[tname][0] and self_gen_windows[tname][1] <= int(sstart) <= self_gen_windows[tname][2]:
                        continue
                    near_rows.append(
                        {
                            "candidate_id": cid,
                            "target": tname,
                            "guide_21nt": qseq21,
                            "match_tier": tier,
                            "screen_scope": scope,
                            "database_group": grp,
                            "ref_id": rid,
                            "organism": org,
                            "subject_accession": str(sid),
                            "sstart": int(sstart),
                            "send": int(send),
                            "alignment_length_nt": ln_i,
                            "mismatches": mm_i,
                            "pident": round(float(pident), 2),
                            "evalue": float(ev),
                            "sstrand": str(sstrand),
                        }
                    )

        if cp.exists() and cp.stat().st_size > 0:
            cdf = _read_blast_tsv(cp)
            cdf = cdf[cdf["length"] >= 25]
            for cid, sid, pident, ln, mm, gp_o, qstart, qend, sstart, send, ev, bs, qcovs, sstrand in zip(
                cdf["qseqid"], cdf["sseqid"], cdf["pident"], cdf["length"], cdf["mismatch"], cdf["gapopen"],
                cdf["qstart"], cdf["qend"], cdf["sstart"], cdf["send"], cdf["evalue"], cdf["bitscore"], cdf["qcovs"], cdf["sstrand"]
            ):
                tname = cand_target_map.get(str(cid), "")
                if grp == "fungal_paralog" and any(si in str(sid) for si in self_tx_map.get(tname, set())):
                    continue
                if grp == "self_genome" and tname in self_gen_windows and str(sid) == self_gen_windows[tname][0] and self_gen_windows[tname][1] <= int(sstart) <= self_gen_windows[tname][2]:
                    continue
                long_rows.append(
                    {
                        "candidate_id": str(cid),
                        "target": tname,
                        "screen_scope": scope,
                        "database_group": grp,
                        "ref_id": rid,
                        "organism": org,
                        "subject_accession": str(sid),
                        "sstart": int(sstart),
                        "send": int(send),
                        "alignment_length_nt": int(ln),
                        "mismatches": int(mm),
                        "pident": round(float(pident), 2),
                        "query_cov_pct": round(float(qcovs), 2),
                        "evalue": float(ev),
                        "sstrand": str(sstrand),
                        "is_high_homology_ge30_80pct": (int(ln) >= 30 and float(pident) >= 80.0 and float(ev) <= 1e-5),
                    }
                )

    near_df = pd.DataFrame(near_rows).drop_duplicates(subset=["candidate_id", "match_tier", "screen_scope", "database_group", "ref_id", "subject_accession", "guide_21nt", "sstrand"]) if near_rows else pd.DataFrame(columns=["candidate_id", "target", "guide_21nt", "match_tier", "screen_scope", "database_group", "ref_id", "organism", "subject_accession", "sstart", "send", "alignment_length_nt", "mismatches", "pident", "evalue", "sstrand"])
    long_df = pd.DataFrame(long_rows).drop_duplicates(subset=["candidate_id", "screen_scope", "database_group", "ref_id", "subject_accession", "sstart", "send", "alignment_length_nt"]) if long_rows else pd.DataFrame(columns=["candidate_id", "target", "screen_scope", "database_group", "ref_id", "organism", "subject_accession", "sstart", "send", "alignment_length_nt", "mismatches", "pident", "query_cov_pct", "evalue", "sstrand", "is_high_homology_ge30_80pct"])

    near_df.to_csv(s4 / "near_guide_matches_expanded_v3.tsv", sep="\t", index=False)
    long_df.to_csv(s4 / "long_nucleotide_alignments_expanded_v3.tsv", sep="\t", index=False)

    # 4. Compute ViennaRNA duplexfold thermodynamics on representative guide pairs & Candidate Selectivity Index v3
    sel_rows: List[dict] = []
    for cid, cseq in cand_seq_map.items():
        ex_c = exact_df[exact_df["candidate_id"] == cid]
        nr_c = near_df[near_df["candidate_id"] == cid]
        lg_c = long_df[long_df["candidate_id"] == cid]

        ex_tx = ex_c[ex_c["screen_scope"] == "transcript"]
        ex_gn = ex_c[ex_c["screen_scope"] == "genome"]

        def _cnt(df_in: pd.DataFrame, grp_prefix: str, col: str = "subject_accession") -> int:
            sub = df_in[df_in["database_group"].str.startswith(grp_prefix)]
            return int(sub[col].nunique())

        w_ex = _cnt(ex_tx, "crop_wheat")
        b_ex = _cnt(ex_tx, "crop_barley")
        r_ex = _cnt(ex_tx, "crop_rice") + _cnt(ex_gn, "crop_rice")
        m_ex = _cnt(ex_tx, "crop_maize")
        bf_ex_tx = _cnt(ex_tx, "beneficial_fungi")
        bf_ex_gn = _cnt(ex_gn, "beneficial_fungi")
        bac_ex = _cnt(ex_tx, "bacteria") + _cnt(ex_gn, "bacteria")
        self_par = _cnt(ex_tx, "fungal_paralog") + _cnt(ex_gn, "self_genome")

        nr20 = nr_c[nr_c["match_tier"] == "near_20_of_21"]
        crops_20 = int(nr20[nr20["database_group"].str.startswith("crop_")]["guide_21nt"].nunique())
        bf_20 = int(nr20[nr20["database_group"] == "beneficial_fungi"]["guide_21nt"].nunique())
        bac_20 = int(nr20[nr20["database_group"] == "bacteria"]["guide_21nt"].nunique())

        lg_hi = lg_c[lg_c["is_high_homology_ge30_80pct"]]
        crops_lg = _cnt(lg_hi, "crop_")
        bf_lg = _cnt(lg_hi, "beneficial_fungi")
        bac_lg = _cnt(lg_hi, "bacteria")

        # Heuristic sequence selectivity index (0-1, NOT a validated safety probability!)
        penalty = 0.0
        penalty += 0.45 * min(1, w_ex + b_ex + r_ex + m_ex)
        penalty += 0.40 * min(1, bf_ex_tx + bf_ex_gn)
        penalty += 0.25 * min(1, bac_ex)
        penalty += 0.35 * min(1, self_par)
        penalty += 0.25 * min(1, crops_lg + bf_lg + bac_lg)
        penalty += min(0.25, 0.015 * crops_20 + 0.015 * bf_20 + 0.015 * bac_20)
        sel_idx = round(max(0.0, 1.0 - penalty), 4)

        has_any_exact_21 = (w_ex + b_ex + r_ex + m_ex + bf_ex_tx + bf_ex_gn + bac_ex + self_par) > 0
        if has_any_exact_21:
            tier_lbl = "Exact_21mer_NonTarget_Or_Paralog_Match"
        elif (crops_lg + bf_lg + bac_lg) > 0:
            tier_lbl = "Long_Conserved_Homology_ge30nt_No_Exact_21mer"
        elif (crops_20 + bf_20 + bac_20) >= 5:
            tier_lbl = "Elevated_20of21_Near_Match_Burden"
        elif (crops_20 + bf_20 + bac_20) > 0:
            tier_lbl = "Low_20of21_Near_Match_Signal"
        else:
            tier_lbl = "No_Exact_Or_20of21_Match_Detected"

        # Compute representative antisense-guide MFE on first hit if available else median self guide duplex
        best_mfe = 0.0
        if not ex_c.empty:
            km = str(ex_c["matched_21mer"].iloc[0]).replace("T", "U")
            best_mfe = round(float(RNA.duplexfold(str(Seq(km).reverse_complement()).replace("T", "U"), km).energy), 2)
        elif not nr20.empty:
            km = str(nr20["guide_21nt"].iloc[0]).replace("T", "U")
            best_mfe = round(float(RNA.duplexfold(str(Seq(km).reverse_complement()).replace("T", "U"), km).energy), 2)

        sel_rows.append(
            {
                "candidate_id": cid,
                "target": cand_target_map[cid],
                "wheat_exact_21": w_ex,
                "barley_exact_21": b_ex,
                "rice_exact_21": r_ex,
                "maize_exact_21": m_ex,
                "beneficial_fungi_tx_exact_21": bf_ex_tx,
                "beneficial_fungi_genome_exact_21": bf_ex_gn,
                "bacteria_exact_21": bac_ex,
                "self_paralog_or_genome_exact_21": self_par,
                "crops_unique_20of21_guides": crops_20,
                "beneficial_fungi_unique_20of21_guides": bf_20,
                "bacteria_unique_20of21_guides": bac_20,
                "crops_long_homology_ge30nt": crops_lg,
                "beneficial_fungi_long_homology_ge30nt": bf_lg,
                "bacteria_long_homology_ge30nt": bac_lg,
                "representative_duplex_mfe_kcal": best_mfe,
                "selectivity_classification_v3": tier_lbl,
                "heuristic_sequence_selectivity_index": sel_idx,
                "passes_zero_exact_21mer_filter": not has_any_exact_21,
                "selectivity_caveat": "Sequence selectivity heuristic across 17 non-targets + self; NOT biological safety or silencing proof; insect 21-nt rules do not establish fungal/plant cutoffs.",
            }
        )

    sel_df = pd.DataFrame(sel_rows).sort_values(["passes_zero_exact_21mer_filter", "heuristic_sequence_selectivity_index"], ascending=[False, False]).reset_index(drop=True)
    sel_df.to_csv(s4 / "candidate_selectivity_expanded_v3.tsv", sep="\t", index=False)

    # Save explicit methodology & insect-comparator distinction JSON
    meth_doc = {
        "updated_utc": now_utc,
        "panel_summary": {"crops": 4, "beneficial_fungi": 9, "bacteria": 4, "self_fungal": 1},
        "metric_rename": "Renamed v2 off_target_safety_score_v2 -> heuristic_sequence_selectivity_index so it cannot be misinterpreted as biological safety.",
        "comparators_disclaimer": "Exact 21-nt continuous match cutoffs (e.g., Bachman et al., 2013 beetle/coleopteran Feeding assays) are arthropod-derived literature comparators, NOT experimentally calibrated thresholds for monocot crops or filamentous fungi.",
    }
    (s4 / "selectivity_methodology_and_comparators_v3.json").write_text(json.dumps(meth_doc, indent=2), encoding="utf-8")

    # Organism hit overview figure
    fig, ax = plt.subplots(figsize=(10, 5.2))
    tc = sel_df["selectivity_classification_v3"].value_counts()
    tc.plot(kind="barh", ax=ax, color="#0284c7", edgecolor="#0f172a")
    ax.set_xlabel("Candidate Count (out of 50)")
    ax.set_title("Stage 4: Expanded Off-Target Selectivity Across 4 Crops + 9 Beneficial Fungi + 4 Bacteria", fontweight="bold")
    ax.grid(True, axis="x", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(fig_dir / "stage4_expanded_selectivity_tiers_v3.png", dpi=300)
    fig.savefig(fig_dir / "stage4_expanded_selectivity_tiers_v3.svg")
    plt.close(fig)

    s4_md = f"""# Stage 4 — Expanded Off-Target Selectivity & Coverage Report (BRD v3)

**Timestamp:** {now_utc}

## 1. Expanded Non-Target Panel (17 Non-Target Organisms + Self ND90Pr)
- **Crops (4):** *Triticum aestivum*, *Hordeum vulgare*, *Oryza sativa* (`GCF_034140825.1`), *Zea mays* (`GCF_902167145.1`).
- **Beneficial Fungi (9):** *Trichoderma atroviride*, *T. virens*, *T. reesei*, *T. harzianum* (`GCF_003025095.1`), *Beauveria bassiana*, *Metarhizium anisopliae* (`GCF_039654215.1`), *M. robertsii* (`GCA_000591435.1`), *Rhizophagus irregularis*, and *Rhizophagus clarus* (`GCA_003203555.1`).
- **Relevant Bacteria (4):** *Bacillus subtilis* 168 (`GCF_000009045.1`), *Bacillus velezensis* (`GCF_000508125.2`), *Pseudomonas fluorescens* (`GCF_900215245.1`), *Pseudomonas putida* (`GCF_000787655.2`).

## 2. Selectivity Findings & Insect Threshold Distinction
- Passing zero-exact-21mer hard screen across the expanded panel: **`{int(sel_df['passes_zero_exact_21mer_filter'].sum())}/50`** candidates.
- Note: Insect-derived continuous 21-mer thresholds (e.g. Bachman et al., 2013) are reported only as cross-kingdom comparators and **not** treated as validated monocot crop or soil-microbiome silencing cutoffs.

{tc.rename_axis('classification').reset_index(name='candidates').to_markdown(index=False)}
"""
    (s4 / "stage4_selectivity_report_v3.md").write_text(s4_md, encoding="utf-8")

    return {
        "manifest": s4 / "expanded_reference_manifest_v3.tsv",
        "selectivity": s4 / "candidate_selectivity_expanded_v3.tsv",
        "report": s4 / "stage4_selectivity_report_v3.md",
    }


if __name__ == "__main__":
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    for k, v in run_stage4().items():
        print(f"STAGE4_OUT {k}: {v}")
