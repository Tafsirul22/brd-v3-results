# Broccoli Proteome Recovery & Validation Methods (BRD v3.1.1)

**Verified:** 2026-11-11 09:10:28 UTC

## 1. Investigation Findings (Exact Historical Provenance)
Inspecting `11_post_execution_audit_v3/src/stage2_genome_qc_comparative.py` (lines 155–180) and `11_post_execution_audit_v3.1/src/stage5_broccoli_orthology.py` (lines 52–125) identified all three source proteomes behind the 36,714-protein Broccoli `v2.0.1` run:
1. **`Bsor_ND90Pr.faa` (`12,210` proteins):** Extracted from NCBI RefSeq `GCF_000338995.1` (`Cocsa1`, `protein.faa`, raw SHA-256 `cacc4be7fa1511474f949a7b0eec2ba863faf625ef3a8d49a3feacbabb79f708`).
2. **`Bsor_Sacc._Shoemaker.faa` (`10,755` proteins):** Extracted from NCBI GenBank `GCA_013416765.1` (`ASM1341676v1`, `protein.faa`, raw SHA-256 `5e4d2cf135d516d6f2ab627a0613d0b13365f32a6c181bfeea0e67ca988944b8`).
3. **`Bcin_ASM83294v1.faa` (`13,749` proteins):** Extracted from pinned upstream `Brunxi/LUCID` commit `37910bd6744b2ce6ce4c9b2582819b130e7ea509` (`data/references/Botrytis_cinerea.ASM83294v1.pep.all.fa`, Ensembl Fungi `ASM83294v1`, raw LF SHA-256 `d897d558115de86c1122058d3b58c5cc386f53c1f78dda18376efec407571da7`).

## 2. Deterministic Reconstruction
Writing `all_proteomes_combined_v3.faa` in sorted species order (`>{species}||{protein_id}\n{seq}\n`, CRLF line endings) and splitting per species with description-free `>{protein_id}\n{SEQ}\n` reproduces all four historical SHA-256 checksums **byte-for-byte (`100%`)**:
- `Bcin_ASM83294v1.faa`: `6fccfbf667b2c63695f8d2dbc603c24e413db53a47acdd2799c218cdd948ffa2` (`13,749` proteins, `7,034,175` B)
- `Bsor_ND90Pr.faa`: `2e8783a00d21843ca48f3b5ffc1fc2cf0b4abe486affb8a3cc887ec2acb54c39` (`12,210` proteins, `5,653,448` B)
- `Bsor_Sacc._Shoemaker.faa`: `583e29570f85cc85797b1b3a6773d6c454dcf9c5ad25d17269e6ccbd1cf740a3` (`10,755` proteins, `5,107,801` B)
- `all_proteomes_combined_v3.faa`: `b45191ae98961fef4f1e1c67d8486093f21cc392fc00864b88600a7763a39fd1` (`36,714` proteins, `18,424,497` B)

## 3. Broccoli Clean-Extraction Reproducibility Boundary
- **Proteome inputs:** `100%` reproducibly recoverable from public NCBI + pinned GitHub endpoints (`hydrate_broccoli_proteomes_v3_1_1.py`).
- **Executable requirement:** Re-running Broccoli requires `broccoli` (`v2.0.1`) + `diamond` (`v2.2.8`). Where installed (e.g., workspace `tools/bin/` or Bioconda), a clean-extraction smoke test runs (`status=PASSED`); without `broccoli`/`diamond` binaries on a fresh system, orthology execution remains dependent on installing those external executables.
