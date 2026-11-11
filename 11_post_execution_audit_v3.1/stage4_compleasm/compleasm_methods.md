# Stage 4 — Compleasm Genome Completeness Methodology & Blocker Record (BRD v3.1)

**Evaluated:** 2026-11-11 08:18:46 UTC  
**Compleasm Version:** `v0.2.10` (`Apache-2.0`, `tools/repos/compleasm/compleasm.py`)  
**Lineages Downloaded & Checksum-Verified Locally:**
1. **`pleosporales_odb12`** (`OrthoDB v12.1`, `6,066` BUSCO HMM markers across `76` Pleosporales species; MD5 `795a5af234fb9b7828b9927dfbebff4a`) — selected as primary taxonomically specific lineage because *Bipolaris sorokiniana* belongs to Order *Pleosporales* (`taxid 92860`).
2. **`fungi_odb12`** (`OrthoDB v12.1`, `1,122` BUSCO HMM markers across `144` fungal species; MD5 `79aba415d64ed9a7b1994452b99106cc`) — Kingdom-level comparator dataset.

## Empirical Windows Execution Probe Result
- `compleasm.py -h`, `compleasm.py list --remote`, and `compleasm.py download` executed cleanly (`exit 0`) under native Windows Python `3.13.14`.
- **`compleasm.py run`** was directly tested on `ND90Pr` (`GCF_000338995.1_Cocsa1_genomic.fna`) and failed (`exit 1`) because **`miniprot` (Linux x64 ELF release only) and `hmmsearch` (HMMER `3.4`, POSIX `<unistd.h>`/`./configure` only)** are unavailable natively on Windows, Docker Desktop is not installed, and WSL2 has no Linux distribution installed.
- **No scores fabricated:** `compleasm_summary_v3_1.tsv` explicitly records `NOT_RUN_DEPENDENCY_BLOCKED` and `NA` completeness across all `19` isolates.
- **Retained separate alternative (`compleasm_marker_details_v3_1.tsv`):** Preserves the v3 50-marker tBLASTn frame integrity screen cleanly labelled as `v3_50_marker_tblastn_frame_audit` (NOT BUSCO and NOT Compleasm), which only evaluates `50` conserved genes for frameshift stop artifacts (`WAI2406`, `WAI2411`) and does not replace proteome-wide `6,066`-marker HMM completeness.
