# Stage 1 — Empirical Tool Compatibility Matrix (Tested on Host)

**Verified:** 2026-11-11 08:21:34 UTC | **Host OS:** Windows 11 (`x86_64`, Python `3.13.14`)

| Component | Route A: Native Windows | Route B: Docker (`linux/amd64`) | Route C: Remote Linux SSH | Empirically Verified Status |
|---|---|---|---|---|
| **Biopython `1.86`** | **PASS** (wheel + `SeqIO` + `Align`) | Supported in `Dockerfile` | Supported | Executed cleanly across all 19 assemblies (`exit 0`) |
| **Broccoli `v2.0.1` (`-phylogenies nj`)** | **PASS** (Rust `1.99.1` GNU build + `diamond.exe`) | Supported in `Dockerfile` | Supported | Executed smoke (`4.07s`) + full 36,714-protein run (`255.98s`, `exit 0`) |
| **Broccoli `v2.0.1` (`-phylogenies bionj` + `FastTree`)** | **BLOCKED** (`MicrobesOnline` upstream `FastTree.exe`/`FastTree.c` gated by anti-bot challenge; no `cc1` compiler without full MinGW/MSVC) | Supported in `Dockerfile` (`apt fasttree`) | Supported | Not executed natively; NJ mode validated instead |
| **Compleasm `v0.2.10` CLI (`-h`, `list`, `download`)** | **PASS** (Python `3.13.14` downloaded & verified `pleosporales_odb12` & `fungi_odb12`) | Supported | Supported | Executed `exit 0` |
| **Compleasm `v0.2.10` (`run` / `analyze` / `protein`)** | **FAIL (`exit 1`)** — requires `miniprot` (`x64-linux` only, `-lz -lpthread`) & `hmmsearch` (HMMER 3.4 POSIX `<unistd.h>` / `./configure`) | Defined in `Dockerfile` (**BLOCKED**: Docker missing on host) | Defined in `remote_linux_runner.py` (**BLOCKED**: no remote host) | **NOT_RUN_DEPENDENCY_BLOCKED** (logged in `compleasm_execution_log.txt`) |
| **QUAST / BUSCO / OrthoFinder / RNAhybrid** | **UNAVAILABLE** (Linux/WSL2 required; `wsl -l -v` reports 0 distros) | Containerizable | Supported | **NOT_RUN_DEPENDENCY_BLOCKED** (never fabricated) |
