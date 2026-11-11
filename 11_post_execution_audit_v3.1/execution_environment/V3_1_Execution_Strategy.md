# Stage 2 — Windows 11 Execution Strategy & Route Verification (BRD v3.1)

**Verified:** 2026-11-11 08:21:34 UTC

| route                     | availability                       | executed_tools                                                                                                                                                                             | blocked_tools                                                                                                             |
|:--------------------------|:-----------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------|
| Route_A_Native_Windows_11 | ACTIVE                             | Python 3.13.14, Biopython 1.86, ViennaRNA 2.8.0, BLAST+ 2.17.0+, DIAMOND 2.2.8, MAFFT, NCBI datasets 18.41.0, Broccoli v2.0.1 (-phylogenies nj), compleasm.py v0.2.10 (list/download only) | compleasm run (missing Win32 miniprot/hmmsearch), Broccoli bionj (FastTree blocked), QUAST, BUSCO, OrthoFinder, RNAhybrid |
| Route_B_Docker_Desktop    | UNAVAILABLE (docker CLI missing)   | None (Dockerfile + docker-compose.yml prepared cleanly; explicitly unexecuted)                                                                                                             | All container targets until Docker Desktop is installed                                                                   |
| Route_C_Remote_Linux_SSH  | UNCONFIGURED (no remote host/keys) | None (remote_linux_runner.py + config template prepared; dry-run/blocked guard verified)                                                                                                   | All remote Linux targets until SSH endpoint is supplied                                                                   |

## Enforcement Policy
- **Route A** was executed natively whenever upstream binaries or Rust/Python builds succeeded (`Biopython 1.86` assembly QC, `Broccoli v2.0.1` NJ orthology, `compleasm.py v0.2.10` dataset acquisition).
- **Route B** (`execution_environment/docker/`) and **Route C** (`execution_environment/remote_linux/`) are ready-to-run reproducible blueprints and are **never claimed as executed** on this Windows 11 workstation.
