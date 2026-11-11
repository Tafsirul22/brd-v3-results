# Stage 1 — Verification & Tool Test Plan

| Tool | Help/Version Check | Smoke Dataset Execution | Full Dataset Execution | Result Verification |
|---|---|---|---|---|
| **Biopython** | `import Bio; Bio.__version__` | Synthetic FASTA known N50/L50 + dup/empty/short assertions | All 19 *B. sorokiniana* FASTAs | 16 metrics + warnings + plots |
| **Broccoli v2.0.1** | `broccoli.exe -help` (`rc=0`) | 200 proteins x 3 species (`rc=0`, `4.07s`) | All 36,714 proteins across 3 species (`rc=0`, `255.98s`) | `9,589` OGs, `26,745` pairs, RBH cross-check |
| **Compleasm v0.2.10** | `compleasm.py -v`, `list --remote` (`rc=0`) | Real probe on `ND90Pr` (`rc=1`, missing `miniprot`/`hmmsearch`) | Not run (`NOT_RUN_DEPENDENCY_BLOCKED`) | Lineage tarballs MD5/SHA256 verified + summary parser unit-tested |
