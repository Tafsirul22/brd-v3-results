# Stage 9 — Target-Gene Biological Evidence Audit (BRD v3.1)

**Audited:** 2026-11-11 08:30:11 UTC

- **`CsSTE11`** retains leading biological relevance (`0.9921`, `PHI:4778` KO loss of pathogenicity, clean `Broccoli OG_1786` 1:1:1 ortholog).
- **`CsPKS1`** (`0.5796`) preserves the mandatory `PHI:2260` **unaffected pathogenicity** warning: RNAi/KO affects melanin pigmentation only, not host virulence.
- **`CsSp1`** (`0.5500`) preserves permanent exclusion of the invalid 1-read `3A_I` RPKM spike (`1,438.44`) and reports validated `3B/4B/5B` expression (`max = 139.61` RPKM, `mean = 79.30` RPKM across `75` reads).
- **`CsFUS3`** (`0.9921`) retains its full off-target selectivity exclusion due to widespread beneficial-fungi and crop MAPK 21-mer homology (`0/13` passing candidates).
- Broccoli `v2.0.1` improved orthology resolution by separating `CsFUS3` (`OG_1315`) and `CsPKS1` (`OG_4251`) from large multi-paralog superfamilies grouped together by DIAMOND RBH (`RBH_OG000002`, `RBH_OG000102`).
