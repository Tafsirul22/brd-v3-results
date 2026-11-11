# Stage 7 — Literature Extraction Methodology (BRD v3)

**Timestamp:** 2026-11-11 06:55:13 UTC

1. **API Provenance:** Literature retrieved via Europe PMC REST API and NCBI E-utilities / OpenAlex (`stage7_fetch_oa.py`, `stage7_extra_fetch.py`), storing raw query JSONs and 25 full-text JATS XML articles under `stage7_literature_benchmark/fulltexts/`. License attributions are stored in `.licenses/`.
2. **Non-Fabrication Protocol:** Zero sequences were reconstructed from abstracts or memory.
   - `1_Direct_dsRNA_Sequence`: Reserved exclusively for verbatim full-length dsRNA/siRNA strings published inline in primary full text (`n=0` full long-dsRNA constructs found verbatim inside the JATS XML bodies; authors placed raw FASTAs in external publisher spreadsheets — explicitly audited rather than fabricated).
   - `2_Reconstructed_Sequence`: Verified solely when exact printed primers in an article matched fetched NCBI accessions (`PMC7214474` `EH229082.1 GCS_H` and `EH249914.1 RP_S16`), documenting exact coordinate/length mismatches honestly (`525` vs `420` nt; `444` vs `500` nt).
   - `3_Gene_Level_Only`: Experimental studies verified from primary XML Results/Tables where species, gene, construct length/delivery, and silencing/disease outcome are documented without inline nucleotide FASTA.
   - `4_Computational_Only`: Unvalidated BRD v3 *B. sorokiniana* candidates.
   - `5_Ambiguous_or_Conflicting_Outcome`: Verified primary failures, non-responsive high-expression genes, uptake-gated species failures, and disease-worsening constructs (`SS1G_06305 +129%`).
