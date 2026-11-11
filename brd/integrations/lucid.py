"""LUCID (sifi21_plus + ViennaRNA accessibility + Primer3) baseline integration."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import primer3
import RNA
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq

# Compatibility shim for sifi21_plus on Biopython >= 1.82 where Bio.SeqUtils.GC was removed
if not hasattr(SeqUtils, "GC") and hasattr(SeqUtils, "gc_fraction"):
    SeqUtils.GC = lambda seq: SeqUtils.gc_fraction(seq) * 100.0

import sifi21_plus.free_energy as free_energy


def compute_unpaired_probs(dna_seq: str, temperature: float = 22.0, window: int = 8) -> np.ndarray:
    """
    Compute per-nucleotide unpaired probability profile using ViennaRNA 2.8.0 partition function
    (sliding 80-nt chunks for long transcripts like CsPKS1 to mirror RNAplfold -W 80 -L 40 at 22C).
    """
    rna_seq = dna_seq.upper().replace("T", "U")
    n = len(rna_seq)
    md = RNA.md()
    md.temperature = temperature

    if n <= 350:
        fc = RNA.fold_compound(rna_seq, md)
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        bpp = fc.bpp()
        unpaired = np.ones(n, dtype=float)
        for i in range(1, n + 1):
            paired_prob = sum(bpp[i][j] for j in range(i + 1, n + 1)) + sum(bpp[j][i] for j in range(1, i))
            unpaired[i - 1] = max(0.0, min(1.0, 1.0 - paired_prob))
        return unpaired

    # Sliding overlapping 160-nt windows (step 40) for long transcripts (e.g., CsPKS1 6462 nt)
    sums = np.zeros(n, dtype=float)
    counts = np.zeros(n, dtype=float)
    wlen = 160
    step = 40
    for start in range(0, n, step):
        end = min(n, start + wlen)
        sub = rna_seq[start:end]
        if len(sub) < 20:
            break
        fc = RNA.fold_compound(sub, md)
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        bpp = fc.bpp()
        m = len(sub)
        for i in range(1, m + 1):
            paired_prob = sum(bpp[i][j] for j in range(i + 1, m + 1)) + sum(bpp[j][i] for j in range(1, i))
            sums[start + i - 1] += max(0.0, min(1.0, 1.0 - paired_prob))
            counts[start + i - 1] += 1.0
        if end == n:
            break
    counts[counts == 0] = 1.0
    return sums / counts


def score_sirnas_lucid_sifi21(
    dna_seq: str,
    sirna_size: int = 21,
    end_stability_threshold: float = 1.0,
    accessibility_threshold: float = 0.1,
    accessibility_window: int = 8,
) -> pd.DataFrame:
    """
    Reproduce sifi21_plus efficiency criteria from Brunxi/LUCID:
      1. Strand selection (antisense 5' MFE >= sense 5' MFE)
      2. End stability differential (antisense5 - sense5 >= 1.0)
      3. Terminal S1/S19 base asymmetry rules (lines 528-552 of sifiPipeline.py)
      4. Target-site accessibility (ViennaRNA unpaired probability >= 0.1)
    """
    seq = dna_seq.upper()
    n = len(seq)
    unpaired = compute_unpaired_probs(seq, temperature=22.0, window=accessibility_window)
    rows = []
    sirna_list = [seq[i : i + sirna_size] for i in range(0, n - sirna_size + 1)]

    for pos, s_seq in enumerate(sirna_list):
        if len(s_seq) < sirna_size or any(c not in "ACGT" for c in s_seq):
            continue
        if pos >= 2:
            s_n2 = sirna_list[pos - 2]
            sense_5 = s_seq[0:3]
            sense_c = str(Seq(s_n2).reverse_complement())[sirna_size - 5 : sirna_size - 1]
            anti_5 = str(Seq(s_n2).reverse_complement())[0:3]
            anti_c = s_seq[sirna_size - 5 : sirna_size - 1]
            try:
                g_s = free_energy.calculate_free_energy(sense_5, check=True, strict=True, c_seq=sense_c[::-1], shift=1)
                g_as = free_energy.calculate_free_energy(anti_5, check=True, strict=True, c_seq=anti_c[::-1], shift=1)
            except Exception:
                g_s = free_energy.calculate_free_energy(s_seq[0:3])
                g_as = free_energy.calculate_free_energy(s_seq[sirna_size - 5 : sirna_size - 2])
        else:
            g_s = free_energy.calculate_free_energy(s_seq[0:3])
            g_as = free_energy.calculate_free_energy(s_seq[sirna_size - 5 : sirna_size - 2])

        strand_sel = g_as >= g_s
        end_stab = (g_as - g_s) >= end_stability_threshold
        acc_idx = min(n - 1, pos + accessibility_window)
        acc_val = float(unpaired[acc_idx])
        acc_ok = acc_val >= accessibility_threshold
        thermo_eff = strand_sel and end_stab and acc_ok

        # Exact sifiPipeline.py terminal_check logic
        s19 = s_seq[sirna_size - 3]
        s1 = s_seq[1]
        if s19 in ("A", "T"):
            if s1 in ("A", "T"):
                is_eff = thermo_eff
            else:
                is_eff = acc_ok
        else:
            if s1 in ("G", "C"):
                is_eff = thermo_eff
            else:
                is_eff = False

        rows.append(
            {
                "sirna_pos_1based": pos + 1,
                "sirna_seq": s_seq,
                "sense5_mfe": round(float(g_s), 3),
                "antisense5_mfe": round(float(g_as), 3),
                "delta_end_stability": round(float(g_as - g_s), 3),
                "strand_selection": strand_sel,
                "end_stability_ok": end_stab,
                "accessibility_prob": round(acc_val, 4),
                "accessibility_ok": acc_ok,
                "is_efficient_sifi21": is_eff,
            }
        )
    return pd.DataFrame(rows)


def design_lucid_windows(
    target_name: str,
    dna_seq: str,
    min_len: int = 200,
    max_len: int = 400,
    step_len: int = 50,
    step_pos: int = 25,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run LUCID/sifi21_plus siRNA efficiency profiling and select optimal accessible windows + Primer3 primers."""
    sirna_df = score_sirnas_lucid_sifi21(dna_seq)
    seq_len = len(dna_seq)
    eff_vec = np.zeros(seq_len, dtype=float)
    acc_vec = np.zeros(seq_len, dtype=float)
    for _, r in sirna_df.iterrows():
        p0 = int(r["sirna_pos_1based"]) - 1
        acc_vec[p0] = float(r["accessibility_prob"])
        if r["is_efficient_sifi21"]:
            eff_vec[p0] = 1.0

    # Handle targets shorter than max_len (e.g., CsSp1 270 nt)
    usable_max = min(max_len, seq_len)
    lengths = [l for l in range(min_len, usable_max + 1, step_len)]
    if not lengths and seq_len >= 150:
        lengths = [min(seq_len, max(150, seq_len - 20)), seq_len]

    win_rows = []
    for wlen in sorted(set(lengths)):
        for start in range(0, seq_len - wlen + 1, step_pos):
            end = start + wlen
            sub_eff = sirna_df[(sirna_df["sirna_pos_1based"] >= start + 1) & (sirna_df["sirna_pos_1based"] + 20 <= end)]
            n_tot = len(sub_eff)
            n_eff = int(sub_eff["is_efficient_sifi21"].sum()) if n_tot else 0
            mean_acc = float(sub_eff["accessibility_prob"].mean()) if n_tot else 0.0
            eff_frac = n_eff / max(n_tot, 1)
            subseq = dna_seq[start:end].upper()
            gc = (subseq.count("G") + subseq.count("C")) / len(subseq) * 100.0

            # Primer3 amplification design (LUCID PHASE2_dsrna.sh settings adapted for window length)
            p_min = max(150, wlen - 30)
            p_max = wlen
            primers = primer3.bindings.design_primers(
                seq_args={"SEQUENCE_ID": f"{target_name}_{start+1}_{end}", "SEQUENCE_TEMPLATE": subseq},
                global_args={
                    "PRIMER_OPT_SIZE": 20,
                    "PRIMER_MIN_SIZE": 18,
                    "PRIMER_MAX_SIZE": 25,
                    "PRIMER_OPT_TM": 60.0,
                    "PRIMER_MIN_TM": 57.0,
                    "PRIMER_MAX_TM": 63.0,
                    "PRIMER_MIN_GC": 20.0,
                    "PRIMER_MAX_GC": 80.0,
                    "PRIMER_PRODUCT_SIZE_RANGE": [[p_min, p_max]],
                    "PRIMER_NUM_RETURN": 1,
                },
            )
            lp = primers.get("PRIMER_LEFT_0_SEQUENCE", "")
            rp = primers.get("PRIMER_RIGHT_0_SEQUENCE", "")
            psize = primers.get("PRIMER_PAIR_0_PRODUCT_SIZE", "")

            win_rows.append(
                {
                    "target": target_name,
                    "tool": "LUCID_sifi21_plus_ViennaRNA",
                    "start_1based": start + 1,
                    "end_1based": end,
                    "length_nt": wlen,
                    "gc_percent": round(gc, 2),
                    "total_21mers": n_tot,
                    "efficient_sirnas": n_eff,
                    "efficient_sirna_fraction": round(eff_frac, 4),
                    "mean_accessibility": round(mean_acc, 4),
                    "primer_left": lp or "NA",
                    "primer_right": rp or "NA",
                    "expected_amplicon_bp": psize or "",
                    "sequence": subseq,
                }
            )
    return pd.DataFrame(win_rows), sirna_df
