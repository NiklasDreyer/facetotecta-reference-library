#!/usr/bin/env python3
# =============================================================================
# primer_scan.py: NEW in-silico test of a primer pair against sequences.
#
# For every sequence:
#   - best binding site of the forward primer (fewest mismatches) on the
#     sequence, and of the reverse primer (searched as its reverse complement)
#   - number of mismatches, and mismatches in the 3'-terminal bases
#   - "amplifiable" = both primers with <= --max-mismatch mismatches, no
#     mismatch in the --three-prime 3'-terminal bases, forward site upstream
#     of the reverse site
#   - the amplicon (region between the primers) is written whenever both
#     sites are found in the right order with <= max-mismatch + 4 mismatches,
#     so the benchmark can test short fragments even for poorly matching taxa
# IUPAC codes are handled in primers and sequences; inosine ("I") and N
# match anything. Sequences are assumed to be in sense orientation.
#
# A sequence too short to contain a primer site cannot be scored; its row
# says "site_not_covered". Short Sanger fragments may not span the primer
# sites of every pair, which is why PacBio reads and the Type AE* COX1 are
# included in the reference set.
# =============================================================================

import argparse

import numpy as np

from yfutils import read_fasta, revcomp, write_fasta, write_metrics

# Each base is a 4-bit mask (A=1, C=2, G=4, T=8). Two bases match if their
# masks share a bit.
IUPAC = {"A": 1, "C": 2, "G": 4, "T": 8, "U": 8, "R": 5, "Y": 10, "S": 6, "W": 9,
         "K": 12, "M": 3, "B": 14, "D": 13, "H": 11, "V": 7, "N": 15, "I": 15}


def encode(seq):
    """DNA string -> numpy array of 4-bit masks (gaps and unknowns -> 0)."""
    return np.array([IUPAC.get(c, 0) for c in seq.upper()], dtype=np.uint8)


def best_site(seq_codes, primer_codes):
    """Return (position, mismatch vector) of the best site, or (None, None)."""
    k = len(primer_codes)
    if len(seq_codes) < k:
        return None, None
    # All windows of length k, compared with the primer at once.
    windows = np.lib.stride_tricks.sliding_window_view(seq_codes, k)
    mism = (windows & primer_codes) == 0
    counts = mism.sum(axis=1)
    pos = int(counts.argmin())
    return pos, mism[pos]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--forward", required=True)
    ap.add_argument("--reverse", required=True)
    ap.add_argument("--max-mismatch", type=int, default=3)
    ap.add_argument("--three-prime", type=int, default=3)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-amplicons", required=True)
    ap.add_argument("--metrics")
    ap.add_argument("--label", default="primer")
    args = ap.parse_args()

    fwd = encode(args.forward)
    # The reverse primer binds the sense strand as its reverse complement;
    # its 3' end is then the FIRST bases of the reverse-complement site.
    rev_rc = encode(revcomp(args.reverse.upper().replace("I", "N")))
    kf, kr, k3 = len(fwd), len(rev_rc), args.three_prime

    amplicons = []
    n = {"total": 0, "scored": 0, "amplifiable": 0}
    with open(args.out_table, "w") as out:
        out.write("id\tfwd_pos\tfwd_mm\tfwd_3prime_mm\tfwd_mm_positions\t"
                  "rev_pos\trev_mm\trev_3prime_mm\trev_mm_positions\tstatus\n")
        for name, _, seq in read_fasta(args.input):
            n["total"] += 1
            s = encode(seq.replace("-", ""))
            fp, fm = best_site(s, fwd)
            rp, rm = best_site(s, rev_rc)
            if fp is None or rp is None:
                out.write(f"{name}\tNA\tNA\tNA\t\tNA\tNA\tNA\t\tsite_not_covered\n")
                continue
            n["scored"] += 1
            f_mm, r_mm = int(fm.sum()), int(rm.sum())
            f3 = int(fm[-k3:].sum())          # forward primer 3' end = last bases
            r3 = int(rm[:k3].sum())           # reverse primer 3' end = first bases of its RC
            # Mismatch positions numbered 5'->3' along each primer.
            f_pos = ",".join(str(i + 1) for i in np.where(fm)[0])
            r_pos = ",".join(str(kr - i) for i in np.where(rm)[0])
            ordered = fp + kf <= rp
            ok = ordered and f_mm <= args.max_mismatch and r_mm <= args.max_mismatch and f3 == 0 and r3 == 0
            status = "amplifiable" if ok else ("wrong_order" if not ordered else "mismatch")
            n["amplifiable"] += ok
            out.write(f"{name}\t{fp + 1}\t{f_mm}\t{f3}\t{f_pos}\t{rp + 1}\t{r_mm}\t{r3}\t{r_pos}\t{status}\n")
            # Amplicon for the benchmark (between the primers).
            if ordered and f_mm <= args.max_mismatch + 4 and r_mm <= args.max_mismatch + 4:
                amp = seq.replace("-", "")[fp + kf: rp]
                if len(amp) >= 50:
                    amplicons.append((name, amp))

    write_fasta(amplicons, args.out_amplicons)
    write_metrics(args.metrics, {
        f"{args.label}_sequences": n["total"],
        f"{args.label}_sites_covered": n["scored"],
        f"{args.label}_amplifiable": n["amplifiable"],
        f"{args.label}_amplifiable_fraction_of_covered":
            round(n["amplifiable"] / n["scored"], 4) if n["scored"] else "NA",
        f"{args.label}_amplicons_extracted": len(amplicons),
    })


if __name__ == "__main__":
    main()
