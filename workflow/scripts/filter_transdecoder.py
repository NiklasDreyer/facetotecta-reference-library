#!/usr/bin/env python3
# =============================================================================
# filter_transdecoder.py: keep COX1 ASVs with an (almost) full-length ORF.
#
# TransDecoder.LongOrfs writes the longest ORFs of each ASV to longest_orfs.pep.
# Headers contain e.g. "len:104" and "ILL-ASV1:2-313(+)". An ASV passes
# when its longest ORF on the + strand covers at least --min-fraction of the
# ASV length (in codons). A stop codon or frameshift inside the ASV breaks
# the ORF, so such ASVs fail (manuscript_v1: screen for stop codons and
# premature gaps).
# =============================================================================

import argparse
import re
import sys

from yfutils import read_fasta, write_fasta, write_metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--pep", required=True, help="TransDecoder longest_orfs.pep")
    ap.add_argument("--min-fraction", type=float, default=0.95)
    ap.add_argument("--out", required=True)
    ap.add_argument("--metrics")
    args = ap.parse_args()

    seqs = {n: s for n, _, s in read_fasta(args.fasta)}

    # Longest ORF length (in nucleotides) per ASV, from the coordinates.
    longest = {}
    for _, desc, _ in read_fasta(args.pep):
        m = re.search(r"(\S+):(\d+)-(\d+)\(([+-])\)", desc)
        if not m:
            continue
        asv, a, b, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
        if strand != "+":
            continue
        span = abs(b - a) + 1
        longest[asv] = max(longest.get(asv, 0), span)

    kept = []
    for name, seq in seqs.items():
        # Fraction of the ASV covered by its longest open reading frame.
        frac = longest.get(name, 0) / max(len(seq) - len(seq) % 3, 1)
        if frac >= args.min_fraction:
            kept.append((name, seq))
    write_fasta(kept, args.out)
    print(f"{len(kept)} of {len(seqs)} ASVs pass the translation screen", file=sys.stderr)
    write_metrics(args.metrics, {
        "illumina_COX1_transdecoder_input": len(seqs),
        "illumina_COX1_transdecoder_kept": len(kept),
    })


if __name__ == "__main__":
    main()
