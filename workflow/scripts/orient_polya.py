#!/usr/bin/env python3
# =============================================================================
# orient_polya.py: put every CCS read in transcript orientation.
#
# Oligo-dT-primed cDNA reads come in both orientations: with a poly(A) tail at
# the 3' end (sense) or a poly(T) stretch at the 5' end (antisense).
# Antisense reads are reverse-complemented (with their quality string
# reversed). A read counts as antisense when it has a run of >= --min-run
# T's within the first --window bases and no comparable run of A's in the
# last --window bases. Reads with neither are kept unchanged and counted.
# =============================================================================

import argparse
import gzip
import re

from yfutils import read_fastq, revcomp, write_fastq_record


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--min-run", type=int, default=12)
    args = ap.parse_args()

    poly_t = re.compile("T{%d,}" % args.min_run)
    poly_a = re.compile("A{%d,}" % args.min_run)
    counts = {"sense": 0, "antisense_flipped": 0, "no_tail": 0}

    with gzip.open(args.output, "wt") as out:
        for name, _, seq, qual in read_fastq(args.input):
            head, tail = seq[: args.window].upper(), seq[-args.window:].upper()
            has_t5 = bool(poly_t.search(head))
            has_a3 = bool(poly_a.search(tail))
            if has_t5 and not has_a3:
                # Antisense: reverse-complement the sequence, reverse the qualities.
                seq, qual = revcomp(seq), qual[::-1]
                counts["antisense_flipped"] += 1
            elif has_a3:
                counts["sense"] += 1
            else:
                counts["no_tail"] += 1
            write_fastq_record(out, name, seq, qual)

    with open(args.report, "w") as rep:
        for k, v in counts.items():
            rep.write(f"pacbio_orientation_{k}\t{v}\n")


if __name__ == "__main__":
    main()
