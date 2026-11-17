#!/usr/bin/env python3
# =============================================================================
# prefix_headers.py: tidy FASTA headers and add a platform tag.
#
# - keeps only the specimen/sequence ID (text before the first space or "|")
# - replaces characters that tree programs rewrite (anything other than
#   letters, digits, ".", "_") with "_"
# - adds "<PREFIX>-" (e.g. SAN-Y123)
# - removes alignment gaps and makes the sequence upper case
# - with --keep-longest: if a specimen appears more than once, keeps the
#   longest sequence (e.g. several 18S fragments of one specimen)
# - writes metrics: number of sequences and total nucleotides
#
# Usage:
#   prefix_headers.py --input in.fasta --output out.fasta --prefix SAN \
#       [--keep-longest] [--metrics m.tsv --metrics-label sanger_COX1]
# =============================================================================

import argparse
import re
import sys

from yfutils import read_fasta, write_fasta, write_metrics


def clean_id(desc):
    """Specimen ID from a header: first word, cut at '|', unsafe characters -> '_'."""
    first = desc.split()[0].split("|")[0]
    return re.sub(r"[^A-Za-z0-9._]", "_", first)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--prefix", required=True, help="platform tag, e.g. SAN")
    ap.add_argument("--keep-longest", action="store_true")
    ap.add_argument("--metrics")
    ap.add_argument("--metrics-label", default="seqs")
    args = ap.parse_args()

    kept = {}       # id -> sequence
    n_in = 0
    n_dup = 0
    for _, desc, seq in read_fasta(args.input):
        n_in += 1
        sid = f"{args.prefix}-{clean_id(desc)}"
        seq = seq.upper().replace("-", "").replace("?", "")
        if sid in kept:
            n_dup += 1
            if args.keep_longest:
                # Keep whichever sequence is longer.
                if len(seq) > len(kept[sid]):
                    kept[sid] = seq
                continue
            # Without --keep-longest, number the duplicates.
            k = 2
            while f"{sid}.{k}" in kept:
                k += 1
            sid = f"{sid}.{k}"
        kept[sid] = seq

    write_fasta(kept.items(), args.output)
    print(f"{n_in} sequences read, {n_dup} duplicate IDs, {len(kept)} written", file=sys.stderr)

    write_metrics(args.metrics, {
        f"{args.metrics_label}_sequences_input": n_in,
        f"{args.metrics_label}_sequences": len(kept),
        f"{args.metrics_label}_duplicate_ids": n_dup,
        f"{args.metrics_label}_nucleotides": sum(len(s) for s in kept.values()),
        f"{args.metrics_label}_min_length": min((len(s) for s in kept.values()), default=0),
        f"{args.metrics_label}_max_length": max((len(s) for s in kept.values()), default=0),
    })


if __name__ == "__main__":
    main()
