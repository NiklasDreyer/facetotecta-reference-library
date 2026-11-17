#!/usr/bin/env python3
# =============================================================================
# select_blast_hits.py: extract reads (or read regions) with good BLAST hits.
#
# For each read, the best hit (highest bitscore) with e-value <= --evalue is
# taken. If --tax and --require-taxon are given, the best hit's lineage must
# contain that taxon (e.g. Facetotecta): this is the "top hit" rule.
# The aligned region of the read (qstart..qend) is extracted and, if the hit
# is on the minus strand of the subject (sstart > send), reverse-complemented
# so all output sequences share the reference orientation.
# --exclude: skip reads already present in another FASTA (e.g. COX1 reads that
# were already recovered by mapping).
# =============================================================================

import argparse
import sys

from yfutils import read_fasta, read_taxonomy, revcomp, write_fasta, write_metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--evalue", type=float, required=True)
    ap.add_argument("--tax")
    ap.add_argument("--require-taxon")
    ap.add_argument("--exclude")
    ap.add_argument("--out", required=True)
    ap.add_argument("--metrics")
    ap.add_argument("--label", default="blast_select")
    args = ap.parse_args()

    tax = read_taxonomy(args.tax) if args.tax else {}
    exclude = set()
    if args.exclude:
        exclude = {n for n, _, _ in read_fasta(args.exclude)}

    # Best hit per read.
    best = {}
    with open(args.hits) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            q, s = f[0], f[1].replace("lcl|", "")  # some BLAST versions prefix local IDs
            qs, qe, ss, se = int(f[4]), int(f[5]), int(f[6]), int(f[7])
            ev, bits = float(f[8]), float(f[9])
            if ev > args.evalue:
                continue
            if q not in best or bits > best[q]["bits"]:
                best[q] = dict(s=s, qs=qs, qe=qe, minus=ss > se, bits=bits)

    out = []
    n_tax_fail, n_excl = 0, 0
    for name, _, seq in read_fasta(args.reads):
        if name not in best:
            continue
        if name in exclude:
            n_excl += 1
            continue
        h = best[name]
        if args.require_taxon:
            lineage = tax.get(h["s"], "")
            if args.require_taxon not in lineage.split(";"):
                n_tax_fail += 1
                continue
        region = seq[h["qs"] - 1: h["qe"]]
        if h["minus"]:
            region = revcomp(region)
        out.append((name, region))

    write_fasta(out, args.out)
    print(f"{len(best)} reads with hits; {n_tax_fail} top hit not {args.require_taxon}; "
          f"{n_excl} excluded; {len(out)} written", file=sys.stderr)
    write_metrics(args.metrics, {
        f"{args.label}_reads_with_hit": len(best),
        f"{args.label}_tophit_other_taxon": n_tax_fail,
        f"{args.label}_excluded_already_found": n_excl,
        f"{args.label}_reads_kept": len(out),
    })


if __name__ == "__main__":
    main()
