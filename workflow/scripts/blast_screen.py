#!/usr/bin/env python3
# =============================================================================
# blast_screen.py: keep ASVs whose top BLAST hit belongs to the target taxa.
#
# Contaminant screen of manuscript_v1: ASVs are searched against MIDORI (COX1)
# or SILVA (18S) plus our Sanger library. An ASV is kept when the lineage of
# its best hit (highest bitscore) contains one of the --keep taxa
# (default: Facetotecta, Arthropoda). ASVs without any hit are removed.
#
# Inputs: ASV FASTA, BLAST tabular output
# (qseqid sseqid pident length qstart qend sstart send evalue bitscore),
# taxonomy table of the database.
# Outputs: kept ASVs, a per-ASV report and metrics.
# =============================================================================

import argparse
from collections import defaultdict

from yfutils import read_fasta, read_taxonomy, write_fasta, write_metrics


def best_hits(path):
    """Best hit (highest bitscore) per query from BLAST tabular output."""
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            q, s, pid, bits = f[0], f[1], float(f[2]), float(f[9])
            s = s.replace("lcl|", "")  # some BLAST versions prefix local IDs
            if q not in best or bits > best[q][2]:
                best[q] = (s, pid, bits, float(f[8]))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--tax", required=True)
    ap.add_argument("--keep", nargs="+", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--metrics")
    ap.add_argument("--label", default="screen")
    args = ap.parse_args()

    tax = read_taxonomy(args.tax)
    best = best_hits(args.hits)
    kept = []
    category = defaultdict(int)

    with open(args.report, "w") as rep:
        rep.write("asv\tbest_hit\tidentity\tevalue\tlineage\tkept\n")
        for name, _, seq in read_fasta(args.fasta):
            if name not in best:
                rep.write(f"{name}\tNA\tNA\tNA\tno_hit\tno\n")
                category["no_hit"] += 1
                continue
            s, pid, _, ev = best[name]
            lineage = tax.get(s, "unknown")
            ok = any(k in lineage.split(";") for k in args.keep)
            # Record which kept taxon was matched (the most specific wins).
            if "Facetotecta" in lineage:
                category["Facetotecta"] += 1
            elif ok:
                category["other_kept"] += 1
            else:
                category["removed"] += 1
            rep.write(f"{name}\t{s}\t{pid}\t{ev}\t{lineage}\t{'yes' if ok else 'no'}\n")
            if ok:
                kept.append((name, seq))

    write_fasta(kept, args.out_fasta)
    write_metrics(args.metrics, {
        f"{args.label}_screen_kept": len(kept),
        f"{args.label}_screen_tophit_facetotecta": category["Facetotecta"],
        f"{args.label}_screen_tophit_other_kept": category["other_kept"],
        f"{args.label}_screen_removed": category["removed"],
        f"{args.label}_screen_no_hit": category["no_hit"],
    })


if __name__ == "__main__":
    main()
