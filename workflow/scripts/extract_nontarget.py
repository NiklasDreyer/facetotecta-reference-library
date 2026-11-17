#!/usr/bin/env python3
# =============================================================================
# extract_nontarget.py: NEW non-target sequences for the assignment benchmark.
#
# Takes sequences from the public reference whose lineage contains one of
# the --groups names as a rank (e.g. Copepoda, Cirripedia), excluding any
# Facetotecta, and samples at most --max-per-group per group (random, with a
# fixed seed). These are the crustaceans most likely to be wrongly assigned
# to Facetotecta.
# Outputs a FASTA and a table: id <tab> group.
# =============================================================================

import argparse
import random
import sys
from collections import defaultdict

from yfutils import read_fasta, read_taxonomy, write_fasta


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--tax", required=True)
    ap.add_argument("--groups", nargs="+", required=True)
    ap.add_argument("--max-per-group", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-groups", required=True)
    args = ap.parse_args()

    tax = read_taxonomy(args.tax)
    by_group = defaultdict(list)
    for sid, lineage in tax.items():
        ranks = lineage.split(";")
        if "Facetotecta" in ranks:
            continue
        for g in args.groups:
            if g in ranks:
                by_group[g].append(sid)
                break

    rng = random.Random(args.seed)
    chosen = {}
    for g, ids in by_group.items():
        pick = ids if len(ids) <= args.max_per_group else rng.sample(ids, args.max_per_group)
        for sid in pick:
            chosen[sid] = g
        print(f"{g}: {len(ids)} available, {len(pick)} sampled", file=sys.stderr)

    write_fasta([(n, s) for n, _, s in read_fasta(args.fasta) if n in chosen], args.out_fasta)
    with open(args.out_groups, "w") as out:
        out.write("id\tgroup\n")
        for sid, g in chosen.items():
            out.write(f"{sid}\t{g}\n")


if __name__ == "__main__":
    main()
