#!/usr/bin/env python3
# =============================================================================
# make_subsets.py: split the authentic sequences of one locus by platform.
#
# After pre-clustering, each authentic centroid stands for one or more
# original sequences of a single platform. This script lists, for every
# centroid, the platforms/samples it represents:
#   SAN       Sanger
#   PB        PacBio (pooled Y1+Y2 RNA)
#   ILL_<s>   Illumina sample <s> (from the DADA2 ASV table: count > 0)
# and writes one FASTA per expected tag (subset_<tag>.fasta), even if empty,
# plus a table centroid <tab> comma-separated tags (used for the UpSet plot).
# =============================================================================

import argparse
import csv
import os
from collections import defaultdict

from yfutils import read_fasta, write_fasta


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="authentic centroids")
    ap.add_argument("--members", required=True, help="precluster_members.tsv")
    ap.add_argument("--asv-table", help="DADA2 ASV x sample table (Illumina loci)")
    ap.add_argument("--expected", nargs="*", default=[], help="tags to write subsets for")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--platforms-out", required=True)
    args = ap.parse_args()

    authentic = {n: s for n, _, s in read_fasta(args.fasta)}

    # Samples in which each Illumina ASV occurs.
    asv_samples = defaultdict(set)
    if args.asv_table:
        with open(args.asv_table) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            samples = header[1:]
            for row in reader:
                for s, c in zip(samples, row[1:]):
                    if float(c) > 0:
                        asv_samples[row[0]].add(s)

    # Tags per centroid from its members.
    tags = defaultdict(set)
    with open(args.members) as fh:
        next(fh)  # header
        for line in fh:
            centroid, member, plat = line.rstrip("\n").split("\t")
            if centroid not in authentic:
                continue
            if plat == "ILL":
                for s in asv_samples.get(member, []):
                    tags[centroid].add(f"ILL_{s}")
            else:
                tags[centroid].add(plat)

    # Table of tags per authentic sequence.
    with open(args.platforms_out, "w") as out:
        out.write("sequence\tplatforms\n")
        for n in authentic:
            out.write(f"{n}\t{','.join(sorted(tags.get(n, [])))}\n")

    # One FASTA per tag.
    all_tags = set(args.expected) | {t for ts in tags.values() for t in ts}
    for t in sorted(all_tags):
        if t == "SAN":
            continue
        recs = [(n, s) for n, s in authentic.items() if t in tags.get(n, set())]
        write_fasta(recs, os.path.join(args.outdir, f"subset_{t}.fasta"))
        print(f"subset {t}: {len(recs)} sequences")


if __name__ == "__main__":
    main()
