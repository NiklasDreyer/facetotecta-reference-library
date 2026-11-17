#!/usr/bin/env python3
# =============================================================================
# concat_alignments.py: concatenate per-locus alignments by specimen.
#
# - the specimen key is the sequence ID (e.g. SAN-Y123), the same across loci
# - specimens missing a locus get "?" (missing data) for that locus
# - writes a NEXUS partition file for IQ-TREE: one partition per non-coding
#   locus and one per codon position for coding loci (COX1, H3), as in
#   manuscript_v1. Coding alignments must already be in frame (start at codon
#   position 1); sanger_frame + Gblocks -t=c ensure this.
# - writes a specimen x locus coverage table and metrics (number of specimens,
#   aligned length, mean / min / max non-missing length per specimen)
#
# Usage:
#   concat_alignments.py --alignments COX1.gb.fasta 18S.gb.fasta ... \
#       --loci COX1 18S ... --coding COX1 H3 --out-fasta concat.fasta \
#       --out-nexus partitions.nex --out-coverage cov.tsv --metrics m.tsv
# =============================================================================

import argparse
import statistics

from yfutils import read_fasta, write_fasta, write_metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignments", nargs="+", required=True)
    ap.add_argument("--loci", nargs="+", required=True, help="locus names, same order as --alignments")
    ap.add_argument("--coding", nargs="*", default=[])
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-nexus", required=True)
    ap.add_argument("--out-coverage", required=True)
    ap.add_argument("--metrics")
    args = ap.parse_args()
    assert len(args.alignments) == len(args.loci), "one locus name per alignment"

    # Read every alignment: locus -> {specimen: aligned sequence}.
    data, lengths = {}, {}
    for locus, path in zip(args.loci, args.alignments):
        seqs = {n: s.upper() for n, _, s in read_fasta(path)}
        if not seqs:
            print(f"warning: {path} is empty; locus skipped")
            continue
        data[locus] = seqs
        lengths[locus] = len(next(iter(seqs.values())))

    loci = [l for l in args.loci if l in data]
    specimens = sorted(set().union(*[set(d) for d in data.values()]))

    # Build the concatenated sequences and the coverage table.
    concat = []
    cov_rows = []
    for sp in specimens:
        parts, present = [], []
        for l in loci:
            s = data[l].get(sp)
            if s is None:
                parts.append("?" * lengths[l])
                present.append(0)
            else:
                parts.append(s)
                present.append(len(s.replace("-", "").replace("?", "")))
        concat.append((sp, "".join(parts)))
        cov_rows.append((sp, present))
    write_fasta(concat, args.out_fasta)

    # Coverage table: number of non-gap bases per locus and in total.
    with open(args.out_coverage, "w") as fh:
        fh.write("specimen\t" + "\t".join(loci) + "\ttotal_bp\tn_loci\n")
        for sp, present in cov_rows:
            fh.write(f"{sp}\t" + "\t".join(map(str, present)) +
                     f"\t{sum(present)}\t{sum(1 for p in present if p)}\n")

    # NEXUS partition file with codon positions for coding loci.
    with open(args.out_nexus, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        start = 1
        for l in loci:
            end = start + lengths[l] - 1
            if l in args.coding:
                for pos in range(3):
                    fh.write(f"    charset {l}_pos{pos + 1} = {start + pos}-{end}\\3;\n")
            else:
                fh.write(f"    charset {l} = {start}-{end};\n")
            start = end + 1
        fh.write("end;\n")

    totals = [sum(p) for _, p in cov_rows]
    write_metrics(args.metrics, {
        "concat_specimens": len(specimens),
        "concat_aligned_length": sum(lengths.values()),
        "concat_mean_bp_per_specimen": round(statistics.mean(totals), 1) if totals else 0,
        "concat_min_bp_per_specimen": min(totals) if totals else 0,
        "concat_max_bp_per_specimen": max(totals) if totals else 0,
        **{f"concat_aligned_length_{l}": lengths[l] for l in loci},
    })


if __name__ == "__main__":
    main()
