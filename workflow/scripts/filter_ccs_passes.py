#!/usr/bin/env python3
# =============================================================================
# filter_ccs_passes.py: remove CCS reads with fewer than N full passes.
#
# The number of passes is stored as the "np" tag in PacBio CCS BAM files; it
# is not in FASTQ. With --bam, read names with np >= --min-passes are
# collected from the BAM and only those reads are kept from the FASTQ.
# Without a BAM (--bam ''), reads are copied unchanged and a warning is
# printed: then report that the pass filter could not be applied.
# =============================================================================

import argparse
import gzip
import shutil
import sys

from yfutils import read_fastq, write_fastq_record


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--bam", default="")
    ap.add_argument("--min-passes", type=int, default=10)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if not args.bam:
        print("WARNING: no CCS BAM given; pass filter skipped, reads copied unchanged", file=sys.stderr)
        if args.fastq.endswith(".gz"):
            shutil.copyfile(args.fastq, args.out)
        else:
            # Compress so that every later step can rely on .gz input.
            with open(args.fastq, "rb") as fin, gzip.open(args.out, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        return

    import pysam  # only needed when a BAM is given

    keep = set()
    total = 0
    with pysam.AlignmentFile(args.bam, "rb", check_sq=False) as bam:
        for read in bam:
            total += 1
            if read.has_tag("np") and read.get_tag("np") >= args.min_passes:
                keep.add(read.query_name)
    print(f"{len(keep)} of {total} CCS reads have >= {args.min_passes} passes", file=sys.stderr)

    n = 0
    with gzip.open(args.out, "wt") as out:
        for name, _, seq, qual in read_fastq(args.fastq):
            if name in keep:
                write_fastq_record(out, name, seq, qual)
                n += 1
    print(f"{n} reads written", file=sys.stderr)


if __name__ == "__main__":
    main()
