#!/usr/bin/env python3
# =============================================================================
# rename_reads.py: give reads short IDs and write FASTQ + FASTA copies.
#
# PacBio read names (m54xxx/12345/ccs) contain "/", which IQ-TREE and other
# tree programs rewrite. Reads are renamed <prefix><6-digit number>
# (e.g. PB-000001) and a table maps new names to original names.
# =============================================================================

import argparse
import gzip

from yfutils import read_fastq, write_fastq_record


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--prefix", default="PB-")
    ap.add_argument("--out-fastq", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--table", required=True)
    args = ap.parse_args()

    with gzip.open(args.out_fastq, "wt") as fq, open(args.out_fasta, "w") as fa, \
            open(args.table, "w") as tab:
        tab.write("new_id\toriginal_id\n")
        for i, (name, _, seq, qual) in enumerate(read_fastq(args.input), start=1):
            new = f"{args.prefix}{i:06d}"
            write_fastq_record(fq, new, seq, qual)
            fa.write(f">{new}\n{seq}\n")
            tab.write(f"{new}\t{name}\n")


if __name__ == "__main__":
    main()
