#!/usr/bin/env python3
# =============================================================================
# trim_to_reference.py: extract the part of each mapped read that aligns to
# the reference (manuscript_v1: "we trimmed the 5' and 3' ends to match the
# complete COX1 reference (1,527 bp)").
#
# - only primary alignments (secondary, supplementary, unmapped skipped)
# - reads with mapping quality below --min-mapq are counted as ambiguous and
#   removed (manuscript_v1 removed ambiguous reads); BBMap gives low MAPQ to
#   reads that map equally well to several places
# - soft-clipped ends (read parts beyond the reference) are removed
# - sequences are written in reference orientation (BAM stores reverse-strand
#   reads reverse-complemented), with "refstart=<0-based start>" in the header
#   so check_translation.py can find the reading frame
# =============================================================================

import argparse
import sys

import pysam

from yfutils import write_metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--metrics")
    args = ap.parse_args()

    n = {"primary_mapped": 0, "ambiguous_removed": 0, "short_removed": 0, "written": 0}
    with pysam.AlignmentFile(args.bam, "rb") as bam, open(args.out, "w") as out:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            n["primary_mapped"] += 1
            if r.mapping_quality < args.min_mapq:
                n["ambiguous_removed"] += 1
                continue
            # The aligned part of the read, without soft clips.
            seq = r.query_alignment_sequence
            if seq is None or len(seq) < args.min_len:
                n["short_removed"] += 1
                continue
            out.write(f">{r.query_name} refstart={r.reference_start} "
                      f"refend={r.reference_end}\n{seq}\n")
            n["written"] += 1

    print(n, file=sys.stderr)
    write_metrics(args.metrics, {f"pacbio_COX1_map_{k}": v for k, v in n.items()})


if __name__ == "__main__":
    main()
