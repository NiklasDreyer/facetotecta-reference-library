#!/usr/bin/env python3
# =============================================================================
# summarize_numbers.py: merge every *.metrics.tsv file under the results
# folder into one table (key, value, source file), sorted by key.
#
# This is the single source for every number in the manuscript: counts of
# specimens, sequences and reads at each step, MOTU counts per dataset and
# method, congruence, richness estimates, primer coverage, benchmark rates
# and (after the mining workflow) public-data detections.
# =============================================================================

import argparse
import glob
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows = []
    for path in sorted(glob.glob(os.path.join(args.results, "**", "*.metrics.tsv"), recursive=True)):
        rel = os.path.relpath(path, args.results)
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0]:
                    rows.append((parts[0], parts[1], rel))

    # samtools flagstat files (PacBio COX1 mapping: mapPacBio.sh vs minimap2):
    # take the number of primary mapped reads, or "mapped" in older versions.
    for path in sorted(glob.glob(os.path.join(args.results, "**", "*.flagstat.txt"), recursive=True)):
        rel = os.path.relpath(path, args.results)
        name = os.path.basename(path).replace(".flagstat.txt", "")
        mapped = None
        for line in open(path):
            if "primary mapped (" in line or (mapped is None and " mapped (" in line):
                mapped = line.split()[0]
        if mapped is not None:
            rows.append((f"pacbio_COX1_{name}_mapped_reads", mapped, rel))

    rows.sort(key=lambda r: r[0])
    with open(args.output, "w") as out:
        out.write("key\tvalue\tsource\n")
        for k, v, src in rows:
            out.write(f"{k}\t{v}\t{src}\n")
    print(f"{len(rows)} numbers collected into {args.output}")


if __name__ == "__main__":
    main()
