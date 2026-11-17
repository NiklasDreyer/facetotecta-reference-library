#!/usr/bin/env python3
# =============================================================================
# parse_mptp.py: turn mPTP/PTP output into a membership table.
#
# mPTP writes the delimitation as blocks:
#     Species 1:
#     seqA
#     seqB
#     (blank line)
#     Species 2:
#     ...
# Output: sequence <tab> MOTU (MOTU names "mptp_1", "mptp_2", ...).
# =============================================================================

import argparse
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--motu-prefix", default="motu")
    args = ap.parse_args()

    current = None
    rows = []
    with open(args.input) as fh:
        for line in fh:
            line = line.strip()
            m = re.match(r"Species\s+(\d+):", line)
            if m:
                current = m.group(1)
                continue
            if not line:
                current = None  # a blank line ends a species block
                continue
            # Sequence names appear only inside "Species" blocks.
            if current is not None and not line.startswith(("Number", "Null", "Score", "LRT", "Minimum")):
                rows.append((line, f"{args.motu_prefix}_{current}"))

    with open(args.output, "w") as out:
        out.write("sequence\tmotu\n")
        for seq, motu in rows:
            out.write(f"{seq}\t{motu}\n")


if __name__ == "__main__":
    main()
