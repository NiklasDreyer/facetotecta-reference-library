#!/usr/bin/env python3
# =============================================================================
# parse_uc.py: VSEARCH .uc cluster file -> membership table.
#
# .uc records: "S" = centroid (column 9 = its label), "H" = member (column 9 =
# member label, column 10 = its centroid). Each cluster is named after its
# centroid. Output: sequence <tab> MOTU.
# =============================================================================

import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uc", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows = []
    with open(args.uc) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                rows.append((f[8].split()[0], f[8].split()[0]))
            elif f[0] == "H":
                rows.append((f[8].split()[0], f[9].split()[0]))

    with open(args.output, "w") as out:
        out.write("sequence\tmotu\n")
        for seq, motu in rows:
            out.write(f"{seq}\t{motu}\n")


if __name__ == "__main__":
    main()
