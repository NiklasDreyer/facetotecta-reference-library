#!/usr/bin/env python3
# =============================================================================
# count_motus.py: number of MOTUs and singleton MOTUs per method for one dataset.
# Membership files are named <method>.tsv (sequence <tab> MOTU).
# =============================================================================

import argparse
import os
from collections import Counter

from yfutils import read_membership, write_metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--memberships", nargs="+", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    m = {}
    for path in args.memberships:
        method = os.path.basename(path).replace(".tsv", "")
        members = read_membership(path)
        sizes = Counter(members.values())
        m[f"delim_{args.dataset}_{method}_sequences"] = len(members)
        m[f"delim_{args.dataset}_{method}_motus"] = len(sizes)
        m[f"delim_{args.dataset}_{method}_singletons"] = sum(1 for s in sizes.values() if s == 1)
    write_metrics(args.output, m)


if __name__ == "__main__":
    main()
