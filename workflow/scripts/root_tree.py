#!/usr/bin/env python3
# =============================================================================
# root_tree.py: midpoint-root a Newick tree (PTP and mPTP need a rooted tree).
#
# Support values written by IQ-TREE (e.g. "85.2/97") are kept as node names.
# An empty tree "();" (dataset with fewer than 4 sequences) is copied as is.
# =============================================================================

import argparse
import shutil

from ete3 import Tree


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    text = open(args.input).read().strip()
    if text in ("();", ""):
        shutil.copyfile(args.input, args.output)
        return
    t = Tree(text, format=1)
    # Root at the midpoint of the longest tip-to-tip path.
    t.set_outgroup(t.get_midpoint_outgroup())
    # format=5: branch lengths for all nodes, no internal names. mPTP does
    # not need support values and may not parse them.
    t.write(outfile=args.output, format=5)


if __name__ == "__main__":
    main()
