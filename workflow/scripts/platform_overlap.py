#!/usr/bin/env python3
# =============================================================================
# platform_overlap.py: which platforms found each MOTU (Fig. 4b).
#
# Inputs: MOTU membership of the combined HTS dataset (e.g. hts_COX1_all,
# VSEARCH) and the table of platforms per authentic sequence (from
# make_subsets.py). For each MOTU the union of platforms of its sequences is
# taken (SAN, PB, ILL_Y1, ILL_Y2). Outputs a MOTU table, an UpSet plot and
# metrics: total MOTUs, MOTUs per platform, MOTUs found by Sanger plus at least
# one HTS platform, MOTUs found only by HTS platforms.
# =============================================================================

import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet, from_memberships

from yfutils import read_membership, write_metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--members", required=True)
    ap.add_argument("--platforms", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-pdf", required=True)
    ap.add_argument("--metrics")
    ap.add_argument("--label", default="locus")
    args = ap.parse_args()

    members = read_membership(args.members)
    plat = pd.read_csv(args.platforms, sep="\t", dtype=str).fillna("")
    plat = {r.sequence: set(filter(None, r.platforms.split(","))) for r in plat.itertuples()}

    motu_plats = defaultdict(set)
    for seq, motu in members.items():
        motu_plats[motu] |= plat.get(seq, set())

    rows = [(m, ",".join(sorted(p))) for m, p in motu_plats.items()]
    pd.DataFrame(rows, columns=["motu", "platforms"]).to_csv(args.out_table, sep="\t", index=False)

    all_tags = sorted({t for p in motu_plats.values() for t in p})
    m = {f"overlap_{args.label}_motus": len(motu_plats)}
    for t in all_tags:
        m[f"overlap_{args.label}_motus_with_{t}"] = sum(1 for p in motu_plats.values() if t in p)
    m[f"overlap_{args.label}_motus_sanger_and_hts"] = sum(1 for p in motu_plats.values() if "SAN" in p and len(p) > 1)
    m[f"overlap_{args.label}_motus_hts_only"] = sum(1 for p in motu_plats.values() if p and "SAN" not in p)
    m[f"overlap_{args.label}_motus_sanger_only"] = sum(1 for p in motu_plats.values() if p == {"SAN"})
    write_metrics(args.metrics, m)

    # UpSet plot: one bar per combination of platforms.
    memberships = [sorted(p) for p in motu_plats.values() if p]
    if memberships:
        data = from_memberships(memberships)
        UpSet(data, subset_size="count", show_counts=True, sort_by="cardinality").plot()
        plt.suptitle(f"{args.label}: MOTUs by platform")
    plt.savefig(args.out_pdf)


if __name__ == "__main__":
    main()
