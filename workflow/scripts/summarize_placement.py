#!/usr/bin/env python3
# =============================================================================
# summarize_placement.py: how many NON-TARGET sequences were placed inside
# Facetotecta with LWR >= --min-lwr (false placements), overall and per group.
# =============================================================================

import argparse

import pandas as pd

from yfutils import write_metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--min-lwr", type=float, default=0.9)
    ap.add_argument("--label", default="placement")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    a = pd.read_csv(args.assign, sep="\t", dtype={"query": str})
    g = pd.read_csv(args.groups, sep="\t", dtype=str).set_index("id")["group"]
    a["group"] = a["query"].map(g).fillna("Other")
    a["false_facet"] = a["taxopath"].astype(str).str.contains("Facetotecta") & (a["LWR"] >= args.min_lwr)

    m = {f"{args.label}_queries": len(a),
         f"{args.label}_placed_in_facetotecta": int(a["false_facet"].sum()),
         f"{args.label}_false_placement_rate": round(a["false_facet"].mean(), 4) if len(a) else "NA"}
    for grp, sub in a.groupby("group"):
        m[f"{args.label}_{grp}_placed_in_facetotecta"] = f"{int(sub.false_facet.sum())}/{len(sub)}"
    write_metrics(args.output, m)


if __name__ == "__main__":
    main()
