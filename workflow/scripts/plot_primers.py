#!/usr/bin/env python3
# =============================================================================
# plot_primers.py: Fig. 6b. For each primer pair, the share of facetotectan
# sequences with a mismatch at each primer position (forward and reverse
# primer, numbered 5'->3'), from primer_scan.py tables. Only sequences whose
# primer sites were covered are counted.
# =============================================================================

import argparse
import os
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def position_profile(series):
    """Count mismatches per primer position from comma-separated lists."""
    c = Counter()
    for v in series.dropna():
        for p in str(v).split(","):
            if p:
                c[int(p)] += 1
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tables", nargs="+", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    fig, axes = plt.subplots(len(args.tables), 2, figsize=(9, 2.2 * len(args.tables)), squeeze=False)
    for row, path in enumerate(args.tables):
        name = os.path.basename(os.path.dirname(path))
        t = pd.read_csv(path, sep="\t", dtype=str)
        t = t[t["status"] != "site_not_covered"]
        n = len(t)
        for col, (key, label) in enumerate([("fwd_mm_positions", "forward"), ("rev_mm_positions", "reverse")]):
            prof = position_profile(t[key]) if n else Counter()
            ax = axes[row][col]
            if prof:
                xs = range(1, max(prof) + 1)
                ax.bar(list(xs), [prof.get(x, 0) / n for x in xs], color="#4C72B0")
            ax.set_ylim(0, 1)
            ax.set_title(f"{name}: {label} primer (n = {n})", fontsize=8)
            ax.set_xlabel("Primer position (5'->3')", fontsize=7)
            ax.set_ylabel("Share mismatched", fontsize=7)
            ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(args.output)


if __name__ == "__main__":
    main()
