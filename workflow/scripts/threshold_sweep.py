#!/usr/bin/env python3
# =============================================================================
# threshold_sweep.py: number of VSEARCH clusters (MOTUs) at a range of
# identity thresholds (NEW, Fig. S8). Shows how strongly MOTU counts depend
# on the threshold, and where the chosen threshold sits.
# =============================================================================

import argparse
import os
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")  # no display on cluster nodes
import matplotlib.pyplot as plt


def count_clusters(fasta, ident, threads):
    """Number of clusters from vsearch --cluster_fast at one identity."""
    with tempfile.NamedTemporaryFile(suffix=".uc", delete=False) as tmp:
        uc = tmp.name
    subprocess.run(["vsearch", "--cluster_fast", fasta, "--id", f"{ident:.4f}",
                    "--uc", uc, "--threads", str(threads), "--quiet"], check=True)
    n = sum(1 for line in open(uc) if line.startswith("C\t"))
    os.remove(uc)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--start", type=float, default=0.90)
    ap.add_argument("--end", type=float, default=0.99)
    ap.add_argument("--step", type=float, default=0.005)
    ap.add_argument("--chosen", type=float)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-pdf", required=True)
    args = ap.parse_args()

    ids, counts = [], []
    if os.path.getsize(args.input) > 0:
        x = args.start
        while x <= args.end + 1e-9:
            ids.append(round(x, 4))
            counts.append(count_clusters(args.input, x, args.threads))
            x += args.step

    with open(args.out_tsv, "w") as out:
        out.write("identity\tmotus\n")
        for i, c in zip(ids, counts):
            out.write(f"{i}\t{c}\n")

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(ids, counts, marker="o", ms=3)
    if args.chosen:
        ax.axvline(args.chosen, ls="--", color="grey", label=f"chosen ({args.chosen})")
        ax.legend(frameon=False)
    ax.set_xlabel("VSEARCH identity threshold")
    ax.set_ylabel("Number of MOTUs")
    fig.tight_layout()
    fig.savefig(args.out_pdf)


if __name__ == "__main__":
    main()
