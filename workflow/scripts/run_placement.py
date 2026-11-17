#!/usr/bin/env python3
# =============================================================================
# run_placement.py: phylogenetic placement of query sequences.
#
#   1. MAFFT --addfragments --keeplength: add queries to the fixed reference
#      alignment (columns of the reference are not changed)
#   2. split the result into reference and query alignments
#   3. EPA-ng: place each query on the reference tree (likelihood weight
#      ratio, LWR, for each candidate branch)
#   4. gappa examine assign: turn placements into a taxonomic label using the
#      reference taxonomy (Facetotecta / Outgroup); the best label per query
#      and its LWR are written to the output table
# =============================================================================

import argparse
import csv
import os
import subprocess
import sys

from yfutils import read_fasta, write_fasta


def run(cmd):
    """Run a command and stop on failure."""
    print("$ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--query", required=True)
    ap.add_argument("--ref-aln", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--tax", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)
    jplace = os.path.join(args.workdir, "epa_result.jplace")

    queries = [(n, s.replace("-", "")) for n, _, s in read_fasta(args.query)]
    if not queries:
        # Nothing to place: write empty results.
        open(jplace, "w").close()
        open(args.output, "w").write("query\ttaxopath\tLWR\n")
        return

    # 1. Add queries to the reference alignment.
    q_clean = os.path.join(args.workdir, "query.fasta")
    write_fasta(queries, q_clean)
    combined = os.path.join(args.workdir, "combined.aln.fasta")
    with open(combined, "w") as fh:
        subprocess.run(["mafft", "--addfragments", q_clean, "--keeplength",
                        "--thread", str(args.threads), args.ref_aln], check=True, stdout=fh)

    # 2. Keep only the query rows.
    ref_ids = {n for n, _, _ in read_fasta(args.ref_aln)}
    q_aln = os.path.join(args.workdir, "query.aln.fasta")
    write_fasta([(n, s) for n, _, s in read_fasta(combined) if n not in ref_ids], q_aln)

    # 3. EPA-ng placement.
    run(["epa-ng", "--ref-msa", args.ref_aln, "--tree", args.tree, "--query", q_aln,
         "--model", args.model, "--threads", str(args.threads),
         "--out-dir", args.workdir, "--redo"])

    # 4. gappa assignment (best hit per query).
    gdir = os.path.join(args.workdir, "gappa")
    os.makedirs(gdir, exist_ok=True)
    run(["gappa", "examine", "assign", "--jplace-path", jplace, "--taxon-file", args.tax,
         "--out-dir", gdir, "--per-query-results", "--best-hit", "--allow-file-overwriting"])

    # gappa writes per_query.tsv: name, LWR, fract, aLWR, afract, taxopath.
    per_query = os.path.join(gdir, "per_query.tsv")
    best = {}
    with open(per_query) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            name = r.get("name") or r.get("query") or list(r.values())[0]
            lwr = float(r.get("LWR", 0) or 0)
            if name not in best or lwr > best[name][1]:
                best[name] = (r.get("taxopath", ""), lwr)
    with open(args.output, "w") as out:
        out.write("query\ttaxopath\tLWR\n")
        for n, _ in queries:
            tp, lwr = best.get(n, ("unplaced", 0.0))
            out.write(f"{n}\t{tp}\t{lwr}\n")


if __name__ == "__main__":
    main()
