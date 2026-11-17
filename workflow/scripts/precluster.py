#!/usr/bin/env python3
# =============================================================================
# precluster.py: pre-cluster each platform separately before the
# phylogenetic authentication filter (manuscript_v1: VSEARCH --cluster_fast
# --id 0.97).
#
# Each platform (SAN, PB, ILL) is clustered on its own so that every centroid
# keeps one platform tag. With --id 1.0 the step is plain dereplication
# (recommended; see config). Outputs:
#   - all centroids in one FASTA
#   - a membership table: centroid <tab> member <tab> platform
#   - metrics: sequences in / centroids out per platform
#
# Usage: precluster.py --inputs SAN=a.fasta PB=b.fasta ILL=c.fasta --id 0.97 ...
# =============================================================================

import argparse
import os
import subprocess
import sys
import tempfile

from yfutils import read_fasta, write_metrics


def run_vsearch(fasta, ident, threads, workdir, tag):
    """Cluster one FASTA; return (centroids FASTA path, uc path)."""
    cent = os.path.join(workdir, f"{tag}.centroids.fasta")
    uc = os.path.join(workdir, f"{tag}.uc")
    if ident >= 1.0:
        # Dereplication: identical sequences only (full length).
        cmd = ["vsearch", "--derep_fulllength", fasta, "--output", cent,
               "--uc", uc, "--minseqlength", "1", "--quiet"]
    else:
        cmd = ["vsearch", "--cluster_fast", fasta, "--id", str(ident),
               "--centroids", cent, "--uc", uc, "--threads", str(threads), "--quiet"]
    subprocess.run(cmd, check=True)
    return cent, uc


def parse_uc(uc):
    """Yield (centroid, member) pairs from a VSEARCH .uc file."""
    with open(uc) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":           # centroid record
                yield f[8], f[8]
            elif f[0] == "H":         # hit: member -> centroid
                yield f[9], f[8]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True, help="TAG=path pairs")
    ap.add_argument("--id", type=float, required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-members", required=True)
    ap.add_argument("--metrics")
    ap.add_argument("--label", default="locus")
    args = ap.parse_args()

    metrics = {}
    with tempfile.TemporaryDirectory() as tmp, \
            open(args.out_fasta, "w") as fa, open(args.out_members, "w") as mem:
        mem.write("centroid\tmember\tplatform\n")
        for item in args.inputs:
            tag, path = item.split("=", 1)
            n_in = sum(1 for _ in read_fasta(path))
            if n_in == 0:
                print(f"{tag}: no sequences", file=sys.stderr)
                continue
            cent, uc = run_vsearch(path, args.id, args.threads, tmp, tag)
            n_out = 0
            for name, _, seq in read_fasta(cent):
                fa.write(f">{name}\n{seq}\n")
                n_out += 1
            for centroid, member in parse_uc(uc):
                mem.write(f"{centroid}\t{member}\t{tag}\n")
            metrics[f"precluster_{args.label}_{tag}_input"] = n_in
            metrics[f"precluster_{args.label}_{tag}_centroids"] = n_out
            print(f"{tag}: {n_in} sequences -> {n_out} centroids", file=sys.stderr)
    write_metrics(args.metrics, metrics)


if __name__ == "__main__":
    main()
