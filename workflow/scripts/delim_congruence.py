#!/usr/bin/env python3
# =============================================================================
# delim_congruence.py: NEW agreement between delimitations, and MOTU grades.
#
# 1. Within each dataset, for every pair of methods (VSEARCH, PTP, mPTP, ASAP):
#      - adjusted Rand index (ARI; 1 = identical partitions)
#      - match ratio (Ahrens et al. 2016) = 2 x (MOTUs with identical
#        membership) / (MOTUs in partition 1 + MOTUs in partition 2)
#    computed on the sequences present in both partitions.
# 2. Between the primary (e.g. sanger_COX1) and nuclear (e.g. sanger_18S)
#    datasets, per method, on specimens sequenced for both loci.
# 3. Grades for the MOTUs of the primary dataset (VSEARCH partition):
#      A  identical in every other method AND not split by the nuclear locus
#      B  identical in at least one other method
#      C  found by one method only
#    "Not split by the nuclear locus": all specimens of the MOTU that have
#    nuclear data fall in a single nuclear MOTU (nuclear loci may lump
#    species, so lumping is allowed). MOTUs without nuclear data can reach
#    grade B at most.
# =============================================================================

import argparse
import itertools
import os
from collections import defaultdict

from sklearn.metrics import adjusted_rand_score

from yfutils import read_membership, write_metrics


def groups(members, keep=None):
    """MOTU -> frozenset of sequences (optionally restricted to 'keep')."""
    g = defaultdict(set)
    for s, m in members.items():
        if keep is None or s in keep:
            g[m].add(s)
    return {m: frozenset(v) for m, v in g.items()}


def compare(a, b):
    """(n shared, ARI, match ratio) between two membership dicts."""
    shared = sorted(set(a) & set(b))
    if len(shared) < 2:
        return len(shared), None, None
    ari = adjusted_rand_score([a[s] for s in shared], [b[s] for s in shared])
    ga, gb = set(groups(a, set(shared)).values()), set(groups(b, set(shared)).values())
    match = 2 * len(ga & gb) / (len(ga) + len(gb))
    return len(shared), round(ari, 4), round(match, 4)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--delim-dir", required=True)
    ap.add_argument("--datasets", nargs="+", required=True)
    ap.add_argument("--methods", nargs="+", required=True)
    ap.add_argument("--primary", required=True)
    ap.add_argument("--nuclear", required=True)
    ap.add_argument("--out-pairs", required=True)
    ap.add_argument("--out-grades", required=True)
    ap.add_argument("--metrics")
    args = ap.parse_args()

    def load(ds, method):
        p = os.path.join(args.delim_dir, ds, "membership", f"{method}.tsv")
        return read_membership(p) if os.path.exists(p) else {}

    metrics = {}
    with open(args.out_pairs, "w") as out:
        out.write("comparison\tdataset_a\tmethod_a\tdataset_b\tmethod_b\tshared\tARI\tmatch_ratio\n")
        # 1. Between methods within each dataset.
        for ds in args.datasets:
            for m1, m2 in itertools.combinations(args.methods, 2):
                n, ari, mr = compare(load(ds, m1), load(ds, m2))
                out.write(f"methods\t{ds}\t{m1}\t{ds}\t{m2}\t{n}\t{ari}\t{mr}\n")
                if ds == args.primary:
                    metrics[f"congruence_{ds}_{m1}_vs_{m2}_ARI"] = ari
                    metrics[f"congruence_{ds}_{m1}_vs_{m2}_match_ratio"] = mr
        # 2. Between loci (same method), on shared specimens.
        for m in args.methods:
            n, ari, mr = compare(load(args.primary, m), load(args.nuclear, m))
            out.write(f"loci\t{args.primary}\t{m}\t{args.nuclear}\t{m}\t{n}\t{ari}\t{mr}\n")
            metrics[f"congruence_{args.primary}_vs_{args.nuclear}_{m}_shared"] = n
            metrics[f"congruence_{args.primary}_vs_{args.nuclear}_{m}_ARI"] = ari
            metrics[f"congruence_{args.primary}_vs_{args.nuclear}_{m}_match_ratio"] = mr

    # 3. Grades for the primary dataset's VSEARCH MOTUs.
    base = load(args.primary, "vsearch")
    others = {m: load(args.primary, m) for m in args.methods if m != "vsearch"}
    others = {m: v for m, v in others.items() if v}
    nuclear = load(args.nuclear, "vsearch")
    grade_counts = defaultdict(int)
    with open(args.out_grades, "w") as out:
        out.write("motu\tn_sequences\tmethods_identical\tmethods_total\tnuclear_status\tgrade\tmembers\n")
        for motu, seqs in sorted(groups(base).items(), key=lambda x: -len(x[1])):
            identical = 0
            for m, mem in others.items():
                # Compare on sequences present in both partitions: the MOTU is
                # "identical" if the other method has a MOTU with exactly the
                # same members.
                shared = set(base) & set(mem)
                present = frozenset(seqs & shared)
                if present and present in set(groups(mem, shared).values()):
                    identical += 1
            # Nuclear check: specimens with nuclear data in a single nuclear MOTU.
            nuc = {nuclear[s] for s in seqs if s in nuclear}
            nuc_status = "no_data" if not nuc else ("congruent" if len(nuc) == 1 else "split")
            if others and identical == len(others) and nuc_status == "congruent":
                grade = "A"
            elif identical >= 1:
                grade = "B"
            else:
                grade = "C"
            grade_counts[grade] += 1
            out.write(f"{motu}\t{len(seqs)}\t{identical}\t{len(others)}\t{nuc_status}\t{grade}\t"
                      f"{','.join(sorted(seqs))}\n")
    for g in "ABC":
        metrics[f"grades_{args.primary}_{g}"] = grade_counts[g]
    write_metrics(args.metrics, metrics)


if __name__ == "__main__":
    main()
