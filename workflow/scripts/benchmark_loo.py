#!/usr/bin/env python3
# =============================================================================
# benchmark_loo.py: NEW leave-one-MOTU-out test of taxonomic assignment.
#
# Question: if a facetotectan lineage is MISSING from the reference library,
# is a read from it still recognised as Facetotecta, and how often are reads
# from other crustaceans wrongly called Facetotecta?
#
# For every facetotectan MOTU m (from the Sanger VSEARCH delimitation):
#   reference = all facetotectan sequences EXCEPT those of m
#             + non-target sequences not in this round's test fold
#   queries   = sequences of m (positives)
#             + the test fold of non-target sequences (negatives; each fold
#               is queried once, in the first round that uses it)
#   assigned with
#     BLAST  : top hit (highest bitscore); "Facetotecta" if the top hit is a
#              facetotectan sequence; its percent identity is recorded
#     SINTAX : VSEARCH --sintax, label at the order level ("o:Facetotecta"
#              or "o:<non-target group>"), with its bootstrap confidence
# Then, for every identity threshold (BLAST) or confidence cut-off (SINTAX):
#   true-positive rate  = positives called Facetotecta / all positives
#   false-positive rate = negatives called Facetotecta / all negatives
# The chosen threshold is the lowest one whose false-positive rate is at or
# below --target-fpr. manuscript_v1's cut-offs (92% and 90%) are also scored.
# =============================================================================

import argparse
import os
import random
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from yfutils import read_fasta, read_membership, write_fasta, write_metrics


def blast_best(query, ref, workdir, threads):
    """Best BLAST hit per query: {query: (subject, pident, qcov)}."""
    db = os.path.join(workdir, "refdb")
    subprocess.run(["makeblastdb", "-in", ref, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL)
    out = os.path.join(workdir, "blast.tsv")
    subprocess.run(["blastn", "-query", query, "-db", db, "-max_target_seqs", "10",
                    "-num_threads", str(threads), "-evalue", "1e-5",
                    # stitle = full subject header, which carries the ";tax=" label
                    # (sseqid formatting differs between BLAST versions).
                    "-outfmt", "6 qseqid stitle pident length evalue bitscore qlen",
                    "-out", out], check=True)
    best = {}
    for line in open(out):
        q, s, pid, length, ev, bits, qlen = line.rstrip("\n").split("\t")
        bits = float(bits)
        if q not in best or bits > best[q][3]:
            best[q] = (s, float(pid), int(length) / int(qlen), bits)
    return {q: v[:3] for q, v in best.items()}


def sintax(query, ref, workdir, threads, seed):
    """SINTAX order-level label and confidence per query: {query: (label, conf)}."""
    out = os.path.join(workdir, "sintax.tsv")
    subprocess.run(["vsearch", "--sintax", query, "--db", ref, "--tabbedout", out,
                    "--sintax_cutoff", "0.0", "--randseed", str(seed),
                    "--threads", str(threads), "--quiet"], check=True)
    res = {}
    for line in open(out):
        f = line.rstrip("\n").split("\t")
        label, conf = "unassigned", 0.0
        # Column 2 looks like "k:Arthropoda(1.00),o:Facetotecta(0.87)".
        for part in f[1].split(","):
            if part.startswith("o:"):
                label = part[2:part.index("(")]
                conf = float(part[part.index("(") + 1:-1])
        res[f[0]] = (label, conf)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--facet", required=True)
    ap.add_argument("--nontarget", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--members", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--folds", type=int, default=10)
    ap.add_argument("--tmin", type=float, default=80)
    ap.add_argument("--tmax", type=float, default=100)
    ap.add_argument("--target-fpr", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out-hits", required=True)
    ap.add_argument("--out-roc", required=True)
    ap.add_argument("--out-pdf", required=True)
    ap.add_argument("--out-thresholds", required=True)
    ap.add_argument("--metrics")
    ap.add_argument("--label", default="bench")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    facet = {n: s.replace("-", "") for n, _, s in read_fasta(args.facet)}
    nontarget = {n: s.replace("-", "") for n, _, s in read_fasta(args.nontarget)}
    groups = pd.read_csv(args.groups, sep="\t", dtype=str).set_index("id")["group"].to_dict()
    members = read_membership(args.members)
    # MOTUs of the sequences we have (amplicons keep their original IDs).
    motus = {}
    for n in facet:
        motus.setdefault(members.get(n, n), []).append(n)
    if not facet or not nontarget:
        sys.exit("no facetotectan or no non-target sequences")

    # Random folds for the non-target sequences.
    rng = random.Random(args.seed)
    nt_ids = list(nontarget)
    rng.shuffle(nt_ids)
    folds = [nt_ids[i::args.folds] for i in range(args.folds)]

    rows = []
    for i, (motu, seqs) in enumerate(sorted(motus.items())):
        fold = folds[i % args.folds]
        query_nt = fold if i < args.folds else []
        held_out = set(seqs) | set(fold)
        d = os.path.join(args.workdir, "round")
        os.makedirs(d, exist_ok=True)
        ref = os.path.join(d, "ref.fasta")
        qry = os.path.join(d, "query.fasta")
        # Reference with SINTAX labels in the headers (BLAST ignores them).
        ref_recs = [(f"{n};tax=k:Arthropoda,o:Facetotecta;", s) for n, s in facet.items() if n not in held_out]
        ref_recs += [(f"{n};tax=k:Arthropoda,o:{groups.get(n, 'Other')};", s)
                     for n, s in nontarget.items() if n not in held_out]
        if not any("o:Facetotecta" in r[0] for r in ref_recs):
            continue  # a MOTU containing every facetotectan sequence: nothing left to test
        write_fasta(ref_recs, ref)
        write_fasta([(n, facet[n]) for n in seqs] + [(n, nontarget[n]) for n in query_nt], qry)
        bh = blast_best(qry, ref, d, args.threads)
        sx = sintax(qry, ref, d, args.threads, args.seed)
        for n in seqs + query_nt:
            hit = bh.get(n)
            rows.append(dict(
                query=n, truth="Facetotecta" if n in facet else groups.get(n, "Other"),
                positive=n in facet, held_out_motu=motu,
                blast_hit=hit[0].split(";")[0] if hit else "none",
                blast_hit_facet=bool(hit and "o:Facetotecta" in hit[0]) if hit else False,
                blast_pident=hit[1] if hit else 0.0, blast_qcov=round(hit[2], 3) if hit else 0.0,
                sintax_label=sx.get(n, ("unassigned", 0))[0], sintax_conf=sx.get(n, ("", 0.0))[1]))
        print(f"round {i + 1}/{len(motus)} ({motu}): {len(seqs)} positives, {len(query_nt)} negatives",
              file=sys.stderr)

    hits = pd.DataFrame(rows)
    hits.to_csv(args.out_hits, sep="\t", index=False)
    pos, neg = hits[hits.positive], hits[~hits.positive]

    # ROC tables.
    roc = []
    for t in np.arange(args.tmin, args.tmax + 0.001, 0.5):
        tp = ((pos.blast_hit_facet) & (pos.blast_pident >= t)).mean() if len(pos) else np.nan
        fp = ((neg.blast_hit_facet) & (neg.blast_pident >= t)).mean() if len(neg) else np.nan
        roc.append(("blast_identity", round(t, 1), tp, fp))
    for c in np.arange(0, 1.0001, 0.05):
        tp = ((pos.sintax_label == "Facetotecta") & (pos.sintax_conf >= c)).mean() if len(pos) else np.nan
        fp = ((neg.sintax_label == "Facetotecta") & (neg.sintax_conf >= c)).mean() if len(neg) else np.nan
        roc.append(("sintax_confidence", round(c, 2), tp, fp))
    roc = pd.DataFrame(roc, columns=["method", "threshold", "tpr", "fpr"])
    roc.to_csv(args.out_roc, sep="\t", index=False)

    # Chosen thresholds: lowest threshold with FPR <= target (highest TPR).
    out = []
    for method, sub in roc.groupby("method"):
        ok = sub[sub.fpr <= args.target_fpr].sort_values("threshold")
        if len(ok):
            r = ok.iloc[0]
            out.append((method, "chosen_fpr<=%.3g" % args.target_fpr, r.threshold, r.tpr, r.fpr))
    for t in (90.0, 92.0):
        r = roc[(roc.method == "blast_identity") & (roc.threshold == t)]
        if len(r):
            out.append(("blast_identity", "manuscript_v1", t, r.iloc[0].tpr, r.iloc[0].fpr))
    th = pd.DataFrame(out, columns=["method", "criterion", "threshold", "tpr", "fpr"])
    th.to_csv(args.out_thresholds, sep="\t", index=False)

    # Plots: ROC curves and rates against BLAST identity.
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for method, sub in roc.groupby("method"):
        axes[0].plot(sub.fpr, sub.tpr, marker=".", label=method)
    axes[0].set_xlabel("False-positive rate")
    axes[0].set_ylabel("True-positive rate (unsampled MOTUs)")
    axes[0].legend(frameon=False, fontsize=8)
    b = roc[roc.method == "blast_identity"]
    axes[1].plot(b.threshold, b.tpr, label="TPR")
    axes[1].plot(b.threshold, b.fpr, label="FPR")
    axes[1].set_xlabel("BLAST identity threshold (%)")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle(args.label)
    fig.tight_layout()
    fig.savefig(args.out_pdf)

    m = {f"{args.label}_positives": len(pos), f"{args.label}_negatives": len(neg)}
    for r in th.itertuples():
        key = f"{args.label}_{r.method}_{r.criterion.split('_')[0]}_{r.threshold}"
        m[f"{key}_tpr"] = round(r.tpr, 4)
        m[f"{key}_fpr"] = round(r.fpr, 4)
    write_metrics(args.metrics, m)


if __name__ == "__main__":
    main()
