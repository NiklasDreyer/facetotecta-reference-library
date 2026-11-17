#!/usr/bin/env python3
# =============================================================================
# authenticate_filter.py: iterative phylogenetic filter for high-throughput
# sequences (manuscript_v1, "Sequence authenticity & filtering").
#
# Each round:
#   1. align ingroup + outgroups with MAFFT (--auto --maxiterate 1000)
#   2. trim with Gblocks (less stringent settings; run_gblocks.sh)
#   3. ML tree with IQ-TREE 2: ModelFinder (-m MFP), UFBoot with -bnni,
#      SH-aLRT (manuscript_v1: 10,000 replicates each)
#   4. root the tree on the outgroups
#   5. remove HIGH-THROUGHPUT sequences (PB-, ILL-, PUB-) that
#      (1) fall outside Facetotecta, defined as the smallest clade containing
#          all Sanger (SAN-) sequences
#      (2) sit on a terminal branch more than --long-branch times longer than
#          the branches leading to their sister's tips (Yang & Smith 2014)
#      (3) belong to a long-branch clade (subtending branch more than
#          --long-branch times the median internal branch) that contains no
#          Sanger sequence and reads from only one high-throughput platform
#      (4) optional (--require-cross-platform): whose closest well-supported
#          clade (SH-aLRT >= --min-shalrt and UFBoot >= --min-ufboot) contains
#          no sequence from any other platform
#   Sanger sequences are never removed (they were proofed by hand).
# Rounds repeat until nothing is removed or --max-iter is reached.
# The support of the Facetotecta clade is reported each round; manuscript_v1
# requires > 80% SH-aLRT and > 95% UFBoot.
# =============================================================================

import argparse
import os
import shutil
import statistics
import subprocess
import sys

from ete3 import Tree

from yfutils import platform_of, read_fasta, write_fasta, write_metrics

SCRIPTS = os.path.dirname(os.path.abspath(__file__))


def run(cmd, log):
    """Run a shell command, appending its output to the round's log."""
    with open(log, "a") as fh:
        fh.write(f"\n$ {cmd}\n")
        fh.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=fh, stderr=fh)


def support(node):
    """(SH-aLRT, UFBoot) from an IQ-TREE node label such as '85.2/97'."""
    try:
        parts = [float(x) for x in str(node.name).split("/") if x]
    except ValueError:
        return None, None
    if len(parts) == 2:
        return parts[0], parts[1]
    if len(parts) == 1:
        return None, parts[0]
    return None, None


def root_on_outgroups(tree, outgroups):
    """Root the tree on the outgroup clade (or on one outgroup if they are
    not monophyletic). Returns True if the outgroups form a clade."""
    leaves = [l for l in tree.get_leaves() if l.name in outgroups]
    if not leaves:
        sys.exit("no outgroup sequences found in the tree")
    # Root on one outgroup first so the outgroup MRCA is well defined.
    tree.set_outgroup(leaves[0])
    if len(leaves) == 1:
        return True
    mrca = tree.get_common_ancestor(leaves)
    if mrca is not tree:
        tree.set_outgroup(mrca)
    mono = set(mrca.get_leaf_names()) <= set(outgroups)
    return mono


def round_filter(tree, outgroups, args):
    """Apply the removal criteria to one rooted tree.
    Returns (dict id -> criterion, facet support tuple, warnings list)."""
    removed, warnings = {}, []
    san = [l for l in tree.get_leaves() if platform_of(l.name) == "SAN"]
    if not san:
        sys.exit("no Sanger (SAN-) sequences in the tree")
    facet = tree.get_common_ancestor(san) if len(san) > 1 else san[0].up
    facet_leaves = set(facet.get_leaf_names())
    if facet_leaves & set(outgroups):
        warnings.append("Sanger sequences are not monophyletic with respect to the outgroups")
    facet_sh, facet_ub = support(facet)
    if facet_sh is not None and facet_ub is not None and \
            (facet_sh < args.min_shalrt or facet_ub < args.min_ufboot):
        warnings.append(f"Facetotecta clade support below threshold: SH-aLRT {facet_sh}, UFBoot {facet_ub}")

    hts = lambda name: platform_of(name) in {"PB", "ILL", "PUB"}

    # (1) HTS sequences outside Facetotecta.
    for leaf in tree.get_leaves():
        if hts(leaf.name) and leaf.name not in facet_leaves:
            removed[leaf.name] = "outside_facetotecta"

    # Branch-length scales inside Facetotecta (to avoid dividing by ~0).
    tips = [l.dist for l in facet.get_leaves()]
    floor = max(statistics.median(tips) if tips else 0, 1e-6)
    internal = [n.dist for n in facet.traverse() if not n.is_leaf() and n is not facet and n.dist > 0]
    int_median = max(statistics.median(internal) if internal else 0, 1e-6)

    # (2) Long terminal branches compared with the sister lineage.
    for leaf in facet.get_leaves():
        if not hts(leaf.name) or leaf.name in removed:
            continue
        sisters = leaf.get_sisters()
        if not sisters:
            continue
        ref = []
        for s in sisters:
            # Mean path length from the shared parent to the sister's tips.
            ref += [s.dist + s.get_distance(t) for t in s.get_leaves()] if not s.is_leaf() else [s.dist]
        ref_len = max(statistics.mean(ref), floor)
        if leaf.dist > args.long_branch * ref_len:
            removed[leaf.name] = "long_terminal_branch"

    # (3) Long-branch clades without Sanger sequences and from one HTS platform.
    for node in facet.traverse():
        if node.is_leaf() or node is facet:
            continue
        if node.dist <= args.long_branch * int_median:
            continue
        names = node.get_leaf_names()
        plats = {platform_of(n) for n in names}
        if "SAN" not in plats and len(plats) == 1:
            for n in names:
                removed.setdefault(n, "long_branch_single_platform_clade")

    # (4) Optional: no other platform in the closest supported clade.
    if args.require_cross_platform:
        for leaf in facet.get_leaves():
            if not hts(leaf.name) or leaf.name in removed:
                continue
            node = leaf.up
            while node is not None and node is not facet:
                sh, ub = support(node)
                if sh is not None and ub is not None and sh >= args.min_shalrt and ub >= args.min_ufboot:
                    break
                node = node.up
            node = node or facet
            others = {platform_of(n) for n in node.get_leaf_names()} - {platform_of(leaf.name), "OUT"}
            if not others:
                removed[leaf.name] = "no_cross_platform_support"

    return removed, (facet_sh, facet_ub), warnings


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="ingroup sequences (SAN-, PB-, ILL-, PUB-)")
    ap.add_argument("--outgroups", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--min-shalrt", type=float, default=80)
    ap.add_argument("--min-ufboot", type=float, default=95)
    ap.add_argument("--long-branch", type=float, default=10)
    ap.add_argument("--max-iter", type=int, default=10)
    ap.add_argument("--ufboot", type=int, default=10000)
    ap.add_argument("--alrt", type=int, default=10000)
    ap.add_argument("--require-cross-platform", action="store_true")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-removed", required=True)
    ap.add_argument("--out-tree", required=True)
    ap.add_argument("--out-aln", required=True)
    ap.add_argument("--metrics")
    ap.add_argument("--label", default="locus")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    ingroup = {n: s.replace("-", "") for n, _, s in read_fasta(args.input) if platform_of(n) != "OUT"}
    # Outgroup IDs must start with "OUT-"; add the tag if missing.
    outgroups = {}
    for n, _, s in read_fasta(args.outgroups):
        outgroups[n if n.startswith("OUT-") else f"OUT-{n}"] = s.replace("-", "")
    start_counts = {}
    for n in ingroup:
        start_counts[platform_of(n)] = start_counts.get(platform_of(n), 0) + 1

    all_removed = []          # (id, round, criterion)
    converged = False
    last = None
    facet_support = (None, None)
    for it in range(1, args.max_iter + 1):
        d = os.path.join(args.workdir, f"round{it:02d}")
        os.makedirs(d, exist_ok=True)
        log = os.path.join(d, "round.log")
        seqs = os.path.join(d, "input.fasta")
        write_fasta(list(ingroup.items()) + list(outgroups.items()), seqs)

        # 1-3: align, trim, tree.
        aln = os.path.join(d, "aln.fasta")
        gb = os.path.join(d, "aln.gb.fasta")
        run(f"mafft --auto --maxiterate 1000 --thread {args.threads} {seqs} > {aln}", log)
        run(f"bash {SCRIPTS}/run_gblocks.sh {aln} {gb} d", log)
        prefix = os.path.join(d, "iqtree")
        run(f"iqtree2 -s {gb} -m MFP -B {args.ufboot} -alrt {args.alrt} -bnni "
            f"-T {args.threads} --prefix {prefix} -redo", log)

        # 4: root and filter.
        tree = Tree(f"{prefix}.treefile", format=1)
        mono = root_on_outgroups(tree, set(outgroups))
        removed, facet_support, warnings = round_filter(tree, set(outgroups), args)
        if not mono:
            warnings.append("outgroups are not monophyletic")
        for w in warnings:
            print(f"round {it}: WARNING: {w}", file=sys.stderr)
        print(f"round {it}: {len(ingroup)} ingroup sequences, {len(removed)} removed; "
              f"Facetotecta support SH-aLRT/UFBoot = {facet_support}", file=sys.stderr)
        last = (prefix, gb, tree)

        for name, crit in removed.items():
            all_removed.append((name, it, crit))
            ingroup.pop(name, None)
        if not removed:
            converged = True
            break

    # Outputs from the last round.
    prefix, gb, tree = last
    write_fasta(ingroup.items(), args.out_fasta)
    tree.write(outfile=args.out_tree, format=1)
    shutil.copyfile(gb, args.out_aln)
    with open(args.out_removed, "w") as fh:
        fh.write("id\tround\tcriterion\n")
        for name, it, crit in all_removed:
            fh.write(f"{name}\t{it}\t{crit}\n")

    m = {f"auth_{args.label}_rounds": it, f"auth_{args.label}_converged": "yes" if converged else "no",
         f"auth_{args.label}_facet_shalrt": facet_support[0], f"auth_{args.label}_facet_ufboot": facet_support[1]}
    for plat, n in start_counts.items():
        kept = sum(1 for x in ingroup if platform_of(x) == plat)
        m[f"auth_{args.label}_{plat}_input"] = n
        m[f"auth_{args.label}_{plat}_kept"] = kept
    for crit in sorted({c for _, _, c in all_removed}):
        m[f"auth_{args.label}_removed_{crit}"] = sum(1 for _, _, c in all_removed if c == crit)
    write_metrics(args.metrics, m)


if __name__ == "__main__":
    main()
