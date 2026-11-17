#!/usr/bin/env python3
# =============================================================================
# assign_two_pass.py: assign ASVs of a public dataset twice.
#
#   pass 1 "public":    BLAST against the public reference only
#                       (what could be found WITHOUT our library)
#   pass 2 "augmented": BLAST against the public reference + our library
#
# For each ASV the best hit (highest bitscore) of each pass is recorded with
# its lineage, identity and e-value. An ASV is a CANDIDATE facetotectan
# sequence when, in pass 2, its top hit is Facetotecta with identity
# >= --identity and e-value <= --evalue (manuscript_v1 cut-offs; replace
# them with the benchmark thresholds). The "gain" columns record what each
# candidate was in pass 1: unassigned (no hit), assigned to another taxon, or
# already Facetotecta (sequences of the few published facetotectans).
# Candidates are written as "<dataset>-ASV<n>" for placement and filtering.
# =============================================================================

import argparse
import os
import subprocess

from yfutils import read_fasta, read_taxonomy, write_fasta, write_metrics


def blast(query, db, threads, out):
    """BLAST query against db; return {query: (subject, pident, evalue)} best by bitscore."""
    subprocess.run(["blastn", "-query", query, "-db", db, "-max_target_seqs", "20",
                    "-num_threads", str(threads), "-evalue", "1e-5",
                    "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
                    "-out", out], check=True)
    best = {}
    for line in open(out):
        q, s, pid, _, ev, bits = line.rstrip("\n").split("\t")
        s = s.replace("lcl|", "")  # some BLAST versions prefix local IDs
        bits = float(bits)
        if q not in best or bits > best[q][3]:
            best[q] = (s, float(pid), float(ev), bits)
    return best


def short_lineage(lineage, depth=4):
    """Last few ranks of a lineage, for readable tables."""
    parts = [p for p in lineage.split(";") if p]
    return ";".join(parts[-depth:]) if parts else "none"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--db-public", required=True)
    ap.add_argument("--tax-public", required=True)
    ap.add_argument("--db-augmented", required=True)
    ap.add_argument("--tax-augmented", required=True)
    ap.add_argument("--identity", type=float, required=True)
    ap.add_argument("--evalue", type=float, required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--metrics")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    asvs = {n: s for n, _, s in read_fasta(args.fasta)}
    tax_pub, tax_aug = read_taxonomy(args.tax_public), read_taxonomy(args.tax_augmented)
    p1 = blast(args.fasta, args.db_public, args.threads, os.path.join(args.workdir, "blast_public.tsv"))
    p2 = blast(args.fasta, args.db_augmented, args.threads, os.path.join(args.workdir, "blast_augmented.tsv"))

    candidates = []
    gain = {"was_unassigned": 0, "was_other_taxon": 0, "was_facetotecta": 0}
    with open(args.out_table, "w") as out:
        out.write("asv\tpublic_hit\tpublic_identity\tpublic_lineage\t"
                  "augmented_hit\taugmented_identity\taugmented_evalue\taugmented_lineage\t"
                  "candidate\tbefore\n")
        for n in asvs:
            h1, h2 = p1.get(n), p2.get(n)
            lin1 = tax_pub.get(h1[0], "") if h1 else ""
            lin2 = tax_aug.get(h2[0], "") if h2 else ""
            is_facet = bool(h2) and "Facetotecta" in lin2.split(";")
            cand = is_facet and h2[1] >= args.identity and h2[2] <= args.evalue
            before = "NA"
            if cand:
                if not h1:
                    before = "was_unassigned"
                elif "Facetotecta" in lin1.split(";"):
                    before = "was_facetotecta"
                else:
                    before = "was_other_taxon"
                gain[before] += 1
                candidates.append((n, asvs[n]))
            out.write(f"{n}\t{h1[0] if h1 else 'none'}\t{h1[1] if h1 else 'NA'}\t{short_lineage(lin1)}\t"
                      f"{h2[0] if h2 else 'none'}\t{h2[1] if h2 else 'NA'}\t{h2[2] if h2 else 'NA'}\t"
                      f"{short_lineage(lin2)}\t{'yes' if cand else 'no'}\t{before}\n")

    write_fasta(candidates, args.out_candidates)
    m = {f"mining_{args.dataset}_asvs": len(asvs), f"mining_{args.dataset}_candidates": len(candidates)}
    m.update({f"mining_{args.dataset}_candidates_{k}": v for k, v in gain.items()})
    write_metrics(args.metrics, m)


if __name__ == "__main__":
    main()
