#!/usr/bin/env python3
# =============================================================================
# mining_summary.py: results of the public-data mining.
#
# A public ASV counts as a CONFIRMED facetotectan detection when it
#   1. passed the BLAST cut-offs with a Facetotecta top hit (assign_two_pass.py)
#   2. was placed inside Facetotecta with LWR >= --min-lwr (EPA-ng + gappa)
#   3. survived the phylogenetic authentication filter
# Confirmed ASVs are grouped into MOTUs together with our library; a MOTU
# with no library sequence is a lineage NEW to the library.
#
# Outputs:
#   gain_by_dataset.tsv  ASVs, candidates, confirmed detections, what they
#                        were before our library (unassigned / other taxon /
#                        already Facetotecta), MOTUs and new MOTUs per dataset
#   detections.tsv       one row per confirmed ASV, with coordinates
#   detections_map.pdf   simple longitude/latitude plot of the datasets
#   mining.metrics.tsv   totals per marker
# =============================================================================

import argparse
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from yfutils import platform_of, read_fasta, read_membership, write_metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assign", nargs="+", required=True)
    ap.add_argument("--placement", nargs="+", required=True)
    ap.add_argument("--authentic", nargs="+", required=True)
    ap.add_argument("--members", nargs="+", required=True)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--min-lwr", type=float, default=0.9)
    ap.add_argument("--out-gain", required=True)
    ap.add_argument("--out-detections", required=True)
    ap.add_argument("--out-pdf", required=True)
    ap.add_argument("--metrics")
    args = ap.parse_args()

    ds = pd.read_csv(args.datasets, sep="\t", comment="#", dtype=str).set_index("dataset")

    # Placement result per public ASV (IDs "PUB-<dataset>-ASV<n>").
    placed = {}
    for p in args.placement:
        a = pd.read_csv(p, sep="\t", dtype={"query": str})
        for r in a.itertuples():
            placed[r.query] = ("Facetotecta" in str(r.taxopath)) and float(r.LWR) >= args.min_lwr
    # Authentic sequences and MOTU membership.
    authentic = {}
    for p in args.authentic:
        authentic.update({n: s for n, _, s in read_fasta(p)})
    members = {}
    for p in args.members:
        members.update(read_membership(p))
    # MOTUs containing at least one library (non-PUB) sequence.
    library_motus = {m for s, m in members.items() if platform_of(s) != "PUB"}

    detections, gain_rows = [], []
    for path in args.assign:
        dataset = os.path.basename(os.path.dirname(path))
        t = pd.read_csv(path, sep="\t", dtype=str)
        cand = t[t["candidate"] == "yes"]
        confirmed = []
        for r in cand.itertuples():
            pid = f"PUB-{r.asv}"
            if placed.get(pid) and pid in authentic:
                motu = members.get(pid, "NA")
                confirmed.append(dict(
                    dataset=dataset, marker=ds.loc[dataset, "marker"], region=ds.loc[dataset, "region"],
                    lat=ds.loc[dataset, "lat"], lon=ds.loc[dataset, "lon"], asv=r.asv,
                    identity=r.augmented_identity, before=r.before, motu=motu,
                    new_lineage="yes" if motu not in library_motus else "no",
                    sequence=authentic[pid].replace("-", "")))
        detections += confirmed
        c = pd.DataFrame(confirmed)
        before_counts = c["before"].value_counts().to_dict() if len(c) else {}
        gain_rows.append(dict(
            dataset=dataset, marker=ds.loc[dataset, "marker"], region=ds.loc[dataset, "region"],
            asvs=len(t), candidates=len(cand), confirmed=len(c),
            was_unassigned=before_counts.get("was_unassigned", 0),
            was_other_taxon=before_counts.get("was_other_taxon", 0),
            was_facetotecta=before_counts.get("was_facetotecta", 0),
            motus=c["motu"].nunique() if len(c) else 0,
            new_motus=c.loc[c.new_lineage == "yes", "motu"].nunique() if len(c) else 0))

    gain = pd.DataFrame(gain_rows)
    gain.to_csv(args.out_gain, sep="\t", index=False)
    det = pd.DataFrame(detections, columns=["dataset", "marker", "region", "lat", "lon", "asv",
                                            "identity", "before", "motu", "new_lineage", "sequence"])
    det.to_csv(args.out_detections, sep="\t", index=False)

    # Map: one point per dataset, filled if Facetotecta were detected.
    fig, ax = plt.subplots(figsize=(8, 4))
    for r in gain.itertuples():
        lat, lon = float(ds.loc[r.dataset, "lat"]), float(ds.loc[r.dataset, "lon"])
        ax.scatter(lon, lat, s=20 + 10 * r.motus, facecolors="#C44E52" if r.confirmed else "none",
                   edgecolors="#C44E52" if r.confirmed else "grey")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Public datasets screened (filled = Facetotecta detected; size = MOTUs)")
    fig.tight_layout()
    fig.savefig(args.out_pdf)

    m = {"mining_datasets_screened": len(gain),
         "mining_datasets_with_detections": int((gain.confirmed > 0).sum()) if len(gain) else 0}
    for marker, sub in det.groupby("marker"):
        m[f"mining_{marker}_confirmed_asvs"] = len(sub)
        m[f"mining_{marker}_motus"] = sub["motu"].nunique()
        m[f"mining_{marker}_new_motus"] = sub.loc[sub.new_lineage == "yes", "motu"].nunique()
        m[f"mining_{marker}_was_unassigned"] = int((sub.before == "was_unassigned").sum())
        m[f"mining_{marker}_was_other_taxon"] = int((sub.before == "was_other_taxon").sum())
        m[f"mining_{marker}_regions"] = ";".join(sorted(sub.region.unique()))
    write_metrics(args.metrics, m)


if __name__ == "__main__":
    main()
