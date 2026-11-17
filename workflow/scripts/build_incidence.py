#!/usr/bin/env python3
# =============================================================================
# build_incidence.py: iNEXT input tables from MOTU membership + specimen data.
#
# Abundance table: one row per MOTU, one column per site plus "All"; values =
#   number of specimens of that MOTU (each Sanger specimen is one individual).
# Incidence-frequency table: first row = number of sampling events at the
#   site; following rows = number of events in which each MOTU was found.
#   (iNEXT datatype "incidence_freq".)
# Sequence IDs are "SAN-<specimen_id>"; specimen_id must be in specimens.tsv.
# =============================================================================

import argparse
import sys
from collections import defaultdict

import pandas as pd

from yfutils import read_membership


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--members", required=True)
    ap.add_argument("--specimens", required=True)
    ap.add_argument("--out-abundance", required=True)
    ap.add_argument("--out-incidence", required=True)
    args = ap.parse_args()

    members = read_membership(args.members)
    spec = pd.read_csv(args.specimens, sep="\t", comment="#", dtype=str).set_index("specimen_id")

    abund = defaultdict(lambda: defaultdict(int))       # motu -> site -> count
    events = defaultdict(set)                           # site -> sampling events
    found = defaultdict(lambda: defaultdict(set))       # site -> motu -> events
    missing = 0
    for seq, motu in members.items():
        sid = seq.split("-", 1)[1] if seq.startswith("SAN-") else seq
        if sid not in spec.index:
            missing += 1
            site, event = "unknown", "unknown"
        else:
            site = spec.loc[sid, "site"]
            event = spec.loc[sid, "sampling_event"]
        abund[motu][site] += 1
        abund[motu]["All"] += 1
        for s in (site, "All"):
            events[s].add(f"{site}|{event}")
            found[s][motu].add(f"{site}|{event}")
    if missing:
        print(f"WARNING: {missing} sequences not found in {args.specimens}", file=sys.stderr)

    sites = sorted(s for s in events if s != "All") + ["All"]
    motus = sorted(abund)
    a = pd.DataFrame([[abund[m][s] for s in sites] for m in motus], index=motus, columns=sites)
    a.index.name = "motu"
    a.to_csv(args.out_abundance, sep="\t")

    rows = [[len(events[s]) for s in sites]]
    rows += [[len(found[s][m]) for s in sites] for m in motus]
    inc = pd.DataFrame(rows, index=["n_sampling_events"] + motus, columns=sites)
    inc.index.name = "motu"
    inc.to_csv(args.out_incidence, sep="\t")


if __name__ == "__main__":
    main()
