#!/usr/bin/env python3
# =============================================================================
# metamate_tools.py: helpers for run_metamate.sh.
#
#   readmap : convert the DADA2 ASV x sample table (TSV) into the CSV read map
#             metaMATE expects (one row per library, one column per ASV).
#   choose  : choose a threshold set from metaMATE's results.csv. With an
#             integer --index, that index is used. With "auto", the row with
#             the highest accuracy score is chosen among rows that keep all
#             verified-authentic ASVs. Column names are matched loosely
#             because they differ between metaMATE versions; the chosen row
#             is printed to stderr so you can check it.
# =============================================================================

import argparse
import sys

import csv


def readmap(args):
    """Transpose ASV x sample counts to sample x ASV."""
    with open(args.table) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header, data = rows[0], rows[1:]
    samples = header[1:]
    asvs = [r[0] for r in data]
    with open(args.output, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow([""] + asvs)
        for j, s in enumerate(samples):
            w.writerow([s] + [r[j + 1] for r in data])


def choose(args):
    """Pick the index of the threshold set to apply."""
    if args.index != "auto":
        print(args.index)
        return
    with open(args.results) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit("results.csv is empty")
    cols = rows[0].keys()
    # Column with the overall accuracy score.
    acc_col = next((c for c in cols if "accuracy" in c.lower()), None)
    # Column counting verified-authentic ASVs that were rejected.
    va_rej = next((c for c in cols if "verifiedauthentic" in c.lower().replace("_", "").replace("-", "")
                   and ("reject" in c.lower() or "removed" in c.lower())), None)
    idx_col = next((c for c in cols if c.lower() in ("resultindex", "result_index", "index")), None)
    candidates = rows
    if va_rej:
        safe = [r for r in rows if float(r[va_rej] or 0) == 0]
        candidates = safe or rows
    if acc_col:
        best = max(candidates, key=lambda r: float(r[acc_col] or 0))
    else:
        best = candidates[0]
        print("WARNING: no accuracy column found; using the first row", file=sys.stderr)
    print(f"chosen row: {best}", file=sys.stderr)
    print(best[idx_col] if idx_col else rows.index(best))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="mode", required=True)
    r = sub.add_parser("readmap")
    r.add_argument("--table", required=True)
    r.add_argument("--output", required=True)
    c = sub.add_parser("choose")
    c.add_argument("--results", required=True)
    c.add_argument("--index", default="auto")
    args = ap.parse_args()
    readmap(args) if args.mode == "readmap" else choose(args)


if __name__ == "__main__":
    main()
