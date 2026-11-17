#!/usr/bin/env python3
# =============================================================================
# check_translation.py: translation checks for protein-coding loci.
#
# Modes:
#   report : translate each (ungapped) sequence in all three forward frames
#            and report the frame with fewest stop codons. Sequences with
#            stops in every frame are flagged (possible NUMT, frameshift or
#            sequencing error). Nothing is removed.
#   frame  : for an ALIGNMENT, find the column offset (0, 1, 2) that gives
#            the fewest stop codons across all sequences, and trim the
#            alignment so it starts at codon position 1 and ends on a
#            complete codon. Needed for codon-aware Gblocks and codon partitions.
#   filter : keep only sequences with no internal stop codon. With
#            --frame-from-header the frame is read from a "refstart=<n>" field
#            in the header (0-based start on a reference that begins at codon
#            position 1, written by trim_to_reference.py). Otherwise the best
#            of the three forward frames is used.
#
# Genetic code tables are NCBI numbers (5 = invertebrate mitochondrial,
# 1 = standard). Biopython does the translation.
# =============================================================================

import argparse
import re
import sys

from Bio.Seq import Seq

from yfutils import read_fasta, write_fasta


def stops_in(seq, table, frame):
    """Count internal stop codons of seq translated from 'frame' (0-2).
    A stop in the last codon is allowed (complete gene ends)."""
    sub = seq[frame:]
    sub = sub[: len(sub) - len(sub) % 3]
    if len(sub) < 3:
        return 0
    # Ambiguous bases (N, R, Y...) translate to X and never count as stops.
    prot = str(Seq(sub).translate(table=table))
    return prot[:-1].count("*")


def best_frame(seq, table):
    """Frame (0-2) with fewest internal stops, and that number."""
    counts = [stops_in(seq, table, f) for f in range(3)]
    f = min(range(3), key=lambda i: counts[i])
    return f, counts[f]


def mode_report(args):
    """Write one line per sequence: id, length, best frame, stops."""
    n_bad = 0
    with open(args.output, "w") as out:
        out.write("id\tlength\tbest_frame\tstops_best_frame\tflag\n")
        for name, _, seq in read_fasta(args.input):
            seq = seq.upper().replace("-", "")
            f, n = best_frame(seq, args.table)
            flag = "STOP_IN_ALL_FRAMES" if n > 0 else "ok"
            n_bad += n > 0
            out.write(f"{name}\t{len(seq)}\t{f}\t{n}\t{flag}\n")
    print(f"{n_bad} sequences with stop codons in every frame", file=sys.stderr)


def mode_frame(args):
    """Find the alignment's reading frame and trim to whole codons."""
    recs = [(n, s.upper()) for n, _, s in read_fasta(args.input)]
    length = max(len(s) for _, s in recs)
    best = None
    for offset in range(3):
        total = 0
        for _, s in recs:
            # Read codons column-wise; skip codons containing gaps.
            for i in range(offset, length - 2, 3):
                codon = s[i:i + 3]
                if len(codon) < 3 or "-" in codon or "?" in codon:
                    continue
                if str(Seq(codon).translate(table=args.table)) == "*":
                    total += 1
        print(f"offset {offset}: {total} stop codons", file=sys.stderr)
        if best is None or total < best[1]:
            best = (offset, total)
    offset = best[0]
    end = offset + (length - offset) // 3 * 3
    write_fasta([(n, s.ljust(length, "-")[offset:end]) for n, s in recs], args.output)
    print(f"chosen offset {offset}; alignment trimmed to columns {offset + 1}-{end}", file=sys.stderr)


def mode_filter(args):
    """Keep sequences without internal stop codons."""
    kept, removed = [], 0
    rep = open(args.report, "w") if args.report else None
    if rep:
        rep.write("id\tframe\tstops\tkept\n")
    for name, desc, seq in read_fasta(args.input):
        seq = seq.upper().replace("-", "")
        if args.frame_from_header:
            m = re.search(r"refstart=(\d+)", desc)
            start = int(m.group(1)) if m else 0
            # Bases to skip so the first full codon begins on a codon boundary.
            frame = (3 - start % 3) % 3
            n = stops_in(seq, args.table, frame)
        else:
            frame, n = best_frame(seq, args.table)
        ok = n == 0
        if rep:
            rep.write(f"{name}\t{frame}\t{n}\t{'yes' if ok else 'no'}\n")
        if ok:
            kept.append((name, seq))
        else:
            removed += 1
    if rep:
        rep.close()
    write_fasta(kept, args.output)
    print(f"kept {len(kept)}, removed {removed} with internal stops", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["report", "frame", "filter"])
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--table", type=int, required=True, help="NCBI genetic code")
    ap.add_argument("--report", help="filter mode: per-sequence table")
    ap.add_argument("--frame-from-header", action="store_true")
    args = ap.parse_args()
    {"report": mode_report, "frame": mode_frame, "filter": mode_filter}[args.mode](args)


if __name__ == "__main__":
    main()
