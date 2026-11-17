#!/usr/bin/env python3
# =============================================================================
# run_asap.py: run ASAP (Puillandre et al. 2021) and read its best partition.
#
# The ASAP command-line program is not on conda (see install.sh). The command
# template comes from config.yaml (delimitation.asap.command), with {fasta}
# and {outdir} filled in. ASAP writes one file per ranked partition; the file
# matching --partition-glob (the best-ranked one) is read. It is expected to
# have two columns: sequence name and group number (comma, tab or space
# separated). CHECK the file names and format of your ASAP version on the
# first run and adjust the glob if needed.
# Output: sequence <tab> MOTU ("asap_<group>").
# =============================================================================

import argparse
import glob
import os
import re
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--command", required=True)
    ap.add_argument("--partition-glob", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    if os.path.getsize(args.fasta) == 0:
        open(args.output, "w").write("sequence\tmotu\n")
        return

    cmd = args.command.format(fasta=args.fasta, outdir=args.outdir)
    print(f"$ {cmd}", file=sys.stderr)
    subprocess.run(cmd, shell=True, check=True)

    files = sorted(glob.glob(os.path.join(args.outdir, args.partition_glob)))
    if not files:
        sys.exit(f"no ASAP partition file matching {args.partition_glob} in {args.outdir}")
    print(f"reading {files[0]}", file=sys.stderr)

    rows = []
    with open(files[0]) as fh:
        for line in fh:
            parts = [p for p in re.split(r"[,\t ]+", line.strip()) if p]
            if len(parts) < 2 or not parts[-1].isdigit():
                continue  # skip headers and comments
            rows.append((parts[0], f"asap_{parts[-1]}"))

    with open(args.output, "w") as out:
        out.write("sequence\tmotu\n")
        for seq, motu in rows:
            out.write(f"{seq}\t{motu}\n")


if __name__ == "__main__":
    main()
