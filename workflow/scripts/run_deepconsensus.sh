#!/usr/bin/env bash
# =============================================================================
# run_deepconsensus.sh: OPTIONAL DeepConsensus error correction (Baid et al. 2023).
#
# DeepConsensus re-calls each CCS read from its SUBREADS, so it needs the
# subread BAM. It was trained on Sequel II data; check that your chemistry
# is supported before using it. Steps (from the DeepConsensus quick start):
#   1. ccs: CCS reads from the subreads (BAM, with the per-read information)
#   2. actc: align subreads to their own CCS read
#   3. deepconsensus run: polish
# Only reads that survived our cleaning steps (keep list) are written out.
#
# Usage: run_deepconsensus.sh <subreads.bam> <cleaned_reads.fastq.gz>
#                             <checkpoint> <workdir> <threads> <out.fastq.gz>
# =============================================================================

set -euo pipefail

subreads=$1; keep=$2; ckpt=$3; work=$4; threads=$5; out=$6
mkdir -p "$work"

# 1. CCS from subreads (BAM output keeps per-read information).
ccs --min-rq=0.88 -j "$threads" "$subreads" "$work/ccs.bam"

# 2. Align subreads to their CCS read.
actc -j "$threads" "$subreads" "$work/ccs.bam" "$work/subreads_to_ccs.bam"

# 3. Polish.
deepconsensus run --subreads_to_ccs="$work/subreads_to_ccs.bam" \
    --ccs_bam="$work/ccs.bam" --checkpoint="$ckpt" \
    --output="$work/deepconsensus.fastq" --cpus "$threads"

# 4. Keep only reads that passed our filters (names from the cleaned FASTQ).
zcat -f "$keep" | awk 'NR%4==1{sub(/^@/,""); print $1}' > "$work/keep.txt"
python - "$work/deepconsensus.fastq" "$work/keep.txt" "$out" <<'EOF'
# Filter the polished FASTQ to the names in keep.txt and gzip the result.
import gzip, sys
fq, keepf, out = sys.argv[1:]
keep = set(l.strip() for l in open(keepf))
with open(fq) as fh, gzip.open(out, "wt") as o:
    while True:
        h = fh.readline()
        if not h:
            break
        s, p, q = fh.readline(), fh.readline(), fh.readline()
        if h[1:].split()[0] in keep:
            o.write(h + s + p + q)
EOF
