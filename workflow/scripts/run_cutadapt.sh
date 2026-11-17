#!/usr/bin/env bash
# =============================================================================
# run_cutadapt.sh: remove the original study's primers from public reads.
#
# Paired-end: -g forward primer on R1, -G reverse primer on R2; pairs in which
# either primer is missing are discarded (--discard-untrimmed), so only real
# amplicons of the expected marker remain.
# Single-end: forward primer only; an empty R2 file is created so the
# workflow has the same outputs for both layouts (DADA2 then treats the
# dataset as single-end).
#
# Usage: run_cutadapt.sh <PE|SE> <fwd> <rev> <threads> <out_R1> <out_R2> <in_R1> [in_R2]
# =============================================================================

set -euo pipefail

layout=$1; fwd=$2; rev=$3; threads=$4; o1=$5; o2=$6; i1=$7; i2=${8:-}

if [ "$layout" = "PE" ]; then
    cutadapt -g "$fwd" -G "$rev" --discard-untrimmed -j "$threads" \
        -o "$o1" -p "$o2" "$i1" "$i2"
else
    cutadapt -g "$fwd" --discard-untrimmed -j "$threads" -o "$o1" "$i1"
    # Empty placeholder for the second read file (size 0, which dada2_asv.R
    # recognises as "single-end").
    : > "$o2"
fi
