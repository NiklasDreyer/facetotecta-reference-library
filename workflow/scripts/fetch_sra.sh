#!/usr/bin/env bash
# =============================================================================
# fetch_sra.sh: download one SRA run and write gzipped FASTQ.
#
#   prefetch      downloads the .sra file
#   fasterq-dump  converts it to FASTQ (--split-files: one file per read)
# Output names: <outdir>/<run>_R1.fastq.gz (and _R2 for paired-end).
#
# Usage: fetch_sra.sh <run accession> <outdir> <PE|SE>
# Runs on the login node (needs internet).
# =============================================================================

set -euo pipefail

run=$1; outdir=$2; layout=$3
mkdir -p "$outdir"
tmp=$(mktemp -d)

prefetch "$run" -O "$tmp"
fasterq-dump --split-files --outdir "$tmp" --temp "$tmp" "$tmp/$run"

if [ "$layout" = "PE" ]; then
    pigz -c "$tmp/${run}_1.fastq" > "$outdir/${run}_R1.fastq.gz"
    pigz -c "$tmp/${run}_2.fastq" > "$outdir/${run}_R2.fastq.gz"
else
    # Single-end runs are written as <run>.fastq or <run>_1.fastq depending on version.
    f=$(ls "$tmp/${run}"*.fastq | head -n 1)
    pigz -c "$f" > "$outdir/${run}_R1.fastq.gz"
fi
rm -rf "$tmp"
