#!/usr/bin/env bash
# =============================================================================
# run_gblocks.sh: run Gblocks 0.91b with the "less stringent" settings.
#
# Less stringent settings (as offered on the Gblocks web server):
#   -b1 = minimum number of sequences for a conserved position = n/2 + 1 (default)
#   -b2 = minimum number for a flanking position = b1 (less strict flanks)
#   -b3 = 8  maximum number of contiguous non-conserved positions (default)
#   -b4 = 5  minimum block length (smaller final blocks allowed)
#   -b5 = h  gap positions allowed if present in up to half the sequences
#
# Usage: run_gblocks.sh <in.aln.fasta> <out.fasta> <d|c|p>
#   d = DNA, c = codons (removes whole codons; keeps the reading frame), p = protein
# =============================================================================

set -euo pipefail

in=$1
out=$2
type=${3:-d}

# Number of sequences in the alignment.
n=$(grep -c '^>' "$in")
# b1 = floor(n/2) + 1, and b2 = b1 for less strict flanking positions.
b1=$(( n / 2 + 1 ))

# Gblocks writes "<input>-gb" next to the input; work on a copy in a temp dir.
tmp=$(mktemp -d)
cp "$in" "$tmp/aln.fasta"

# Gblocks exits with status 1 even when it succeeds, so do not stop on it.
Gblocks "$tmp/aln.fasta" -t="$type" -b1="$b1" -b2="$b1" -b3=8 -b4=5 -b5=h -e=-gb || true

if [ ! -s "$tmp/aln.fasta-gb" ]; then
    echo "Gblocks produced no output for $in" >&2
    exit 1
fi

# Gblocks writes sequences in blocks with spaces: remove the spaces.
sed '/^>/! s/ //g' "$tmp/aln.fasta-gb" > "$out"

# Keep the HTML report next to the output for inspection.
[ -f "$tmp/aln.fasta-gb.htm" ] && cp "$tmp/aln.fasta-gb.htm" "${out}.htm"
rm -rf "$tmp"
echo "Gblocks done: $n sequences, b1=b2=$b1, type=$type"
