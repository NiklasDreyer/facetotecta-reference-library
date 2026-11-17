#!/usr/bin/env bash
# =============================================================================
# run_metamate.sh: NUMT / error filtering of COX1 ASVs with metaMATE
# (Andujar et al. 2021, Molecular Ecology Resources).
#
# What metaMATE does: it labels ASVs as "verified authentic" (match a curated
# reference, here our Sanger COX1 library) or "verified non-authentic"
# (wrong length or stop codons), then tests a grid of read-abundance
# thresholds (applied per library, per clade, in total...) and estimates for
# each threshold set how many authentic and non-authentic ASVs it removes.
# One threshold set is then chosen and applied with "metamate dump".
#
# Usage:
#   run_metamate.sh <asvs.fasta> <asv_table.tsv> <reference.fasta> <specs.txt>
#                   <expected_length> <index|auto> <outdir> <threads>
#                   <out.fasta> <metrics.tsv>
#
# CHECK on the first run: metaMATE's command-line flags and file formats
# have changed between versions. Run "metamate find -h" and "metamate dump -h"
# in the metamate environment and adjust the two commands below if needed.
# =============================================================================

set -euo pipefail

asvs=$1; table=$2; ref=$3; specs=$4; length=$5; index=$6
outdir=$7; threads=$8; out=$9; metrics=${10}
mkdir -p "$outdir"
here=$(dirname "$0")

# 1. Read map: metaMATE expects read counts per library (rows) x ASV (columns),
#    i.e. the transpose of our ASV x sample table.
python "$here/metamate_tools.py" readmap --table "$table" --output "$outdir/readmap.csv"

# 2. Copy the curated Sanger COX1 library as the metaMATE reference set.
cp "$ref" "$outdir/references.fasta"

# 3. Test the threshold grid.
metamate find -A "$asvs" -M "$outdir/readmap.csv" -S "$specs" \
    -R "$outdir/references.fasta" -L "$length" --table 5 \
    -o "$outdir/find" -t "$threads" --overwrite

# 4. Choose a threshold set (or use the index given in the config).
chosen=$(python "$here/metamate_tools.py" choose --results "$outdir/find/results.csv" --index "$index")
echo "metaMATE threshold set chosen: $chosen"

# 5. Apply it.
metamate dump -A "$asvs" -C "$outdir/find/resultcache" -o "$outdir/dump" -i "$chosen" --overwrite

# 6. The kept ASVs: metaMATE writes a FASTA in the dump folder.
kept=$(find "$outdir/dump" -name "*.fasta" | head -n 1)
cp "$kept" "$out"

# 7. Metrics.
n_in=$(grep -c '^>' "$asvs")
n_out=$(grep -c '^>' "$out")
printf "illumina_COX1_metamate_input\t%s\nillumina_COX1_metamate_kept\t%s\nillumina_COX1_metamate_result_index\t%s\n" \
    "$n_in" "$n_out" "$chosen" > "$metrics"
