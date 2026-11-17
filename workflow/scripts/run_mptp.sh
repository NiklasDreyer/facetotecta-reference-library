#!/usr/bin/env bash
# =============================================================================
# run_mptp.sh: species delimitation with mPTP (Kapli et al. 2017).
#
#   --multi  = mPTP (each species has its own coalescent rate)
#   --single = PTP  (one rate for all species; the model behind bPTP,
#                    Zhang et al. 2013, which manuscript_v1 used)
#
# Steps:
#   1. --minbr_auto: estimate the minimum branch length from the alignment
#      (branches shorter than this are treated as zero, e.g. identical
#      sequences or sequencing noise). Falls back to 0.0001.
#   2a. "ml":   maximum-likelihood delimitation -> <prefix>.txt
#   2b. "mcmc": MCMC runs for support values  -> <prefix>.* files
#
# Usage:
#   run_mptp.sh ml   <rooted.nwk> <aln.fasta> <--multi|--single> <prefix>
#   run_mptp.sh mcmc <rooted.nwk> <aln.fasta> <--multi|--single> <prefix> <generations> <runs> <sample>
# =============================================================================

set -euo pipefail

mode=$1; tree=$2; aln=$3; model=$4; prefix=$5
mkdir -p "$(dirname "$prefix")"

# Empty tree (fewer than 4 sequences): nothing to delimit.
if [ "$(tr -d ' \n' < "$tree")" = "();" ]; then
    : > "${prefix}.txt"
    exit 0
fi

# 1. Minimum branch length from the alignment.
minbr_out=$(mptp --minbr_auto "$aln" --tree_file "$tree" --output_file "${prefix}.minbr" 2>&1 || true)
echo "$minbr_out"
# Take the last number printed on a line mentioning minbr.
minbr=$(echo "$minbr_out" | grep -i "minbr" | grep -Eo '[0-9]+\.?[0-9]*([eE]-?[0-9]+)?' | tail -n 1 || true)
if [ -z "$minbr" ]; then
    minbr=0.0001
    echo "Could not read --minbr_auto result; using $minbr"
fi
echo "minbr = $minbr"

if [ "$mode" = "ml" ]; then
    # 2a. Maximum-likelihood delimitation.
    mptp --ml "$model" --tree_file "$tree" --minbr "$minbr" --output_file "$prefix"
else
    # 2b. MCMC for support values.
    gens=$6; runs=$7; sample=$8
    mptp --mcmc "$gens" "$model" --mcmc_sample "$sample" --mcmc_runs "$runs" \
         --mcmc_log --tree_file "$tree" --minbr "$minbr" --output_file "$prefix"
fi
