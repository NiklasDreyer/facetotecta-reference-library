#!/usr/bin/env bash
# =============================================================================
# run_toy_tests.sh: run the stand-alone scripts that need only Python
# (pandas, numpy, scikit-learn) on tiny toy data in tests/data/, and check
# their outputs. External tools (MAFFT, IQ-TREE, BLAST, ...) are not needed.
#
# Run from the repository root:  bash tests/run_toy_tests.sh
# =============================================================================

set -euo pipefail

S=workflow/scripts
D=tests/data
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

# Fail with a message unless two values are equal.
expect() { if [ "$1" != "$2" ]; then echo "FAIL: $3 (got '$1', expected '$2')"; exit 1; else echo "ok: $3"; fi; }

# Reference database formatting (MIDORI and SILVA headers).
python $S/build_reference_db.py public --input $D/midori.fasta --format midori --out-fasta $T/pub.fa --out-tax $T/pub.tax 2>/dev/null
expect "$(grep -c '^pub' $T/pub.tax)" 2 "MIDORI reformatting"
python $S/build_reference_db.py public --input $D/silva.fasta --format silva --out-fasta $T/s.fa --out-tax $T/s.tax 2>/dev/null
expect "$(grep -c '^pub' $T/s.tax)" 1 "SILVA keeps eukaryotes only"

# Header cleaning, one sequence per specimen.
python $S/prefix_headers.py --input $D/sanger_cox1.fasta --output $T/c.fa --prefix SAN --keep-longest 2>/dev/null
expect "$(grep -c '^>SAN-' $T/c.fa)" 2 "header cleaning and duplicate removal"

# Concatenation with codon partitions.
python $S/concat_alignments.py --alignments $D/aln_cox1.fasta $D/aln_18s.fasta --loci COX1 18S --coding COX1 \
    --out-fasta $T/cat.fa --out-nexus $T/p.nex --out-coverage $T/cov.tsv
expect "$(grep -c charset $T/p.nex)" 4 "partition file (3 codon positions + 18S)"

# Primer scanning.
python $S/primer_scan.py --input $D/primer_targets.fasta --forward GGTACC --reverse GAATTC \
    --max-mismatch 1 --three-prime 2 --out-table $T/pt.tsv --out-amplicons $T/amp.fa
expect "$(grep -c amplifiable $T/pt.tsv)" 2 "primer scan finds two amplifiable sequences"

# Membership parsing, counting and congruence.
mkdir -p $T/D/sanger_COX1/membership $T/D/sanger_18S/membership
python $S/parse_uc.py --uc $D/clusters.uc --output $T/D/sanger_COX1/membership/vsearch.tsv
python $S/parse_mptp.py --input $D/mptp.txt --output $T/D/sanger_COX1/membership/mptp.tsv --motu-prefix mptp
cp $D/ptp.tsv $T/D/sanger_COX1/membership/ptp.tsv
for m in vsearch mptp ptp; do cp $D/nuclear.tsv $T/D/sanger_18S/membership/$m.tsv; done
python $S/delim_congruence.py --delim-dir $T/D --datasets sanger_COX1 sanger_18S --methods vsearch mptp ptp \
    --primary sanger_COX1 --nuclear sanger_18S --out-pairs $T/pairs.tsv --out-grades $T/grades.tsv
expect "$(awk -F'\t' '$1=="SAN-C"{print $6}' $T/grades.tsv)" A "MOTU grading"

# Subsets per platform.
python $S/make_subsets.py --fasta $D/authentic.fasta --members $D/precluster_members.tsv \
    --asv-table $D/asv_table.tsv --expected PB ILL_Y1 ILL_Y2 --outdir $T --platforms-out $T/plat.tsv >/dev/null
expect "$(grep -c '>' $T/subset_ILL_Y2.fasta)" 2 "Illumina sample subsets"

# Summary table.
printf "k1\t1\nk2\t2\n" > $T/x.metrics.tsv
python $S/summarize_numbers.py --results $T --output $T/sum.tsv >/dev/null
expect "$(grep -c '^k[0-9]' $T/sum.tsv)" 2 "summary table"

echo "All toy tests passed."
