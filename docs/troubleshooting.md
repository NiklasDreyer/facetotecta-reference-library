# Troubleshooting

## Steps to check on the first run

These tools change their command-line flags or output formats between versions. The pipeline handles the common cases, but check these on the first run:

| Step | What to check | Where |
|---|---|---|
| metaMATE | Flags (`metamate find -h`) and read-map orientation | `workflow/scripts/run_metamate.sh` |
| mPTP | How `--minbr_auto` reports the value (falls back to 0.0001; the log says so) | `workflow/scripts/run_mptp.sh` |
| ASAP | Output file names (`partition_glob`) and flags | `config/config.yaml`, `workflow/scripts/run_asap.py` |
| TransDecoder | Location of `longest_orfs.pep` | `rules/illumina.smk` |
| gappa | Column names of `per_query.tsv` | `workflow/scripts/run_placement.py` |

## Common problems

| Symptom | Fix |
|---|---|
| `MissingInputException` in a dry run | An input file is missing; the message names it. See [inputs.md](inputs.md). |
| Jobs fail immediately on SLURM | Check `slurm_account` and `slurm_partition`; look in `.snakemake/slurm_logs/`. |
| Conda environment build fails | Run `bash install.sh` again on the login node; check `conda config --show channel_priority` (should be `strict`). |
| `Gblocks produced no output` | The alignment is too short or too gappy; check the MAFFT output. |
| IQ-TREE runs out of memory | Increase `mem_mb` for the rule in `profiles/slurm/config.yaml`. |
| No reads left after BBDuk | `trimq=33` is too strict for the data; lower it and report the change. |
| DADA2 merges almost no reads | The amplicon is too long for the read length; set `merge_pairs: false`. |
| Empty `subset_*.fasta` | No authentic sequences from that platform; the delimitation steps then write empty results. |

## Getting help

Open an [issue](https://github.com/NiklasDreyer/facetotecta-reference-library/issues) with the failing rule, its log from `results/logs/`, and the output of `snakemake --version`.
