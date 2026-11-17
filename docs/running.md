# Running the workflow

Activate the environment first:

```bash
conda activate yfacet-snakemake
```

## Main workflow

```bash
snakemake --profile profiles/slurm -n    # dry run: lists every job
snakemake --profile profiles/slurm       # run
```

Snakemake must keep running while it submits jobs. Start it inside `screen` or `tmux`, or submit it as a small long-running SLURM job:

```bash
sbatch --time=7-00:00:00 --mem=4G --cpus-per-task=1 --wrap \
  "source \$(conda info --base)/etc/profile.d/conda.sh && conda activate yfacet-snakemake && snakemake --profile profiles/slurm"
```

## Public-data mining (after the main workflow)

```bash
snakemake -s workflow/Snakefile_mining --profile profiles/slurm -n
snakemake -s workflow/Snakefile_mining --profile profiles/slurm
```

Downloads run on the login node; everything else runs on SLURM.

## Partial runs

| Goal | Command |
|---|---|
| Sanger library only | `snakemake --profile profiles/slurm results/sanger/concat/iqtree/concat.treefile` |
| Check the first delimitation | `snakemake --profile profiles/slurm results/delim/sanger_COX1/motu_counts.metrics.tsv` |
| Illumina up to ASVs | `snakemake --profile profiles/slurm --until dada2` |
| PacBio cleaning only | `snakemake --profile profiles/slurm results/pacbio/stats/cleaning.metrics.tsv` |
| Summary table only | `snakemake --profile profiles/slurm results/summary_numbers.tsv` |

## Without SLURM (workstation)

```bash
snakemake --sdm conda --cores 16
```

## Reruns and inspection

| Goal | Command |
|---|---|
| See why jobs would rerun | `snakemake --profile profiles/slurm -n --reason` |
| Rerun one rule | `snakemake --profile profiles/slurm -R delim_vsearch` |
| Draw the workflow graph | `snakemake --rulegraph \| dot -Tpdf > rulegraph.pdf` |
| HTML report of runtimes | `snakemake --report report.html` |

Logs for every job are in `results/logs/`, and SLURM logs in `.snakemake/slurm_logs/`.

## Suggested first run

1. Run the Sanger part alone (`results/delim/sanger_COX1/motu_counts.metrics.tsv`) and compare the COX1 MOTU counts with manuscript_v1.
2. Run up to DADA2 and check the read tracking in `results/illumina/*/dada2/track_reads.tsv`.
3. Run PacBio cleaning and check `results/pacbio/stats/cleaning.metrics.tsv`.
4. Then run everything.
