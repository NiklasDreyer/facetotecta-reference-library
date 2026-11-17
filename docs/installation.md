# Installation

All software is installed with conda. Run the installer once, on a **login node** (compute nodes often have no internet access).

## Requirements

- Linux with a SLURM cluster (the pipeline also runs on a single workstation with `--cores N`)
- conda or mamba. If neither is available, install [Miniforge](https://github.com/conda-forge/miniforge):
  ```bash
  curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
  source $HOME/miniforge3/etc/profile.d/conda.sh
  ```
  On many clusters this is available as a module (for example `module load miniforge3`).
- About 30 GB of disk space for the conda environments and reference databases

## One-command installation

```bash
git clone git@github.com:NiklasDreyer/facetotecta-reference-library.git
cd facetotecta-reference-library
bash install.sh
```

`install.sh` does the following:

1. Creates the `yfacet-snakemake` environment with Snakemake ≥ 8 and `snakemake-executor-plugin-slurm`.
2. Builds every per-rule environment in `workflow/envs/` in advance, so jobs never download anything.
3. Copies the default metaMATE specifications file to `resources/metamate/`.
4. Writes the exact package versions of every environment to `envs_lock/`, for the Methods section and for reproducibility. Commit these files.

## Software installed

| Environment | Tools |
|---|---|
| `qc` | FastQC, MultiQC, NanoPlot, NanoQC, seqkit |
| `trim` | Trimmomatic, BBTools, seqkit, pigz |
| `dada2` | R, DADA2, optparse |
| `blast` | BLAST+, Biopython, pandas |
| `align` | MAFFT, Gblocks, IQ-TREE 2, seqkit |
| `authenticate` | MAFFT, Gblocks, IQ-TREE 2, VSEARCH, ete3 |
| `vsearch` | VSEARCH |
| `transdecoder` | TransDecoder |
| `metamate` | metaMATE (installed from GitHub), MAFFT, FastTree, R |
| `mapping` | BBTools (mapPacBio.sh), minimap2, samtools, pysam |
| `mptp` | mPTP |
| `inext` | R, iNEXT, ggplot2 |
| `placement` | MAFFT, RAxML-NG, EPA-ng, gappa, IQ-TREE 2 |
| `benchmark` | BLAST+, VSEARCH, numpy, pandas, matplotlib |
| `pytools` | pandas, scikit-learn, matplotlib, UpSetPlot, ete3 |
| `sra`, `cutadapt` | SRA Toolkit, cutadapt (mining workflow) |
| `deepconsensus` | DeepConsensus, pbccs, actc (only when enabled) |

## ASAP (manual step)

ASAP is not distributed through conda:

1. Download the command-line version from the [ASAP website](https://bioinfo.mnhn.fr/abi/public/asap/).
2. Compile it and copy the binary to `bin/asap` in this repository.
3. Run `bin/asap -h` and check its flags against `delimitation.asap.command` in `config/config.yaml`.
4. Set `delimitation.asap.enabled: true`.

Until then, the other delimitation methods run without it. ASAP can also be run on its web server, using the alignments in `results/delim/*/aln.fasta`.

## Updating

```bash
git pull
bash install.sh      # rebuilds only the environments that changed
```
