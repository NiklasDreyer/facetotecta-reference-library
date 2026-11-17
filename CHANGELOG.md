# Changelog

All notable changes to this project are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] - 2026-11-17

### Added
- Complete Snakemake workflow, rebuilt from the raw Illumina and PacBio data with the methods of manuscript_v1:
  - Sanger reference library, concatenation and IQ-TREE analysis
  - Illumina processing with Trimmomatic, BBDuk, DADA2, BLAST screen and TransDecoder
  - PacBio processing with read cleaning, mapPacBio.sh mapping and BLAST extraction
  - iterative phylogenetic authentication filter
  - species delimitation with VSEARCH, PTP, mPTP and ASAP
- New analyses for the *Molecular Ecology Resources* submission:
  - metaMATE NUMT filtering
  - PacBio COX1 BLAST rescue
  - threshold sweeps
  - congruence and MOTU grading
  - iNEXT richness estimation
  - platform-overlap plots
  - in-silico primer coverage
  - leave-one-MOTU-out assignment benchmark
  - phylogenetic placement
- Separate workflow for mining public metabarcoding datasets, with gain tables and Darwin Core export.
- `results/summary_numbers.tsv`, which collects every number for the manuscript.
- SLURM profile for Snakemake 8 and a one-command installer.
- Documentation, tests and GitHub Actions checks.
