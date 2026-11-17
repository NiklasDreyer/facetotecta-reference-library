# Outputs

All results are written to `results/` (set by `outdir` in the config).

## The summary table

**`results/summary_numbers.tsv`** collects every `*.metrics.tsv` file written by the steps (and the mapping statistics) into one table:

| key | value | source |
|---|---|---|
| `sanger_COX1_sequences` | … | `sanger/COX1/COX1.clean.metrics.tsv` |
| `delim_sanger_concat_vsearch_motus` | … | `delim/sanger_concat/motu_counts.metrics.tsv` |

Fill every number in the manuscript from this file.

## Files by figure and table

| Manuscript item | File(s) |
|---|---|
| Table 1 (sequences and MOTUs per locus) | `delim/*/motu_counts.metrics.tsv`, `sanger/*/*.clean.metrics.tsv`, `sanger/concat/concat.metrics.tsv` |
| Table S1 (specimens) | `config/specimens.tsv`, `sanger/concat/specimen_coverage.tsv` |
| Table S3 (sequencing statistics) | `illumina/stats/*.seqkit.tsv`, `illumina/*/dada2/track_reads.tsv`, `pacbio/stats/cleaning.seqkit.tsv` |
| Table S5 (public datasets) | `mining/summary/gain_by_dataset.tsv` |
| Table S7 (PTP vs. mPTP) | `delim/*/membership/{ptp,mptp}.tsv` |
| Fig. 3 (library tree, distances, MOTU counts) | `sanger/concat/iqtree/concat.treefile`, `delim/congruence/` |
| Fig. 4 (Sesoko lineages across platforms) | `authenticate/COX1/final.treefile`, `richness/overlap_COX1_upset.pdf` |
| Fig. 5 (rarefaction) | `richness/inext_curves.pdf`, `richness/inext_estimates.tsv` |
| Fig. 6a (benchmark) | `benchmark/{marker}/{region}/roc.pdf`, `thresholds.tsv` |
| Fig. 6b (primers) | `primers/primer_mismatch_summary.pdf` |
| Fig. 7 (public detections) | `mining/summary/detections.tsv`, `detections_map.pdf` |
| Fig. S8 (threshold sweep) | `delim/*/sweep/threshold_sweep.pdf` |

## Folders

| Folder | Content |
|---|---|
| `db/` | Reference FASTA, taxonomy tables, BLAST databases |
| `sanger/{locus}/` | Cleaned sequences, translation reports, alignments, Gblocks output |
| `sanger/concat/` | Concatenated matrix, partitions, coverage, IQ-TREE results |
| `illumina/qc/` | FastQC and MultiQC reports |
| `illumina/{marker}/` | DADA2 ASVs and tables, screen report, TransDecoder and metaMATE results, `final.fasta` |
| `pacbio/qc/`, `pacbio/clean/` | Read QC and each cleaning step |
| `pacbio/{locus}/` | Extracted COX1, 18S and 16S reads |
| `authenticate/{locus}/` | `authentic.fasta`, `removed.tsv` (with reasons), final tree and alignment, per-platform subsets |
| `delim/{dataset}/` | Alignment, tree, membership tables per method, sweep |
| `delim/congruence/` | Pairwise agreement, MOTU grades |
| `richness/` | iNEXT tables and plots, platform overlap |
| `primers/` | Mismatch tables and amplicons per primer pair |
| `benchmark/` | Leave-one-out hits, ROC tables and plots, chosen thresholds |
| `placement/` | Reference trees, EPA-ng results, assignments |
| `mining/` | Raw and trimmed public reads, ASVs, assignments, authentication, summary, Darwin Core |
| `logs/` | One log per job |

## Membership tables

Every delimitation method writes `results/delim/{dataset}/membership/{method}.tsv`:

```
sequence    motu
SAN-Y0123   SAN-Y0123
SAN-Y0456   SAN-Y0123
```
