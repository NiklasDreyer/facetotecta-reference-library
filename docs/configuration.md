# Configuration

All settings are in [`config/config.yaml`](../config/config.yaml), with comments for each value. Values reproduce manuscript_v1 unless marked **NEW** (an added analysis) or **CHECK** (manuscript_v1 contradicts itself or leaves the value out).

## Values to set before the first run

| Setting | Why |
|---|---|
| `delimitation.vsearch_id.COX1` | **CHECK.** The Methods say 95% and the Results report ~97% from ASAP. |
| `authenticate.precluster_id` | **CHECK.** manuscript_v1 used 0.97. 1.0 (dereplication only) keeps the within-species variation that ASAP and PTP need. |
| `illumina.trimmomatic_steps` | **CHECK.** `MINLEN:3` looks like a typo; BBDuk removes reads under 100 bp afterwards anyway. |
| `illumina.bbduk` | **CHECK.** `trimq=33` is very strict for MiSeq data; report how many reads remain. |
| `illumina.dada2.merge_pairs` | **CHECK.** 2 × 150 bp reads cannot merge amplicons longer than about 270 bp. |
| `illumina.metamate.expected_length` | **CHECK.** Set from your COX1 primers. |
| `pacbio.adapters` | **CHECK.** Use the oligos you actually used (SMARTer and/or NEBNext TSO). |
| `pacbio.deepconsensus.enabled` | **CHECK.** Needs subreads; leave it `false` with CCS reads only. |
| `authenticate.require_cross_platform` | **CHECK.** The intended meaning of criterion (4) in manuscript_v1 is unclear. |

## Main parameter groups

| Group | Contents |
|---|---|
| `sanger` | Loci; genetic codes for COX1 (5) and H3 (1) |
| `references` | MIDORI2 and SILVA paths |
| `illumina` | Trimmomatic, BBDuk and DADA2 settings; BLAST screen; TransDecoder; metaMATE |
| `pacbio` | Pass and quality filters, oligo trimming, entropy filter, mapPacBio.sh, minimap2, rRNA BLAST |
| `authenticate` | Support thresholds (SH-aLRT 80, UFBoot 95), long-branch factor (10×), replicates |
| `delimitation` | VSEARCH thresholds per locus, mPTP MCMC settings, ASAP command, sweep range, grading |
| `richness` | Dataset and method for iNEXT, bootstrap replicates |
| `benchmark` | Non-target groups, folds, threshold range, target false-positive rate |
| `primers` | Primer table, allowed mismatches |
| `mining` | Identity and e-value cut-offs, minimum placement LWR |

## SLURM profile

[`profiles/slurm/config.yaml`](../profiles/slurm/config.yaml) sets the account, partition, default resources and the larger resources for heavy rules (IQ-TREE, authentication, DADA2, mapping). Set at least:

```yaml
default-resources:
  slurm_account: "YOUR_ACCOUNT"
  slurm_partition: "YOUR_PARTITION"
```

If your cluster needs extra options (QoS, constraints), add `slurm_extra: "'--qos=normal'"` to `default-resources`.
