# Changes from manuscript_v1

These differences between the pipeline and the methods text of manuscript_v1 should be stated in the revised Methods.

## Corrections

| Step | manuscript_v1 | Pipeline | Reason |
|---|---|---|---|
| mapPacBio.sh | `qin=42`, `k=script` | `qin=33`, `k` default | `qin` is the quality offset (33 or 64); `k=script` is a typo |
| BBDuk PacBio trimming | `ktrim=r ktrim=l` in one call | Two passes (right, then left) | BBDuk accepts one `ktrim` value per call |
| Gblocks for COX1 and H3 | DNA mode | Codon mode after putting the alignment in frame | Keeps codon-position partitions correct |
| IQ-TREE partitions | `-spp` | `-p` | Same model (edge-proportional); IQ-TREE 2 syntax |
| PacBio orientation | FASTX reverse-complement | Per-read orientation by poly-A/poly-T | Same aim, done per read |
| bPTP | bPTP web/Python | PTP (`mptp --single`) and mPTP (`--multi`) with MCMC support | PTP is the bPTP model; mPTP corrects oversplitting |

## Additions

| Addition | Purpose |
|---|---|
| metaMATE on Illumina COX1 ASVs | NUMT and error filtering beyond stop codons |
| COX1 BLAST rescue of PacBio reads | Recovers lineages too divergent to map to one reference |
| Per-platform subsets and UpSet plots | Quantifies agreement between genomic and transcript data |
| Threshold sweep | Shows how MOTU counts depend on the threshold |
| Congruence and MOTU grades | Agreement between methods and loci |
| iNEXT rarefaction, Chao1/Chao2, coverage | Turns observed counts into richness estimates |
| Leave-one-MOTU-out benchmark | Measures false- and true-positive rates of assignment |
| In-silico primer tests | Explains past non-detection; guides primer choice |
| Two-pass public-data assignment | Measures what the library adds |
| Darwin Core export | Publishes detections to GBIF/OBIS |

## Settings left to decide

See the **CHECK** items in [configuration.md](configuration.md).
