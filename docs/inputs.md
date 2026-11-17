# Inputs

| Input | Location | Format |
|---|---|---|
| Curated Sanger sequences | `data/sanger/{locus}.fasta` | One FASTA per locus: COX1, 12S, 16S, 18S, 28S, H3 |
| Raw Illumina reads | paths in `config/illumina_samples.tsv` | Paired FASTQ (gzip) |
| PacBio CCS reads | `pacbio.ccs_fastq` in the config | FASTQ (gzip) |
| PacBio CCS BAM (optional) | `pacbio.ccs_bam` | BAM with `np` tags, used for the ≥10-pass filter |
| PacBio subreads (optional) | `pacbio.deepconsensus.subreads_bam` | Only for DeepConsensus |
| Specimen metadata | `config/specimens.tsv` | Table S1 |
| Outgroups | `resources/outgroups/{COX1,18S,16S}.fasta` | FASTA, headers starting `OUT-` |
| Type AE* COX1 | `resources/reference/typeAE_COX1.fasta` | Complete COX1, 1,527 bp |
| PacBio oligos | `resources/pacbio/primers_tso.fasta` | SMARTer primers, oligo-dT, TSO |
| MIDORI2 CO1 and lrRNA | `resources/db/` | "LONGEST" nucleotide FASTA |
| SILVA SSU Ref NR99 | `resources/db/` | `tax_silva` FASTA, unzipped |
| Primer pairs | `config/insilico_primers.tsv` | name, marker, forward, reverse, reference |
| Public datasets | `config/public_datasets.tsv` | For the mining workflow only |

## Sanger FASTA files

Chromatogram trimming, assembly and proofreading in Geneious Prime stay manual. Export one FASTA per locus. Each header must **start with the specimen ID**; anything after the first space or `|` is ignored:

```
>Y0123 Sesoko 2018-10-02
ACTTTATATTTTATTTTTGG...
```

The pipeline adds the platform tag (`SAN-Y0123`). If a specimen has more than one sequence for a locus, it keeps the longest.

## `config/illumina_samples.tsv`

| sample | marker | r1 | r2 |
|---|---|---|---|
| Y1 | COX1 | data/raw/illumina/Y1_COX1_R1.fastq.gz | …_R2.fastq.gz |

## `config/specimens.tsv`

Columns: `specimen_id`, `site`, `region`, `sampling_event`, `date`, `lat`, `lon`, `stage`, `voucher`, `new`.

`sampling_event` is one plankton haul or collection day. The richness estimates use it as the sampling unit.

## `config/public_datasets.tsv`

Columns: `dataset`, `marker`, `runs` (SRA accessions separated by `;`), `layout` (PE/SE), `fwd`, `rev`, `region`, `lat`, `lon`, `citation`.

## Reference databases

| Database | Where to get it |
|---|---|
| MIDORI2 | www.reference-midori.info: "LONGEST" nucleotide files for CO1 and lrRNA |
| SILVA | www.arb-silva.de: `SILVA_<release>_SSURef_NR99_tax_silva.fasta.gz` |

Set the paths under `references:` in the config.
