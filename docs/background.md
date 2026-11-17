# Background

## Dark taxa and reference databases

Most marine invertebrate species are still undescribed. A large part of this gap lies in **dark taxa**: groups that are diverse, hard to collect or identify, and taxonomically neglected. DNA metabarcoding of environmental samples (eDNA) and bulk plankton can now survey whole communities across ocean basins, but every such survey depends on a **reference database**. A read can only be named if it matches a correctly identified reference sequence. Lineages without references are reported as "unassigned" or are assigned to the wrong taxon, so they disappear from the results even when they are abundant.

## Facetotecta: the y-larva problem

Facetotecta (Pancrustacea: Thecostraca) are known almost only from their planktonic larvae, the **y-larvae**, first described by Hansen in 1899. They are found:

- from the Arctic to tropical coral reefs, and at depths from the intertidal to more than 8,700 m
- in high local diversity at some sites, especially Sesoko Island, Okinawa

Yet the **adult stage has never been found**. The life cycle has 5–7 naupliar instars (y-nauplii), a y-cyprid that probably attaches to a host, and an ypsigon stage that may invade it. This closely resembles the parasitic rhizocephalan barnacles, which control the physiology and behaviour of their hosts. Facetotectans are therefore thought to be **endoparasites of still-unknown hosts**, and a diverse group of them could matter ecologically.

Several facts keep the group dark:

- **16 species are formally described**, only two of them with DNA data
- wild-caught larvae **cannot be identified reliably to species**, so every reference sequence must come from an individually reared, imaged and vouchered larva
- before this work, GenBank held very few facetotectan sequences, and **no metabarcoding or eDNA study had reported Facetotecta**

Published estimates of diversity are hypotheses that this study tests:

| Hypothesis | Source |
|---|---|
| At least 100 species worldwide | Hansen 1899 |
| At least 40 morphotypes at Sesoko Island | Glenner et al. 2008 |
| The four species described from Sesoko represent < 10% of its diversity | Olesen et al. 2022 |

## Study design

### 1. A voucher-linked reference library

Between 2017 and 2022, plankton was collected with hand nets at Sesoko Island (Japan), around Keelung and on Green Island (Taiwan), in the Azores, and at the White Sea Biological Station (Russia), with further specimens from Malaysia and Antarctica. Each y-larva was:

1. reared individually and followed through development with video microscopy
2. sorted into morphotypes from the images
3. extracted for DNA, with its moulted exoskeleton (exuviae) mounted on a slide as the morphological voucher
4. Sanger-sequenced for up to six loci: COX1, 12S, 16S, 18S, 28S and histone H3

Each sequence is therefore tied to a physical voucher and a developmental image series.

### 2. Three sequencing approaches at one site

At Sesoko Island, two bulk samples of hand-sorted y-larvae (Y1: October 2018, about 600 larvae; Y2: June 2019, about 3,000 larvae) were used for:

- **Illumina MiSeq amplicons** of COX1 and 18S (genomic DNA). These are cheap and deep, but sensitive to NUMTs, PCR errors and chimeras.
- **PacBio Iso-Seq** of full-length cDNA from the pooled RNA. These are long, accurate reads derived from transcripts, so they largely avoid nuclear mitochondrial pseudogenes (NUMTs).

Combined with the Sanger library, a lineage supported by both genomic and transcript data is unlikely to be an artefact.

### 3. From sequences to lineages

All high-throughput sequences go through an iterative phylogenetic filter. Sequences are removed when they fall outside a supported Facetotecta clade, sit on suspiciously long branches, or form long-branch clades seen by only one platform. The authentic sequences are then delimited into molecular operational taxonomic units (MOTUs) with distance methods (VSEARCH, ASAP) and tree-based methods (PTP, mPTP). The agreement between methods and between mitochondrial and nuclear loci is measured, and each MOTU is graded by how well it is supported.

### 4. Testing the library as a resource

A reference library is only useful if it assigns reads correctly. The workflow therefore:

- runs **leave-one-MOTU-out cross-validation** against other crustaceans (barnacles, copepods, ostracods and others) to measure true- and false-positive rates at every identity threshold
- tests whether **standard metabarcoding primers** (Leray COX1, 18S V4, 18S V9) match facetotectan sequences
- places non-target sequences on the reference tree to check that they are not mistaken for Facetotecta

### 5. Finding Facetotecta in public data

Finally, public metabarcoding datasets are re-analysed twice: once with the public reference alone, and once with the new library added. The difference is the measurable gain from the library: sequences that move from "unassigned" or "another taxon" to Facetotecta. These are confirmed by phylogenetic placement and the authentication filter, then exported as Darwin Core records for GBIF/OBIS.

## Why this matters

Making Facetotecta detectable turns a group known only from scattered larvae into one that existing and future eDNA surveys can record. That opens the way to mapping where facetotectan lineages occur, testing whether they are mostly local or widely dispersed, and finding their hosts, and with them the adult y-organism. The same workflow applies to other dark taxa that can only be identified after rearing.

## Key references (short form; full citations are in the manuscript)

| Reference | Relevance |
|---|---|
| Hansen 1899 | First description of y-larvae; hypothesis of at least 100 species |
| Glenner et al. 2008 | Y-larva diversity at Sesoko Island; at least 40 morphotypes |
| Pérez-Losada et al. 2009 | Position of Facetotecta within Thecostraca |
| Olesen et al. 2022 | Integrative taxonomy of y-larvae from Sesoko; single-specimen rearing workflow |
| Dreyer et al. 2023a, b; 2024; 2025 | Global distribution review, single-specimen sequencing pipeline, phylogenomics |
| Andújar et al. 2021 | metaMATE, NUMT and error filtering of metabarcoding data |
| Kapli et al. 2017 | mPTP species delimitation |
| Puillandre et al. 2021 | ASAP species delimitation |
| Hsieh, Ma & Chao 2016 | iNEXT rarefaction and extrapolation |
