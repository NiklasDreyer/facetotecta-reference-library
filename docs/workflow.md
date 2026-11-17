# Workflow in detail

Each module is a rule file in `workflow/rules/`. Each step calls a commented script in `workflow/scripts/` that can also be run on its own (`-h` shows its options).

## Platform tags

Every sequence ID starts with a tag that records where it came from:

| Tag | Meaning |
|---|---|
| `SAN-` | Sanger, single vouchered specimen |
| `ILL-` | Illumina ASV |
| `PB-` | PacBio CCS read |
| `PUB-` | Public dataset (mining workflow) |
| `OUT-` | Outgroup |

## 1. Sanger reference library (`sanger.smk`)

1. `prefix_headers.py` adds the `SAN-` tag, keeps one sequence per specimen and writes the metrics.
2. `check_translation.py report` translates COX1 (code 5) and H3 (code 1) and reports stop codons.
3. MAFFT aligns each locus (`--auto --maxiterate 1000`).
4. `check_translation.py frame` puts coding alignments in frame.
5. Gblocks trims with the less stringent settings (codon mode for COX1 and H3), via `run_gblocks.sh`.
6. `concat_alignments.py` concatenates by specimen and writes a partition file with codon positions.
7. IQ-TREE 2 builds the tree: `-m MFP+MERGE -p`, 10,000 UFBoot replicates with `-bnni`, and 10,000 SH-aLRT replicates.

## 2. Reference databases (`references.smk`)

- **public:** MIDORI2 (COX1, lrRNA) or SILVA (eukaryotic 18S), reformatted with short IDs and a taxonomy table.
- **augmented:** public plus the Sanger library, labelled Facetotecta.

## 3. Illumina amplicons (`illumina.smk`)

1. FastQC on the raw reads.
2. Trimmomatic: `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:3`.
3. BBDuk: `qtrim=rl trimq=33 minlen=100 ktrim=r k=23 mink=11 hdist=1 tpe tbo`.
4. FastQC on the trimmed reads, then MultiQC.
5. DADA2 (`dada2_asv.R`): filtering, error learning, denoising, merging and chimera removal.
6. BLAST screen against the augmented database: ASVs whose top hit is Facetotecta or Arthropoda are kept.
7. COX1 only: TransDecoder translation screen (invertebrate mitochondrial code), then metaMATE NUMT and error filtering (NEW).

## 4. PacBio Iso-Seq (`pacbio.smk`)

1. NanoPlot, NanoQC and BBDuk histograms.
2. ≥ 10 passes (when the CCS BAM is available).
3. BBDuk removes the SMARTer primers, oligo-dT and TSO (right end, then left end); mean quality must be ≥ Q40.
4. Orientation: reads with a 5′ poly-T are reverse-complemented.
5. Entropy filter (`entropy=0.5 entropywindow=50 entropyk=5`).
6. Optional DeepConsensus (needs subreads).
7. **COX1:** mapPacBio.sh to the Type AE* COX1 (minimap2 for comparison). Primary alignments are trimmed to the reference, and reads with internal stops (frameshifts) are removed. A NEW BLAST rescue recovers divergent reads.
8. **18S and 16S:** BLAST against the augmented database; reads whose top hit is Facetotecta (e-value ≤ 1e-50) are kept, as the aligned region.

## 5. Authentication (`authenticate.smk`)

1. Each platform is pre-clustered separately with VSEARCH.
2. The iterative filter (`authenticate_filter.py`) aligns the sequences with outgroups, trims them with Gblocks, builds an IQ-TREE tree and roots it on the outgroups. It then removes high-throughput sequences that:
   - fall outside the Sanger-defined Facetotecta clade
   - sit on a terminal branch more than 10 times longer than their sister's
   - form long-branch clades from one platform with no Sanger sequence
   - (optional) have no other platform in their closest supported clade

   The rounds repeat until nothing is removed. Sanger sequences are never removed.
3. The authentic set is split by platform and sample (`make_subsets.py`).

## 6. Species delimitation (`delimitation.smk`)

The delimitation runs on every Sanger locus, on the concatenated matrix, on all authentic HTS sequences per locus, and on each platform's subset:

- **VSEARCH** `--cluster_fast` at the locus threshold.
- **PTP** (`mptp --single`) and **mPTP** (`--multi`): maximum-likelihood partitions plus MCMC support, using `--minbr_auto`.
- **ASAP** (optional).
- **Threshold sweep** from 90 to 99% (NEW).
- **Congruence** (NEW): adjusted Rand index and match ratio between methods and between COX1 and 18S. MOTUs are graded A, B or C.

## 7. Richness (`richness.smk`, NEW)

- iNEXT abundance-based curves (specimens), with Chao1 and sample coverage.
- Incidence-based curves (sampling events), with Chao2.
- UpSet plots of MOTUs shared between Sanger, PacBio, Illumina Y1 and Illumina Y2.

## 8. Primers (`primers.smk`, NEW)

`primer_scan.py` finds the best binding site of each primer and counts mismatches, including 3′-end mismatches, with IUPAC codes handled. It also extracts amplicons for the benchmark.

## 9. Benchmark (`benchmark.smk`, `placement.smk`, NEW)

- Non-target crustaceans are sampled from the public reference.
- Leave-one-MOTU-out cross-validation with BLAST (top hit and identity) and SINTAX (confidence) gives ROC curves and the threshold meeting a 1% false-positive rate. The manuscript_v1 cut-offs (90%, 92%) are also scored.
- Repeated on the regions amplified by each primer pair.
- EPA-ng placement of non-targets on the Facetotecta + outgroup tree.

## 10. Public-data mining (`Snakefile_mining`)

1. SRA download, then cutadapt (study primers, untrimmed reads discarded).
2. DADA2.
3. Two-pass BLAST: public reference vs. public plus library. Candidates are top-hit Facetotecta above the identity and e-value cut-offs.
4. Pooled per marker: EPA-ng placement (LWR ≥ 0.9), the authentication filter together with the library, and VSEARCH MOTUs.
5. Gain table per dataset (was unassigned, other taxon, or already Facetotecta), a map, and Darwin Core export.
