# =============================================================================
# sanger.smk: the voucher-linked Sanger reference library (Sections 2.2, 2.6).
#
# Input: one curated FASTA per locus in data/sanger/{locus}.fasta. Chromatogram
# trimming and assembly in Geneious Prime is manual and happens before this.
#
# Steps per locus:
#   1. clean headers, add the "SAN-" platform tag, keep one sequence per specimen
#   2. translation check for COX1 (code 5) and H3 (code 1): report stop codons
#   3. MAFFT alignment (--auto --maxiterate 1000)
#   4. for coding loci: put the alignment in reading frame
#   5. Gblocks 0.91b with the less stringent settings
# Then all loci are concatenated with codon-position partitions and analysed
# with ModelFinder (MFP+MERGE) and IQ-TREE 2 (UFBoot -bnni, SH-aLRT).
# =============================================================================


rule sanger_clean:
    # Tidy headers (specimen ID only), prefix with "SAN-", keep the longest
    # sequence when a specimen has more than one for a locus.
    input:
        lambda w: os.path.join(config["sanger"]["dir"], f"{w.locus}.fasta")
    output:
        fasta=f"{OUT}/sanger/{{locus}}/{{locus}}.clean.fasta",
        metrics=f"{OUT}/sanger/{{locus}}/{{locus}}.clean.metrics.tsv"
    log:
        f"{OUT}/logs/sanger/{{locus}}.clean.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "python {SCRIPTS}/prefix_headers.py --input {input} --output {output.fasta} "
        "--prefix SAN --keep-longest --metrics {output.metrics} "
        "--metrics-label sanger_{wildcards.locus} > {log} 2>&1"


rule sanger_translation_check:
    # Translate COX1 and H3 and report sequences with internal stop codons in
    # every reading frame. Nothing is removed: review the report and fix the
    # sequences in Geneious if needed.
    input:
        f"{OUT}/sanger/{{locus}}/{{locus}}.clean.fasta"
    output:
        f"{OUT}/sanger/{{locus}}/{{locus}}.translation_report.tsv"
    params:
        table=lambda w: CODING[w.locus]
    wildcard_constraints:
        locus="|".join(CODING.keys())
    conda:
        "../envs/blast.yaml"
    shell:
        "python {SCRIPTS}/check_translation.py report --input {input} "
        "--table {params.table} --output {output}"


rule mafft_sanger:
    # Align each locus with MAFFT (--auto --maxiterate 1000), as in manuscript_v1.
    input:
        f"{OUT}/sanger/{{locus}}/{{locus}}.clean.fasta"
    output:
        f"{OUT}/sanger/{{locus}}/{{locus}}.aln.fasta"
    threads: 4
    log:
        f"{OUT}/logs/sanger/{{locus}}.mafft.log"
    conda:
        "../envs/align.yaml"
    shell:
        "mafft --auto --maxiterate 1000 --thread {threads} {input} > {output} 2> {log}"


rule sanger_frame:
    # For protein-coding loci, find the reading frame of the alignment (the
    # offset giving fewest stop codons) and trim it to whole codons, so that
    # codon-aware Gblocks and codon-position partitions are correct.
    input:
        f"{OUT}/sanger/{{locus}}/{{locus}}.aln.fasta"
    output:
        f"{OUT}/sanger/{{locus}}/{{locus}}.aln.inframe.fasta"
    params:
        table=lambda w: CODING[w.locus]
    wildcard_constraints:
        locus="|".join(CODING.keys())
    log:
        f"{OUT}/logs/sanger/{{locus}}.frame.log"
    conda:
        "../envs/align.yaml"
    shell:
        "python {SCRIPTS}/check_translation.py frame --input {input} "
        "--table {params.table} --output {output} > {log} 2>&1"


def gblocks_input(wildcards):
    """In-frame alignment for coding loci, plain alignment otherwise."""
    if wildcards.locus in CODING:
        return f"{OUT}/sanger/{wildcards.locus}/{wildcards.locus}.aln.inframe.fasta"
    return f"{OUT}/sanger/{wildcards.locus}/{wildcards.locus}.aln.fasta"


rule gblocks_sanger:
    # Remove ambiguously aligned positions with Gblocks (less stringent
    # settings). Coding loci use codon mode (-t=c) so whole codons are removed
    # and the reading frame is kept.
    input:
        gblocks_input
    output:
        f"{OUT}/sanger/{{locus}}/{{locus}}.gb.fasta"
    params:
        seqtype=lambda w: "c" if w.locus in CODING else "d"
    log:
        f"{OUT}/logs/sanger/{{locus}}.gblocks.log"
    conda:
        "../envs/align.yaml"
    shell:
        "bash {SCRIPTS}/run_gblocks.sh {input} {output} {params.seqtype} > {log} 2>&1"


rule sanger_concat:
    # Concatenate all loci by specimen. Missing loci are filled with "?".
    # Writes a NEXUS partition file with one partition per non-coding locus
    # and one per codon position for COX1 and H3.
    input:
        expand(f"{OUT}/sanger/{{locus}}/{{locus}}.gb.fasta", locus=SANGER_LOCI)
    output:
        fasta=f"{OUT}/sanger/concat/concat.fasta",
        nexus=f"{OUT}/sanger/concat/partitions.nex",
        coverage=f"{OUT}/sanger/concat/specimen_coverage.tsv",
        metrics=f"{OUT}/sanger/concat/concat.metrics.tsv"
    params:
        loci=" ".join(SANGER_LOCI),
        coding=" ".join(CODING.keys())
    log:
        f"{OUT}/logs/sanger/concat.log"
    conda:
        "../envs/align.yaml"
    shell:
        "python {SCRIPTS}/concat_alignments.py --alignments {input} "
        "--loci {params.loci} --coding {params.coding} "
        "--out-fasta {output.fasta} --out-nexus {output.nexus} "
        "--out-coverage {output.coverage} --metrics {output.metrics} > {log} 2>&1"


rule iqtree_sanger_concat:
    # Maximum-likelihood tree of the concatenated library, as in manuscript_v1:
    # ModelFinder with partition merging (-m MFP+MERGE), edge-proportional
    # partition rates (-p; "-spp" in IQ-TREE 1 syntax), 10,000 UFBoot replicates
    # with NNI optimisation (-bnni) and 10,000 SH-aLRT replicates.
    input:
        fasta=f"{OUT}/sanger/concat/concat.fasta",
        nexus=f"{OUT}/sanger/concat/partitions.nex"
    output:
        tree=f"{OUT}/sanger/concat/iqtree/concat.treefile",
        report=f"{OUT}/sanger/concat/iqtree/concat.iqtree"
    params:
        prefix=f"{OUT}/sanger/concat/iqtree/concat",
        ufboot=config["phylogeny"]["ufboot"],
        alrt=config["phylogeny"]["alrt"]
    threads: 16
    log:
        f"{OUT}/logs/sanger/concat.iqtree.log"
    conda:
        "../envs/align.yaml"
    shell:
        "iqtree2 -s {input.fasta} -p {input.nexus} -m MFP+MERGE "
        "-B {params.ufboot} -alrt {params.alrt} -bnni "
        "-T {threads} --prefix {params.prefix} -redo > {log} 2>&1"
