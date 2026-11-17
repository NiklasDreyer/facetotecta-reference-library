# =============================================================================
# pacbio.smk: PacBio Iso-Seq CCS reads from the pooled Y1+Y2 RNA (Sections 2.5, 2.8).
#
#   QC (NanoPlot, NanoQC, BBDuk histograms)
#   -> optional removal of reads with < 10 passes (needs the CCS BAM)
#   -> BBDuk: remove primers/oligo-dT/TSO (right, then left), mean quality >= Q40
#   -> orient reads so poly(A) is at the 3' end (reverse-complement poly(T) reads)
#   -> BBDuk entropy filter
#   -> optional DeepConsensus (needs subreads)
#   -> rename reads PB-000001 ...
#   COX1: mapPacBio.sh to Type AE* COX1 (+ minimap2 for comparison), trim reads
#         to the reference, discard reads with frameshifts/stops,
#         NEW: rescue divergent COX1 reads by BLAST against the Sanger library
#   18S/16S: BLAST against SILVA/MIDORI + Sanger, keep reads whose top hit is
#            Facetotecta (e-value <= 1e-50)
# =============================================================================

PB = f"{OUT}/pacbio"


rule pacbio_qc_raw:
    # Read length / quality plots (NanoPlot, NanoQC) and counts of the raw CCS reads.
    input:
        config["pacbio"]["ccs_fastq"]
    output:
        nanoplot=directory(f"{PB}/qc/nanoplot"),
        nanoqc=directory(f"{PB}/qc/nanoqc"),
        stats=f"{PB}/qc/raw.seqkit.tsv"
    threads: 4
    log:
        f"{OUT}/logs/pacbio/qc_raw.log"
    conda:
        "../envs/qc.yaml"
    shell:
        "NanoPlot --fastq {input} --N50 --plots dot kde -t {threads} -o {output.nanoplot} > {log} 2>&1; "
        "nanoQC {input} -o {output.nanoqc} >> {log} 2>&1; "
        "seqkit stats -a -T {input} > {output.stats}"


rule pacbio_bbduk_hist:
    # BBDuk histograms (base, quality, GC, average quality, length), as in
    # manuscript_v1. Kept separate because it needs the trim environment.
    input:
        config["pacbio"]["ccs_fastq"]
    output:
        touch(f"{PB}/qc/bbduk_hist.done")
    params:
        d=f"{PB}/qc/bbduk_hist"
    log:
        f"{OUT}/logs/pacbio/bbduk_hist.log"
    conda:
        "../envs/trim.yaml"
    shell:
        "mkdir -p {params.d}; bbduk.sh in={input} bhist={params.d}/bhist.txt "
        "qhist={params.d}/qhist.txt gchist={params.d}/gchist.txt "
        "aqhist={params.d}/aqhist.txt lhist={params.d}/lhist.txt gcbins=auto > {log} 2>&1"


rule pacbio_passes:
    # Remove reads with fewer than min_passes full passes. The number of passes
    # ("np" tag) exists only in the CCS BAM; without it this step copies the
    # reads unchanged and says so in the log.
    input:
        fastq=config["pacbio"]["ccs_fastq"]
    output:
        f"{PB}/clean/01_passes.fastq.gz"
    params:
        bam=config["pacbio"]["ccs_bam"],
        minp=config["pacbio"]["min_passes"]
    log:
        f"{OUT}/logs/pacbio/passes.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "python {SCRIPTS}/filter_ccs_passes.py --fastq {input.fastq} "
        "--bam '{params.bam}' --min-passes {params.minp} --out {output} > {log} 2>&1"


rule pacbio_trim:
    # Remove SMARTer primers, oligo-dT and template-switching oligos from both
    # ends (ktrim=r, then ktrim=l, k=21 mink=11 hdist=2), then remove reads with
    # mean quality < Q40 (maq=40).
    input:
        reads=f"{PB}/clean/01_passes.fastq.gz",
        adapters=config["pacbio"]["adapters"]
    output:
        f"{PB}/clean/02_trimmed.fastq.gz"
    params:
        opts=config["pacbio"]["bbduk_trim"],
        maq=config["pacbio"]["min_mean_q"],
        tmp=f"{PB}/clean/02_tmp_right.fastq.gz"
    threads: 4
    log:
        f"{OUT}/logs/pacbio/trim.log"
    conda:
        "../envs/trim.yaml"
    shell:
        "bbduk.sh in={input.reads} out={params.tmp} ref={input.adapters} ktrim=r "
        "{params.opts} threads={threads} > {log} 2>&1; "
        "bbduk.sh in={params.tmp} out={output} ref={input.adapters} ktrim=l "
        "{params.opts} maq={params.maq} threads={threads} >> {log} 2>&1; "
        "rm -f {params.tmp}"


rule pacbio_orient:
    # Reads carrying a 5' poly(T) are reverse-complemented so that every read
    # is in transcript (sense) orientation with poly(A) at the 3' end.
    # (manuscript_v1 used FASTX; this script does the same per read.)
    input:
        f"{PB}/clean/02_trimmed.fastq.gz"
    output:
        reads=f"{PB}/clean/03_oriented.fastq.gz",
        report=f"{PB}/clean/03_orientation.metrics.tsv"
    log:
        f"{OUT}/logs/pacbio/orient.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "python {SCRIPTS}/orient_polya.py --input {input} --output {output.reads} "
        "--report {output.report} > {log} 2>&1"


rule pacbio_entropy:
    # Remove low-complexity reads (e.g. long homopolymers) with BBDuk's
    # entropy filter (entropy=0.5 entropywindow=50 entropyk=5).
    input:
        f"{PB}/clean/03_oriented.fastq.gz"
    output:
        f"{PB}/clean/04_entropy.fastq.gz"
    params:
        opts=config["pacbio"]["bbduk_entropy"]
    log:
        f"{OUT}/logs/pacbio/entropy.log"
    conda:
        "../envs/trim.yaml"
    shell:
        "bbduk.sh in={input} out={output} {params.opts} > {log} 2>&1"


rule pacbio_deepconsensus:
    # OPTIONAL. DeepConsensus re-calls each CCS read from its subreads. It
    # cannot work from CCS reads alone. Runs only when
    # pacbio.deepconsensus.enabled is true and subreads_bam is set.
    # The corrected reads replace step 04 for everything downstream.
    input:
        subreads=config["pacbio"]["deepconsensus"]["subreads_bam"] or [],
        keep=f"{PB}/clean/04_entropy.fastq.gz"
    output:
        f"{PB}/clean/05_deepconsensus.fastq.gz"
    params:
        ckpt=config["pacbio"]["deepconsensus"]["checkpoint"],
        work=f"{PB}/deepconsensus"
    threads: 16
    log:
        f"{OUT}/logs/pacbio/deepconsensus.log"
    conda:
        "../envs/deepconsensus.yaml"
    shell:
        "bash {SCRIPTS}/run_deepconsensus.sh {input.subreads} {input.keep} "
        "{params.ckpt} {params.work} {threads} {output} > {log} 2>&1"


def pacbio_clean_reads(wildcards):
    """Reads after cleaning: DeepConsensus output if enabled, else step 04."""
    if config["pacbio"]["deepconsensus"]["enabled"]:
        return f"{PB}/clean/05_deepconsensus.fastq.gz"
    return f"{PB}/clean/04_entropy.fastq.gz"


rule pacbio_rename:
    # Give reads short IDs (PB-000001, ...). PacBio read names contain "/"
    # characters that tree software rewrites. The table keeps the original names.
    input:
        pacbio_clean_reads
    output:
        fastq=f"{PB}/clean/reads.fastq.gz",
        fasta=f"{PB}/clean/reads.fasta",
        table=f"{PB}/clean/read_names.tsv"
    log:
        f"{OUT}/logs/pacbio/rename.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "python {SCRIPTS}/rename_reads.py --input {input} --prefix PB- "
        "--out-fastq {output.fastq} --out-fasta {output.fasta} "
        "--table {output.table} > {log} 2>&1"


rule pacbio_read_counts:
    # Read counts after each cleaning step (Table S3).
    input:
        raw=config["pacbio"]["ccs_fastq"],
        s1=f"{PB}/clean/01_passes.fastq.gz",
        s2=f"{PB}/clean/02_trimmed.fastq.gz",
        s3=f"{PB}/clean/03_oriented.fastq.gz",
        final=f"{PB}/clean/reads.fastq.gz"
    output:
        stats=f"{PB}/stats/cleaning.seqkit.tsv",
        metrics=f"{PB}/stats/cleaning.metrics.tsv"
    conda:
        "../envs/qc.yaml"
    shell:
        "seqkit stats -a -T {input.raw} {input.s1} {input.s2} {input.s3} {input.final} > {output.stats}; "
        "awk -F'\\t' 'NR>1{{split(\"raw passes trimmed oriented final\",n,\" \"); "
        "print \"pacbio_reads_\"n[NR-1]\"\\t\"$4; print \"pacbio_bases_\"n[NR-1]\"\\t\"$5}}' "
        "{output.stats} > {output.metrics}"


# ----------------------------- COX1 -----------------------------------------

rule pacbio_map_cox1:
    # Map the cleaned reads to the full-length Type AE* COX1 with mapPacBio.sh
    # ("high sensitivity" settings from manuscript_v1).
    input:
        reads=f"{PB}/clean/reads.fastq.gz",
        ref=config["pacbio"]["cox1_reference"]
    output:
        bam=f"{PB}/COX1/mappacbio.sorted.bam",
        flagstat=f"{PB}/COX1/mappacbio.flagstat.txt"
    params:
        opts=config["pacbio"]["mappacbio"],
        mem=config["pacbio"]["mappacbio_mem"],
        sam=f"{PB}/COX1/mappacbio.sam"
    threads: 16
    log:
        f"{OUT}/logs/pacbio/mappacbio.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        # nodisk: build the index in memory (no ref/ folder in the working directory).
        "mapPacBio.sh -Xmx{params.mem} in={input.reads} ref={input.ref} nodisk "
        "out={params.sam} {params.opts} threads={threads} > {log} 2>&1; "
        "samtools sort -@ {threads} -o {output.bam} {params.sam}; "
        "samtools index {output.bam}; rm -f {params.sam}; "
        "samtools flagstat {output.bam} > {output.flagstat}"


rule pacbio_minimap2:
    # minimap2 (-ax splice:hq -uf -C5), for comparison with mapPacBio.sh only.
    input:
        reads=f"{PB}/clean/reads.fastq.gz",
        ref=config["pacbio"]["cox1_reference"]
    output:
        flagstat=f"{PB}/COX1/minimap2.flagstat.txt"
    params:
        opts=config["pacbio"]["minimap2"]
    threads: 8
    log:
        f"{OUT}/logs/pacbio/minimap2.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "minimap2 {params.opts} -t {threads} {input.ref} {input.reads} 2> {log} "
        "| samtools flagstat - > {output.flagstat}"


rule pacbio_cox1_extract:
    # Take primary alignments, trim each read to the reference coordinates
    # (1-1527), then discard reads whose translation (code 5) contains stop
    # codons, which is how non-triplet indels (frameshifts) show up.
    input:
        bam=f"{PB}/COX1/mappacbio.sorted.bam"
    output:
        trimmed=f"{PB}/COX1/mapped_trimmed.fasta",
        kept=f"{PB}/COX1/mapped_inframe.fasta",
        report=f"{PB}/COX1/mapped_translation_report.tsv",
        metrics=f"{PB}/COX1/mapped.metrics.tsv"
    params:
        table=config["pacbio"]["cox1_table"]
    log:
        f"{OUT}/logs/pacbio/cox1_extract.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "python {SCRIPTS}/trim_to_reference.py --bam {input.bam} --out {output.trimmed} "
        "--min-len 100 --metrics {output.metrics} > {log} 2>&1; "
        "python {SCRIPTS}/check_translation.py filter --input {output.trimmed} "
        "--table {params.table} --frame-from-header --output {output.kept} "
        "--report {output.report} >> {log} 2>&1"


rule pacbio_cox1_rescue:
    # NEW: BLAST all cleaned reads against the Sanger COX1 library to find
    # divergent COX1 reads that did not map to the single Type AE* reference.
    # Hit regions are extracted, oriented and translation-checked (best frame).
    input:
        reads=f"{PB}/clean/reads.fasta",
        db=f"{OUT}/db/COX1.sanger.blastdb.ok",
        mapped=f"{PB}/COX1/mapped_inframe.fasta"
    output:
        hits=f"{PB}/COX1/rescue_blast.tsv",
        fasta=f"{PB}/COX1/rescued_inframe.fasta",
        metrics=f"{PB}/COX1/rescue.metrics.tsv"
    params:
        db=f"{OUT}/db/COX1.sanger",
        evalue=config["pacbio"]["cox1_blast_evalue"],
        table=config["pacbio"]["cox1_table"],
        tmp=f"{PB}/COX1/rescued_raw.fasta"
    threads: 16
    log:
        f"{OUT}/logs/pacbio/cox1_rescue.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn -query {input.reads} -db {params.db} -evalue {params.evalue} "
        "-max_target_seqs 5 -num_threads {threads} -dust no -soft_masking false "
        "-outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' "
        "> {output.hits} 2> {log}; "
        "python {SCRIPTS}/select_blast_hits.py --reads {input.reads} --hits {output.hits} "
        "--exclude {input.mapped} --evalue {params.evalue} --out {params.tmp} "
        "--metrics {output.metrics} --label pacbio_COX1_rescue >> {log} 2>&1; "
        "python {SCRIPTS}/check_translation.py filter --input {params.tmp} "
        "--table {params.table} --output {output.fasta} >> {log} 2>&1"


def pacbio_cox1_parts(wildcards):
    """Mapped reads, plus rescued reads if the rescue step is on."""
    parts = [f"{PB}/COX1/mapped_inframe.fasta"]
    if config["pacbio"]["cox1_blast_rescue"]:
        parts.append(f"{PB}/COX1/rescued_inframe.fasta")
    return parts


rule pacbio_cox1_final:
    # All PacBio COX1 sequences that go to the authentication filter.
    input:
        pacbio_cox1_parts
    output:
        f"{PB}/COX1/final.fasta"
    shell:
        "cat {input} > {output}"


# ----------------------------- 18S / 16S ------------------------------------

rule pacbio_blast_rrna:
    # BLAST the cleaned reads against SILVA (18S) or MIDORI lrRNA (16S) plus the
    # Sanger library (manuscript_v1: -num_alignments 1000 -outfmt 6 -dust no
    # -soft_masking false). Reads whose top hit is Facetotecta with
    # e-value <= 1e-50 are kept, and the hit region is extracted in the
    # orientation of the reference.
    input:
        reads=f"{PB}/clean/reads.fasta",
        db=f"{OUT}/db/{{locus}}.augmented.blastdb.ok",
        tax=f"{OUT}/db/{{locus}}.augmented.tax.tsv"
    output:
        hits=f"{PB}/{{locus}}/blast.tsv",
        fasta=f"{PB}/{{locus}}/final.fasta",
        metrics=f"{PB}/{{locus}}/blast.metrics.tsv"
    params:
        db=f"{OUT}/db/{{locus}}.augmented",
        evalue=config["pacbio"]["rrna_blast_evalue"],
        nal=config["pacbio"]["rrna_num_alignments"]
    wildcard_constraints:
        locus="|".join(PB_BLAST_LOCI)
    threads: 16
    log:
        f"{OUT}/logs/pacbio/{{locus}}.blast.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn -query {input.reads} -db {params.db} -evalue {params.evalue} "
        "-num_alignments {params.nal} -dust no -soft_masking false -num_threads {threads} "
        "-outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' "
        "> {output.hits} 2> {log}; "
        "python {SCRIPTS}/select_blast_hits.py --reads {input.reads} --hits {output.hits} "
        "--tax {input.tax} --require-taxon Facetotecta --evalue {params.evalue} "
        "--out {output.fasta} --metrics {output.metrics} "
        "--label pacbio_{wildcards.locus} >> {log} 2>&1"
