# =============================================================================
# illumina.smk: Illumina MiSeq amplicons (COX1, 18S) from bulk samples Y1 and Y2
# (Sections 2.4, 2.7).
#
# Per library (sample x marker):
#   FastQC -> Trimmomatic -> BBDuk -> FastQC
# Per marker (all libraries together):
#   DADA2 ASVs -> BLAST contaminant screen (MIDORI/SILVA + Sanger library)
#   COX1 only: TransDecoder translation screen -> metaMATE NUMT filter (NEW)
# =============================================================================


rule fastqc_raw:
    # Quality report of the raw reads.
    input:
        r1=lambda w: UNITS.loc[w.unit, "r1"],
        r2=lambda w: UNITS.loc[w.unit, "r2"]
    output:
        directory(f"{OUT}/illumina/qc/raw/{{unit}}")
    threads: 2
    log:
        f"{OUT}/logs/illumina/{{unit}}.fastqc_raw.log"
    conda:
        "../envs/qc.yaml"
    shell:
        "mkdir -p {output} && fastqc -t {threads} -o {output} {input.r1} {input.r2} > {log} 2>&1"


rule trimmomatic:
    # First trimming pass, as in manuscript_v1:
    # PE -phred33 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:3
    # The TruSeq3-PE adapter file ships with the Trimmomatic conda package;
    # the shell code finds it inside the environment.
    input:
        r1=lambda w: UNITS.loc[w.unit, "r1"],
        r2=lambda w: UNITS.loc[w.unit, "r2"]
    output:
        r1=f"{OUT}/illumina/trim/{{unit}}.trimmomatic.R1.fastq.gz",
        r2=f"{OUT}/illumina/trim/{{unit}}.trimmomatic.R2.fastq.gz",
        u1=f"{OUT}/illumina/trim/{{unit}}.trimmomatic.U1.fastq.gz",
        u2=f"{OUT}/illumina/trim/{{unit}}.trimmomatic.U2.fastq.gz"
    params:
        clip=config["illumina"]["trimmomatic_clip"],
        steps=config["illumina"]["trimmomatic_steps"]
    threads: 4
    log:
        f"{OUT}/logs/illumina/{{unit}}.trimmomatic.log"
    conda:
        "../envs/trim.yaml"
    shell:
        # Locate TruSeq3-PE.fa in the conda environment.
        "ADAPT=$(find $CONDA_PREFIX/share -name TruSeq3-PE.fa | head -n 1); "
        "trimmomatic PE -threads {threads} -phred33 "
        "{input.r1} {input.r2} {output.r1} {output.u1} {output.r2} {output.u2} "
        "ILLUMINACLIP:$ADAPT:{params.clip} {params.steps} > {log} 2>&1"


rule bbduk_illumina:
    # Second trimming pass with BBDuk, as in manuscript_v1:
    # qtrim=rl trimq=33 minlen=100 ktrim=r k=23 mink=11 hdist=1 tpe tbo
    # "ref=adapters" uses the adapter set built into BBTools.
    input:
        r1=f"{OUT}/illumina/trim/{{unit}}.trimmomatic.R1.fastq.gz",
        r2=f"{OUT}/illumina/trim/{{unit}}.trimmomatic.R2.fastq.gz"
    output:
        r1=f"{OUT}/illumina/trim/{{unit}}.bbduk.R1.fastq.gz",
        r2=f"{OUT}/illumina/trim/{{unit}}.bbduk.R2.fastq.gz"
    params:
        opts=config["illumina"]["bbduk"]
    threads: 4
    log:
        f"{OUT}/logs/illumina/{{unit}}.bbduk.log"
    conda:
        "../envs/trim.yaml"
    shell:
        "bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} "
        "ref=adapters {params.opts} threads={threads} > {log} 2>&1"


rule fastqc_trimmed:
    # Quality report after both trimming passes.
    input:
        r1=f"{OUT}/illumina/trim/{{unit}}.bbduk.R1.fastq.gz",
        r2=f"{OUT}/illumina/trim/{{unit}}.bbduk.R2.fastq.gz"
    output:
        directory(f"{OUT}/illumina/qc/trimmed/{{unit}}")
    threads: 2
    log:
        f"{OUT}/logs/illumina/{{unit}}.fastqc_trimmed.log"
    conda:
        "../envs/qc.yaml"
    shell:
        "mkdir -p {output} && fastqc -t {threads} -o {output} {input.r1} {input.r2} > {log} 2>&1"


rule illumina_read_counts:
    # Read counts (and mean quality) at each stage for Table S3.
    input:
        raw1=lambda w: UNITS.loc[w.unit, "r1"],
        tm1=f"{OUT}/illumina/trim/{{unit}}.trimmomatic.R1.fastq.gz",
        bb1=f"{OUT}/illumina/trim/{{unit}}.bbduk.R1.fastq.gz"
    output:
        stats=f"{OUT}/illumina/stats/{{unit}}.seqkit.tsv",
        metrics=f"{OUT}/illumina/stats/{{unit}}.metrics.tsv"
    conda:
        "../envs/qc.yaml"
    shell:
        # seqkit stats -a gives read number, length and average quality.
        "seqkit stats -a -T {input.raw1} {input.tm1} {input.bb1} > {output.stats}; "
        # Turn it into key/value metrics (read pairs at each stage).
        "awk -F'\\t' -v u={wildcards.unit} 'NR==2{{print \"illumina_\"u\"_pairs_raw\\t\"$4}} "
        "NR==3{{print \"illumina_\"u\"_pairs_trimmomatic\\t\"$4}} "
        "NR==4{{print \"illumina_\"u\"_pairs_bbduk\\t\"$4}}' {output.stats} > {output.metrics}"


rule multiqc_illumina:
    # One HTML report combining all FastQC and Trimmomatic outputs.
    input:
        expand(f"{OUT}/illumina/qc/raw/{{unit}}", unit=UNITS.index),
        expand(f"{OUT}/illumina/qc/trimmed/{{unit}}", unit=UNITS.index)
    output:
        f"{OUT}/illumina/qc/multiqc_report.html"
    params:
        outdir=f"{OUT}/illumina/qc"
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc -f -o {params.outdir} {input} {OUT}/logs/illumina"


rule dada2:
    # Infer ASVs with DADA2, following the developer's tutorial: filtering
    # (maxN=0, maxEE, truncQ, PhiX removal), error learning, denoising, read
    # merging, sequence table and chimera removal (removeBimeraDenovo).
    # ASVs are named "ILL-ASV<n>".
    input:
        r1=lambda w: expand(f"{OUT}/illumina/trim/{{unit}}.bbduk.R1.fastq.gz", unit=units_of(w.marker)),
        r2=lambda w: expand(f"{OUT}/illumina/trim/{{unit}}.bbduk.R2.fastq.gz", unit=units_of(w.marker))
    output:
        fasta=f"{OUT}/illumina/{{marker}}/dada2/asvs.fasta",
        table=f"{OUT}/illumina/{{marker}}/dada2/asv_table.tsv",
        track=f"{OUT}/illumina/{{marker}}/dada2/track_reads.tsv",
        metrics=f"{OUT}/illumina/{{marker}}/dada2/dada2.metrics.tsv"
    params:
        # optparse takes one value per option, so file lists are comma-joined.
        r1=lambda w, input: ",".join(input.r1),
        r2=lambda w, input: ",".join(input.r2),
        samples=lambda w: ",".join(samples_of(w.marker)),
        filtdir=f"{OUT}/illumina/{{marker}}/dada2/filtered",
        merge=str(config["illumina"]["dada2"]["merge_pairs"]).upper(),
        trunc_f=config["illumina"]["dada2"]["trunc_len_f"],
        trunc_r=config["illumina"]["dada2"]["trunc_len_r"],
        maxee_f=config["illumina"]["dada2"]["max_ee_f"],
        maxee_r=config["illumina"]["dada2"]["max_ee_r"],
        truncq=config["illumina"]["dada2"]["trunc_q"],
        min_overlap=config["illumina"]["dada2"]["min_overlap"],
        chimera=config["illumina"]["dada2"]["chimera_method"]
    threads: 8
    log:
        f"{OUT}/logs/illumina/{{marker}}.dada2.log"
    conda:
        "../envs/dada2.yaml"
    shell:
        "Rscript {SCRIPTS}/dada2_asv.R "
        "--r1 {params.r1} --r2 {params.r2} --samples {params.samples} "
        "--filtdir {params.filtdir} --merge {params.merge} "
        "--trunc_f {params.trunc_f} --trunc_r {params.trunc_r} "
        "--maxee_f {params.maxee_f} --maxee_r {params.maxee_r} "
        "--truncq {params.truncq} --min_overlap {params.min_overlap} "
        "--chimera {params.chimera} --prefix ILL-ASV "
        "--label illumina_{wildcards.marker} --threads {threads} "
        "--out_fasta {output.fasta} --out_table {output.table} "
        "--out_track {output.track} --out_metrics {output.metrics} > {log} 2>&1"


rule blast_screen:
    # Contaminant screen (manuscript_v1): BLAST the ASVs against MIDORI (COX1)
    # or SILVA (18S) plus the Sanger library; keep ASVs whose top hit is
    # Facetotecta or Arthropoda.
    input:
        fasta=f"{OUT}/illumina/{{marker}}/dada2/asvs.fasta",
        db=f"{OUT}/db/{{marker}}.augmented.blastdb.ok",
        tax=f"{OUT}/db/{{marker}}.augmented.tax.tsv"
    output:
        hits=f"{OUT}/illumina/{{marker}}/screen/blast_hits.tsv",
        fasta=f"{OUT}/illumina/{{marker}}/screen/screened.fasta",
        report=f"{OUT}/illumina/{{marker}}/screen/screen_report.tsv",
        metrics=f"{OUT}/illumina/{{marker}}/screen/screen.metrics.tsv"
    params:
        db=f"{OUT}/db/{{marker}}.augmented",
        keep=" ".join(config["illumina"]["screen_keep_taxa"]),
        evalue=config["illumina"]["screen_evalue"]
    threads: 8
    log:
        f"{OUT}/logs/illumina/{{marker}}.blast_screen.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn -query {input.fasta} -db {params.db} -evalue {params.evalue} "
        "-max_target_seqs 20 -num_threads {threads} "
        "-outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' "
        "> {output.hits} 2> {log}; "
        "python {SCRIPTS}/blast_screen.py --fasta {input.fasta} --hits {output.hits} "
        "--tax {input.tax} --keep {params.keep} --out-fasta {output.fasta} "
        "--report {output.report} --metrics {output.metrics} "
        "--label illumina_{wildcards.marker} >> {log} 2>&1"


rule transdecoder_cox1:
    # Translate COX1 ASVs with TransDecoder (invertebrate mitochondrial code)
    # and keep ASVs with an open reading frame covering (almost) the whole ASV,
    # i.e. no stop codons or frameshifts (manuscript_v1).
    input:
        f"{OUT}/illumina/COX1/screen/screened.fasta"
    output:
        fasta=f"{OUT}/illumina/COX1/transdecoder/translated_ok.fasta",
        metrics=f"{OUT}/illumina/COX1/transdecoder/transdecoder.metrics.tsv"
    params:
        outdir=f"{OUT}/illumina/COX1/transdecoder",
        code=config["illumina"]["transdecoder"]["genetic_code"],
        minlen=config["illumina"]["transdecoder"]["min_protein_len"],
        frac=config["illumina"]["transdecoder"]["min_orf_fraction"]
    log:
        f"{OUT}/logs/illumina/COX1.transdecoder.log"
    conda:
        "../envs/transdecoder.yaml"
    shell:
        # -S: the ASVs are in sense orientation (forward primer = start of COX1 fragment).
        "TransDecoder.LongOrfs -t {input} -G {params.code} -m {params.minlen} -S "
        "-O {params.outdir}/td > {log} 2>&1; "
        "python {SCRIPTS}/filter_transdecoder.py --fasta {input} "
        "--pep $(find {params.outdir}/td -name longest_orfs.pep | head -n 1) "
        "--min-fraction {params.frac} --out {output.fasta} "
        "--metrics {output.metrics} >> {log} 2>&1"


rule metamate_cox1:
    # NEW: metaMATE NUMT/error filtering. ASVs matching the curated Sanger COX1
    # library are the "verified authentic" set; metaMATE tests a grid of
    # abundance thresholds and estimates how many authentic and non-authentic
    # ASVs each removes. The chosen threshold set is then applied.
    input:
        fasta=f"{OUT}/illumina/COX1/transdecoder/translated_ok.fasta",
        table=f"{OUT}/illumina/COX1/dada2/asv_table.tsv",
        ref=f"{OUT}/sanger/COX1/COX1.clean.fasta",
        specs=config["illumina"]["metamate"]["specs"]
    output:
        fasta=f"{OUT}/illumina/COX1/metamate/metamate_kept.fasta",
        metrics=f"{OUT}/illumina/COX1/metamate/metamate.metrics.tsv"
    params:
        outdir=f"{OUT}/illumina/COX1/metamate",
        length=config["illumina"]["metamate"]["expected_length"],
        index=config["illumina"]["metamate"]["result_index"]
    threads: 8
    log:
        f"{OUT}/logs/illumina/COX1.metamate.log"
    conda:
        "../envs/metamate.yaml"
    shell:
        "bash {SCRIPTS}/run_metamate.sh {input.fasta} {input.table} {input.ref} "
        "{input.specs} {params.length} {params.index} {params.outdir} {threads} "
        "{output.fasta} {output.metrics} > {log} 2>&1"


def illumina_final_input(wildcards):
    """Last filtering step for each marker: metaMATE (if on) or TransDecoder
    for COX1; the BLAST screen for rRNA markers."""
    if wildcards.marker == "COX1":
        if config["illumina"]["metamate"]["enabled"]:
            return f"{OUT}/illumina/COX1/metamate/metamate_kept.fasta"
        return f"{OUT}/illumina/COX1/transdecoder/translated_ok.fasta"
    return f"{OUT}/illumina/{wildcards.marker}/screen/screened.fasta"


rule illumina_final:
    # The ASVs that go forward to the combined authentication filter.
    input:
        illumina_final_input
    output:
        f"{OUT}/illumina/{{marker}}/final.fasta"
    shell:
        "cp {input} {output}"
