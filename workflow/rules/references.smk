# =============================================================================
# references.smk: reference databases for screening and classification.
#
# For each high-throughput marker (COX1, 18S, 16S) two BLAST databases are built:
#   public    = public reference only (MIDORI2 for COX1/16S, SILVA for 18S)
#   augmented = public reference + the curated Sanger library of this study
# The "augmented" database is used for the Illumina contaminant screen and the
# PacBio rRNA BLAST (as in manuscript_v1). The "public" database is the baseline
# for measuring what the new library adds in public data (mining workflow).
# =============================================================================


def public_source(wildcards):
    """Public reference FASTA for a marker."""
    m = wildcards.marker
    if m == "COX1":
        return config["references"]["midori_cox1"]
    if m == "16S":
        return config["references"]["midori_16s"]
    if m == "18S":
        return config["references"]["silva_ssu"]
    raise ValueError(f"No public reference defined for marker {m}")


rule ref_public:
    # Reformat the public reference into a clean FASTA plus a taxonomy table
    # (id <tab> lineage). For 18S only eukaryotes are kept (manuscript_v1:
    # "all eukaryotic 18S rDNA sequences from SILVA").
    input:
        public_source
    output:
        fasta=f"{OUT}/db/{{marker}}.public.fasta",
        tax=f"{OUT}/db/{{marker}}.public.tax.tsv"
    params:
        fmt=lambda w: "silva" if w.marker == "18S" else "midori"
    log:
        f"{OUT}/logs/db/{{marker}}.public.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "python {SCRIPTS}/build_reference_db.py public "
        "--input {input} --format {params.fmt} "
        "--out-fasta {output.fasta} --out-tax {output.tax} > {log} 2>&1"


rule ref_augmented:
    # Add the curated Sanger sequences of this study, labelled as Facetotecta,
    # to the public reference.
    input:
        fasta=f"{OUT}/db/{{marker}}.public.fasta",
        tax=f"{OUT}/db/{{marker}}.public.tax.tsv",
        sanger=f"{OUT}/sanger/{{marker}}/{{marker}}.clean.fasta"
    output:
        fasta=f"{OUT}/db/{{marker}}.augmented.fasta",
        tax=f"{OUT}/db/{{marker}}.augmented.tax.tsv"
    log:
        f"{OUT}/logs/db/{{marker}}.augmented.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "python {SCRIPTS}/build_reference_db.py augment "
        "--input {input.fasta} --tax {input.tax} --library {input.sanger} "
        "--out-fasta {output.fasta} --out-tax {output.tax} > {log} 2>&1"


rule makeblastdb:
    # Build a nucleotide BLAST database. The ".ok" file marks completion,
    # because BLAST writes several files whose names depend on its version.
    input:
        f"{OUT}/db/{{marker}}.{{kind}}.fasta"
    output:
        touch(f"{OUT}/db/{{marker}}.{{kind}}.blastdb.ok")
    params:
        prefix=f"{OUT}/db/{{marker}}.{{kind}}"
    wildcard_constraints:
        kind="public|augmented|sanger"
    log:
        f"{OUT}/logs/db/{{marker}}.{{kind}}.makeblastdb.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype nucl -parse_seqids "
        "-out {params.prefix} > {log} 2>&1"


rule ref_sanger_only:
    # A BLAST database of the Sanger library alone (used to rescue divergent
    # PacBio COX1 reads).
    input:
        f"{OUT}/sanger/{{marker}}/{{marker}}.clean.fasta"
    output:
        f"{OUT}/db/{{marker}}.sanger.fasta"
    shell:
        "cp {input} {output}"
