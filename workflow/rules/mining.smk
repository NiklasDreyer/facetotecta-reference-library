# =============================================================================
# mining.smk: find Facetotecta in public metabarcoding datasets (Section 2.14).
# Included only by workflow/Snakefile_mining. Needs the main workflow's results
# (Sanger library, authentic HTS sequences, reference databases).
#
# Per dataset (config/public_datasets.tsv):
#   download SRA runs -> cutadapt (primers) -> DADA2 ASVs ->
#   BLAST against the PUBLIC reference and against PUBLIC + OUR LIBRARY ->
#   keep ASVs whose top hit is Facetotecta and pass the identity / e-value
#   cut-offs (manuscript_v1: COX1 >92%, 1e-20; 16S >92%, 18S >90%, 1e-50)
# Per marker (all datasets pooled):
#   phylogenetic placement on the reference tree (LWR) ->
#   authentication filter together with our library ->
#   VSEARCH MOTUs -> gain table (ASVs that were unassigned or assigned to
#   something else before our library was added), map, Darwin Core export.
# =============================================================================

M = f"{OUT}/mining"

DATASETS = pd.read_csv(config["mining"]["datasets"], sep="\t", comment="#", dtype=str).set_index("dataset", drop=False)
DATASETS["run_list"] = DATASETS["runs"].str.split(";")
MINING_MARKERS = sorted(DATASETS["marker"].unique())
PE_SETS = DATASETS.index[DATASETS["layout"] == "PE"].tolist()
SE_SETS = DATASETS.index[DATASETS["layout"] == "SE"].tolist()


def runs_of(dataset):
    """SRA run accessions of a dataset."""
    return [r.strip() for r in DATASETS.loc[dataset, "run_list"] if r.strip()]


def datasets_of(marker):
    """Datasets for one marker."""
    return DATASETS.index[DATASETS["marker"] == marker].tolist()


# Downloads need internet access, which compute nodes often lack: run them on
# the login node (Snakemake runs local rules on the machine it was started on).
localrules: mining_fetch_pe, mining_fetch_se


rule mining_fetch_pe:
    # Download one paired-end run with prefetch + fasterq-dump.
    output:
        r1=f"{M}/raw/{{dataset}}/{{run}}_R1.fastq.gz",
        r2=f"{M}/raw/{{dataset}}/{{run}}_R2.fastq.gz"
    params:
        d=f"{M}/raw/{{dataset}}"
    wildcard_constraints:
        dataset="|".join(PE_SETS) if PE_SETS else "NONE"
    log:
        f"{OUT}/logs/mining/{{dataset}}.{{run}}.fetch.log"
    conda:
        "../envs/sra.yaml"
    shell:
        "bash {SCRIPTS}/fetch_sra.sh {wildcards.run} {params.d} PE > {log} 2>&1"


rule mining_fetch_se:
    # Download one single-end run.
    output:
        r1=f"{M}/raw/{{dataset}}/{{run}}_R1.fastq.gz"
    params:
        d=f"{M}/raw/{{dataset}}"
    wildcard_constraints:
        dataset="|".join(SE_SETS) if SE_SETS else "NONE"
    log:
        f"{OUT}/logs/mining/{{dataset}}.{{run}}.fetch.log"
    conda:
        "../envs/sra.yaml"
    shell:
        "bash {SCRIPTS}/fetch_sra.sh {wildcards.run} {params.d} SE > {log} 2>&1"


def mining_raw(wildcards):
    """Raw read files of one run (R1 only for single-end)."""
    d = f"{M}/raw/{wildcards.dataset}/{wildcards.run}"
    if DATASETS.loc[wildcards.dataset, "layout"] == "PE":
        return [f"{d}_R1.fastq.gz", f"{d}_R2.fastq.gz"]
    return [f"{d}_R1.fastq.gz"]


rule mining_cutadapt:
    # Remove the study's primers; read pairs without the primer are discarded.
    input:
        mining_raw
    output:
        r1=f"{M}/trimmed/{{dataset}}/{{run}}_R1.fastq.gz",
        r2=f"{M}/trimmed/{{dataset}}/{{run}}_R2.fastq.gz"
    params:
        fwd=lambda w: DATASETS.loc[w.dataset, "fwd"],
        rev=lambda w: DATASETS.loc[w.dataset, "rev"],
        layout=lambda w: DATASETS.loc[w.dataset, "layout"]
    threads: 4
    log:
        f"{OUT}/logs/mining/{{dataset}}.{{run}}.cutadapt.log"
    conda:
        "../envs/cutadapt.yaml"
    shell:
        "bash {SCRIPTS}/run_cutadapt.sh {params.layout} {params.fwd} {params.rev} "
        "{threads} {output.r1} {output.r2} {input} > {log} 2>&1"


rule mining_dada2:
    # ASVs for one dataset (same DADA2 script and settings as our own data).
    input:
        r1=lambda w: expand(f"{M}/trimmed/{w.dataset}/{{run}}_R1.fastq.gz", run=runs_of(w.dataset)),
        r2=lambda w: expand(f"{M}/trimmed/{w.dataset}/{{run}}_R2.fastq.gz", run=runs_of(w.dataset))
    output:
        fasta=f"{M}/asv/{{dataset}}/asvs.fasta",
        table=f"{M}/asv/{{dataset}}/asv_table.tsv",
        track=f"{M}/asv/{{dataset}}/track_reads.tsv",
        metrics=f"{M}/asv/{{dataset}}/dada2.metrics.tsv"
    params:
        r1=lambda w, input: ",".join(input.r1),
        # Single-end datasets: an empty R2 list switches DADA2 to forward reads.
        r2=lambda w, input: ",".join(input.r2) if DATASETS.loc[w.dataset, "layout"] == "PE" else "NONE",
        runs=lambda w: ",".join(runs_of(w.dataset)),
        filtdir=f"{M}/asv/{{dataset}}/filtered",
        merge=lambda w: "TRUE" if DATASETS.loc[w.dataset, "layout"] == "PE" else "FALSE"
    threads: 8
    log:
        f"{OUT}/logs/mining/{{dataset}}.dada2.log"
    conda:
        "../envs/dada2.yaml"
    shell:
        "Rscript {SCRIPTS}/dada2_asv.R --r1 {params.r1} --r2 {params.r2} "
        "--samples {params.runs} --filtdir {params.filtdir} --merge {params.merge} "
        "--prefix {wildcards.dataset}-ASV --label mining_{wildcards.dataset} "
        "--threads {threads} --out_fasta {output.fasta} --out_table {output.table} "
        "--out_track {output.track} --out_metrics {output.metrics} > {log} 2>&1"


rule mining_assign:
    # Two-pass BLAST assignment. Pass 1: public reference only (what the
    # original study could have found). Pass 2: public reference + our library.
    input:
        fasta=f"{M}/asv/{{dataset}}/asvs.fasta",
        db_pub=lambda w: f"{OUT}/db/{DATASETS.loc[w.dataset, 'marker']}.public.blastdb.ok",
        db_aug=lambda w: f"{OUT}/db/{DATASETS.loc[w.dataset, 'marker']}.augmented.blastdb.ok",
        tax_pub=lambda w: f"{OUT}/db/{DATASETS.loc[w.dataset, 'marker']}.public.tax.tsv",
        tax_aug=lambda w: f"{OUT}/db/{DATASETS.loc[w.dataset, 'marker']}.augmented.tax.tsv"
    output:
        table=f"{M}/assign/{{dataset}}/assignments.tsv",
        candidates=f"{M}/assign/{{dataset}}/candidates.fasta",
        metrics=f"{M}/assign/{{dataset}}/assign.metrics.tsv"
    params:
        marker=lambda w: DATASETS.loc[w.dataset, "marker"],
        db_pub=lambda w: f"{OUT}/db/{DATASETS.loc[w.dataset, 'marker']}.public",
        db_aug=lambda w: f"{OUT}/db/{DATASETS.loc[w.dataset, 'marker']}.augmented",
        ident=lambda w: config["mining"]["identity"][DATASETS.loc[w.dataset, "marker"]],
        evalue=lambda w: config["mining"]["evalue"][DATASETS.loc[w.dataset, "marker"]],
        workdir=f"{M}/assign/{{dataset}}"
    threads: 8
    log:
        f"{OUT}/logs/mining/{{dataset}}.assign.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "python {SCRIPTS}/assign_two_pass.py --fasta {input.fasta} "
        "--db-public {params.db_pub} --tax-public {input.tax_pub} "
        "--db-augmented {params.db_aug} --tax-augmented {input.tax_aug} "
        "--identity {params.ident} --evalue {params.evalue} --threads {threads} "
        "--workdir {params.workdir} --dataset {wildcards.dataset} "
        "--out-table {output.table} --out-candidates {output.candidates} "
        "--metrics {output.metrics} > {log} 2>&1"


rule mining_pool:
    # Pool candidate ASVs of all datasets for one marker, renamed
    # "PUB-<dataset>-<ASV>", as the query set for placement.
    input:
        lambda w: expand(f"{M}/assign/{{dataset}}/candidates.fasta", dataset=datasets_of(w.marker))
    output:
        f"{OUT}/placement/{{marker}}/mining/query.fasta"
    wildcard_constraints:
        marker="|".join(MINING_MARKERS)
    shell:
        # Candidate headers are already "<dataset>-ASV<n>"; add the PUB- tag.
        "cat {input} | sed 's/^>/>PUB-/' > {output}"


rule mining_authenticate:
    # Authentication filter on our authentic library + the public candidates
    # (same criteria as for our own Illumina/PacBio sequences).
    input:
        library=f"{OUT}/authenticate/{{marker}}/authentic.fasta",
        candidates=f"{OUT}/placement/{{marker}}/mining/query.fasta",
        outgroups=lambda w: os.path.join(config["authenticate"]["outgroups_dir"], f"{w.marker}.fasta")
    output:
        fasta=f"{M}/authenticate/{{marker}}/authentic.fasta",
        removed=f"{M}/authenticate/{{marker}}/removed.tsv",
        tree=f"{M}/authenticate/{{marker}}/final.treefile",
        aln=f"{M}/authenticate/{{marker}}/final.aln.fasta",
        metrics=f"{M}/authenticate/{{marker}}/authenticate.metrics.tsv"
    params:
        workdir=f"{M}/authenticate/{{marker}}/iterations",
        combined=f"{M}/authenticate/{{marker}}/input.fasta",
        min_shalrt=config["authenticate"]["min_shalrt"],
        min_ufboot=config["authenticate"]["min_ufboot"],
        long_branch=config["authenticate"]["long_branch_factor"],
        max_iter=config["authenticate"]["max_iterations"],
        ufboot=config["authenticate"]["ufboot"],
        alrt=config["authenticate"]["alrt"]
    wildcard_constraints:
        marker="|".join(MINING_MARKERS)
    threads: 16
    log:
        f"{OUT}/logs/mining/{{marker}}.authenticate.log"
    conda:
        "../envs/authenticate.yaml"
    shell:
        "cat {input.library} {input.candidates} > {params.combined}; "
        "python {SCRIPTS}/authenticate_filter.py --input {params.combined} "
        "--outgroups {input.outgroups} --workdir {params.workdir} --threads {threads} "
        "--min-shalrt {params.min_shalrt} --min-ufboot {params.min_ufboot} "
        "--long-branch {params.long_branch} --max-iter {params.max_iter} "
        "--ufboot {params.ufboot} --alrt {params.alrt} "
        "--out-fasta {output.fasta} --out-removed {output.removed} "
        "--out-tree {output.tree} --out-aln {output.aln} "
        "--metrics {output.metrics} --label mining_{wildcards.marker} > {log} 2>&1"


rule mining_motus:
    # VSEARCH MOTUs on library + authentic public sequences, at the locus threshold.
    input:
        f"{M}/authenticate/{{marker}}/authentic.fasta"
    output:
        uc=f"{M}/motus/{{marker}}/clusters.uc",
        members=f"{M}/motus/{{marker}}/membership.tsv"
    params:
        ident=lambda w: config["delimitation"]["vsearch_id"][w.marker],
        seqs=f"{M}/motus/{{marker}}/ungapped.fasta"
    wildcard_constraints:
        marker="|".join(MINING_MARKERS)
    conda:
        "../envs/vsearch.yaml"
    shell:
        "vsearch --fastx_filter {input} --fastaout {params.seqs} --quiet; "
        "vsearch --cluster_fast {params.seqs} --id {params.ident} --uc {output.uc} --quiet; "
        "python {SCRIPTS}/parse_uc.py --uc {output.uc} --output {output.members}"


rule mining_summary:
    # Per-dataset gain table, list of detections with coordinates, map, and
    # Darwin Core files for GBIF/OBIS.
    input:
        assign=expand(f"{M}/assign/{{dataset}}/assignments.tsv", dataset=DATASETS.index),
        placement=expand(f"{OUT}/placement/{{marker}}/mining/assignments.tsv", marker=MINING_MARKERS),
        authentic=expand(f"{M}/authenticate/{{marker}}/authentic.fasta", marker=MINING_MARKERS),
        members=expand(f"{M}/motus/{{marker}}/membership.tsv", marker=MINING_MARKERS),
        datasets=config["mining"]["datasets"]
    output:
        gain=f"{M}/summary/gain_by_dataset.tsv",
        detections=f"{M}/summary/detections.tsv",
        pdf=f"{M}/summary/detections_map.pdf",
        metrics=f"{M}/summary/mining.metrics.tsv"
    params:
        lwr=config["mining"]["min_lwr"]
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/mining_summary.py --assign {input.assign} "
        "--placement {input.placement} --authentic {input.authentic} "
        "--members {input.members} --datasets {input.datasets} --min-lwr {params.lwr} "
        "--out-gain {output.gain} --out-detections {output.detections} "
        "--out-pdf {output.pdf} --metrics {output.metrics}"


rule mining_dwc:
    # Darwin Core occurrence + DNA-derived data extension tables.
    input:
        detections=f"{M}/summary/detections.tsv",
        datasets=config["mining"]["datasets"]
    output:
        occ=f"{M}/dwc/occurrence.tsv",
        dna=f"{M}/dwc/dna_derived_data.tsv"
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/export_dwc.py --detections {input.detections} "
        "--datasets {input.datasets} --out-occurrence {output.occ} --out-dna {output.dna}"
