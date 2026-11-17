# =============================================================================
# benchmark.smk: NEW benchmarking of taxonomic assignment (Section 2.12).
#
#   1. Non-target set: crustacean groups other than Facetotecta, sampled from
#      the public reference (MIDORI2 COX1, SILVA 18S).
#   2. Leave-one-MOTU-out cross-validation with BLAST (top hit) and SINTAX, on
#      full-length references ("full") and on the regions amplified by each
#      primer pair (from primers.smk). Output: true- and false-positive rates
#      for every identity threshold, ROC curves, and the threshold meeting the
#      target false-positive rate.
#   3. Placement of non-target sequences on the Facetotecta + outgroup tree
#      (how often a non-target would be placed inside Facetotecta).
# =============================================================================

B = f"{OUT}/benchmark"
BENCH_MARKERS = config["benchmark"]["markers"]


def bench_regions(marker):
    """'full' plus every primer pair for the marker."""
    return ["full"] + primers_of(marker)


rule extract_nontarget:
    # Sample up to max_per_group sequences of each non-target group.
    input:
        fasta=f"{OUT}/db/{{marker}}.public.fasta",
        tax=f"{OUT}/db/{{marker}}.public.tax.tsv"
    output:
        fasta=f"{B}/{{marker}}/nontarget.fasta",
        groups=f"{B}/{{marker}}/nontarget_groups.tsv"
    params:
        groups=" ".join(config["benchmark"]["nontarget_groups"]),
        maxn=config["benchmark"]["max_per_group"],
        seed=config["benchmark"]["seed"]
    log:
        f"{OUT}/logs/benchmark/{{marker}}.nontarget.log"
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/extract_nontarget.py --fasta {input.fasta} --tax {input.tax} "
        "--groups {params.groups} --max-per-group {params.maxn} --seed {params.seed} "
        "--out-fasta {output.fasta} --out-groups {output.groups} > {log} 2>&1"


def bench_facet(wildcards):
    """Facetotectan query/reference sequences for a region."""
    if wildcards.region == "full":
        return f"{OUT}/sanger/{wildcards.marker}/{wildcards.marker}.clean.fasta"
    return f"{OUT}/primers/{wildcards.region}/facet_amplicons.fasta"


def bench_nontarget(wildcards):
    """Non-target sequences for a region."""
    if wildcards.region == "full":
        return f"{B}/{wildcards.marker}/nontarget.fasta"
    return f"{OUT}/primers/{wildcards.region}/nontarget_amplicons.fasta"


rule benchmark_loo:
    # Leave-one-MOTU-out cross-validation (see benchmark_loo.py).
    input:
        facet=bench_facet,
        nontarget=bench_nontarget,
        groups=f"{B}/{{marker}}/nontarget_groups.tsv",
        members=f"{OUT}/delim/sanger_{{marker}}/membership/vsearch.tsv"
    output:
        hits=f"{B}/{{marker}}/{{region}}/loo_hits.tsv",
        roc=f"{B}/{{marker}}/{{region}}/roc.tsv",
        pdf=f"{B}/{{marker}}/{{region}}/roc.pdf",
        thresholds=f"{B}/{{marker}}/{{region}}/thresholds.tsv",
        metrics=f"{B}/{{marker}}/{{region}}/benchmark.metrics.tsv"
    params:
        workdir=f"{B}/{{marker}}/{{region}}/work",
        folds=config["benchmark"]["folds"],
        tmin=config["benchmark"]["threshold_min"],
        tmax=config["benchmark"]["threshold_max"],
        fpr=config["benchmark"]["target_fpr"],
        seed=config["benchmark"]["seed"]
    threads: 8
    log:
        f"{OUT}/logs/benchmark/{{marker}}.{{region}}.loo.log"
    conda:
        "../envs/benchmark.yaml"
    shell:
        "python {SCRIPTS}/benchmark_loo.py --facet {input.facet} --nontarget {input.nontarget} "
        "--groups {input.groups} --members {input.members} --workdir {params.workdir} "
        "--folds {params.folds} --tmin {params.tmin} --tmax {params.tmax} "
        "--target-fpr {params.fpr} --seed {params.seed} --threads {threads} "
        "--out-hits {output.hits} --out-roc {output.roc} --out-pdf {output.pdf} "
        "--out-thresholds {output.thresholds} --metrics {output.metrics} "
        "--label bench_{wildcards.marker}_{wildcards.region} > {log} 2>&1"


rule benchmark_placement_query:
    # A random sample of 500 non-target sequences to place on the reference tree.
    input:
        f"{B}/{{marker}}/nontarget.fasta"
    output:
        f"{OUT}/placement/{{marker}}/nontarget/query.fasta"
    params:
        seed=config["benchmark"]["seed"]
    conda:
        "../envs/qc.yaml"
    shell:
        # Shuffle with a fixed seed, then take the first 500 (two steps, so no
        # broken pipe under "set -o pipefail").
        "seqkit shuffle -s {params.seed} {input} > {output}.tmp; "
        "seqkit head -n 500 {output}.tmp > {output}; rm -f {output}.tmp"


rule benchmark_placement_summary:
    # How many non-target sequences were placed inside Facetotecta (LWR >= 0.9).
    input:
        assign=f"{OUT}/placement/{{marker}}/nontarget/assignments.tsv",
        groups=f"{B}/{{marker}}/nontarget_groups.tsv"
    output:
        f"{B}/{{marker}}/placement_nontarget.metrics.tsv"
    params:
        lwr=config["mining"]["min_lwr"]
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/summarize_placement.py --assign {input.assign} "
        "--groups {input.groups} --min-lwr {params.lwr} "
        "--label bench_{wildcards.marker}_placement --output {output}"
