# =============================================================================
# summary.smk: collect every number quoted in the manuscript into one table.
#
# Every step writes a small "*.metrics.tsv" file (key <tab> value). This rule
# runs last and merges them into results/summary_numbers.tsv. Fill the
# manuscript from this file only, so text and results cannot disagree.
# =============================================================================


rule summary_numbers:
    input:
        # All final targets (defined in the Snakefile) must exist first.
        ALL_TARGETS
    output:
        f"{OUT}/summary_numbers.tsv"
    params:
        outdir=OUT
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/summarize_numbers.py --results {params.outdir} --output {output}"
