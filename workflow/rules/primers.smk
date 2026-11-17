# =============================================================================
# primers.smk: NEW in-silico primer coverage (Section 2.13).
#
# For each primer pair in config/insilico_primers.tsv, count mismatches against
# every facetotectan reference sequence (Sanger + PacBio + Type AE* COX1) and
# against the non-target crustaceans, and extract the amplicon region (used by
# the benchmark on short fragments).
# =============================================================================

PRIMERS = pd.read_csv(config["primers"]["table"], sep="\t", comment="#", dtype=str).set_index("name", drop=False)
PRIMER_NAMES = PRIMERS.index.tolist()


def primers_of(marker):
    """Primer pairs targeting one marker."""
    return PRIMERS.loc[PRIMERS["marker"] == marker, "name"].tolist()


rule primer_reference_set:
    # All facetotectan sequences available for a marker: the Sanger library,
    # the PacBio reads (long, so they span more primer sites) and, for COX1,
    # the full-length Type AE* COX1.
    input:
        sanger=f"{OUT}/sanger/{{marker}}/{{marker}}.clean.fasta",
        pacbio=lambda w: f"{OUT}/pacbio/{w.marker}/final.fasta" if pacbio_has(w.marker) else [],
        ref=lambda w: config["pacbio"]["cox1_reference"] if w.marker == "COX1" else []
    output:
        f"{OUT}/primers/reference_{{marker}}.fasta"
    shell:
        "cat {input} > {output}"


rule primer_scan:
    # Mismatch table, coverage summary and amplicons for one primer pair,
    # against facetotectan references ("facet") or non-targets ("nontarget").
    input:
        seqs=lambda w: (f"{OUT}/primers/reference_{PRIMERS.loc[w.pp, 'marker']}.fasta"
                        if w.set == "facet"
                        else f"{OUT}/benchmark/{PRIMERS.loc[w.pp, 'marker']}/nontarget.fasta")
    output:
        table=f"{OUT}/primers/{{pp}}/{{set}}_mismatches.tsv",
        amplicons=f"{OUT}/primers/{{pp}}/{{set}}_amplicons.fasta",
        metrics=f"{OUT}/primers/{{pp}}/{{set}}.metrics.tsv"
    params:
        fwd=lambda w: PRIMERS.loc[w.pp, "forward"],
        rev=lambda w: PRIMERS.loc[w.pp, "reverse"],
        maxmm=config["primers"]["max_mismatch"],
        three=config["primers"]["three_prime_exact"]
    wildcard_constraints:
        pp="|".join(PRIMER_NAMES),
        set="facet|nontarget"
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/primer_scan.py --input {input.seqs} --forward {params.fwd} "
        "--reverse {params.rev} --max-mismatch {params.maxmm} --three-prime {params.three} "
        "--out-table {output.table} --out-amplicons {output.amplicons} "
        "--metrics {output.metrics} --label primer_{wildcards.pp}_{wildcards.set}"


rule primer_plot:
    # Fig. 6b: mismatches per primer position, all primer pairs together.
    input:
        expand(f"{OUT}/primers/{{pp}}/facet_mismatches.tsv", pp=PRIMER_NAMES)
    output:
        f"{OUT}/primers/primer_mismatch_summary.pdf"
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/plot_primers.py --tables {input} --output {output}"
