# =============================================================================
# richness.smk: NEW richness estimation and platform overlap (Section 2.11).
#
#   - abundance (specimens per MOTU per site) and incidence (MOTU presence per
#     sampling event) tables from the Sanger library
#   - iNEXT rarefaction/extrapolation, sample coverage and Chao estimates
#   - for each HTS locus: which platforms found each MOTU (UpSet plot)
# =============================================================================

R = f"{OUT}/richness"


rule richness_tables:
    # Build iNEXT input tables from the MOTU membership of the chosen dataset
    # and the specimen table.
    input:
        members=f"{OUT}/delim/{config['richness']['dataset']}/membership/{config['richness']['method']}.tsv",
        specimens=config["specimens"]
    output:
        abundance=f"{R}/abundance_by_site.tsv",
        incidence=f"{R}/incidence_by_site.tsv"
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/build_incidence.py --members {input.members} "
        "--specimens {input.specimens} --out-abundance {output.abundance} "
        "--out-incidence {output.incidence}"


rule richness_inext:
    # Rarefaction/extrapolation curves (q = 0), sample coverage and asymptotic
    # richness with 95% CIs, by site and for the whole library.
    input:
        abundance=f"{R}/abundance_by_site.tsv",
        incidence=f"{R}/incidence_by_site.tsv"
    output:
        estimates=f"{R}/inext_estimates.tsv",
        pdf=f"{R}/inext_curves.pdf",
        metrics=f"{R}/richness.metrics.tsv"
    params:
        nboot=config["richness"]["nboot"],
        factor=config["richness"]["endpoint_factor"]
    log:
        f"{OUT}/logs/richness/inext.log"
    conda:
        "../envs/inext.yaml"
    shell:
        "Rscript {SCRIPTS}/richness_inext.R --abundance {input.abundance} "
        "--incidence {input.incidence} --nboot {params.nboot} "
        "--endpoint_factor {params.factor} --out_estimates {output.estimates} "
        "--out_pdf {output.pdf} --out_metrics {output.metrics} > {log} 2>&1"


rule platform_overlap:
    # For each MOTU of the combined HTS dataset, list the platforms/samples in
    # which it was found (Sanger, PacBio, Illumina Y1, Illumina Y2) and draw an
    # UpSet plot (Fig. 4b).
    input:
        members=f"{OUT}/delim/hts_{{locus}}_all/membership/vsearch.tsv",
        platforms=f"{OUT}/authenticate/{{locus}}/seq_platforms.tsv"
    output:
        table=f"{R}/overlap_{{locus}}.tsv",
        pdf=f"{R}/overlap_{{locus}}_upset.pdf",
        metrics=f"{R}/overlap_{{locus}}.metrics.tsv"
    wildcard_constraints:
        locus="|".join(HTS_LOCI)
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/platform_overlap.py --members {input.members} "
        "--platforms {input.platforms} --out-table {output.table} "
        "--out-pdf {output.pdf} --metrics {output.metrics} --label {wildcards.locus}"
