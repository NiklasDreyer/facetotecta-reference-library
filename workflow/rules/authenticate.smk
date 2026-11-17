# =============================================================================
# authenticate.smk: combine Sanger, Illumina and PacBio sequences per locus and
# keep only authentic facetotectan sequences (Section 2.9).
#
#   1. Pre-cluster each platform separately with VSEARCH (precluster_id).
#   2. Iterative phylogenetic filter: align with outgroups (MAFFT), Gblocks,
#      IQ-TREE (ModelFinder, UFBoot -bnni, SH-aLRT), root on outgroups, remove
#      high-throughput sequences that
#        (1) fall outside a well-supported monophyletic Facetotecta,
#        (2) sit on a branch > 10x longer than their sister,
#        (3) form long-branch clades with reads from one HTS platform only and
#            no Sanger sequence,
#        (4) optionally: have no other platform in their closest supported clade.
#      Repeat until nothing is removed. Sanger sequences are never removed.
#   3. Write per-platform subsets of the authentic set (for per-platform MOTUs).
# =============================================================================


def locus_platform_inputs(wildcards):
    """FASTA files that feed a locus, by platform."""
    l = wildcards.locus
    files = {"SAN": f"{OUT}/sanger/{l}/{l}.clean.fasta"}
    if pacbio_has(l):
        files["PB"] = f"{OUT}/pacbio/{l}/final.fasta"
    if l in ILL_MARKERS:
        files["ILL"] = f"{OUT}/illumina/{l}/final.fasta"
    return files


rule precluster:
    # Pre-cluster each platform separately (so each centroid keeps its platform
    # tag) and record which original sequences each centroid stands for.
    input:
        unpack(locus_platform_inputs)
    output:
        fasta=f"{OUT}/authenticate/{{locus}}/centroids.fasta",
        members=f"{OUT}/authenticate/{{locus}}/precluster_members.tsv",
        metrics=f"{OUT}/authenticate/{{locus}}/precluster.metrics.tsv"
    params:
        ident=config["authenticate"]["precluster_id"],
        inputs=lambda w: " ".join(f"{k}={v}" for k, v in locus_platform_inputs(w).items())
    wildcard_constraints:
        locus="|".join(HTS_LOCI)
    threads: 4
    log:
        f"{OUT}/logs/authenticate/{{locus}}.precluster.log"
    conda:
        "../envs/authenticate.yaml"
    shell:
        "python {SCRIPTS}/precluster.py --inputs {params.inputs} --id {params.ident} "
        "--threads {threads} --out-fasta {output.fasta} --out-members {output.members} "
        "--metrics {output.metrics} --label {wildcards.locus} > {log} 2>&1"


rule authenticate:
    # The iterative phylogenetic filter (runs MAFFT, Gblocks and IQ-TREE in a
    # loop inside the script; see authenticate_filter.py).
    input:
        fasta=f"{OUT}/authenticate/{{locus}}/centroids.fasta",
        outgroups=lambda w: os.path.join(config["authenticate"]["outgroups_dir"], f"{w.locus}.fasta")
    output:
        fasta=f"{OUT}/authenticate/{{locus}}/authentic.fasta",
        removed=f"{OUT}/authenticate/{{locus}}/removed.tsv",
        tree=f"{OUT}/authenticate/{{locus}}/final.treefile",
        aln=f"{OUT}/authenticate/{{locus}}/final.aln.fasta",
        metrics=f"{OUT}/authenticate/{{locus}}/authenticate.metrics.tsv"
    params:
        workdir=f"{OUT}/authenticate/{{locus}}/iterations",
        min_shalrt=config["authenticate"]["min_shalrt"],
        min_ufboot=config["authenticate"]["min_ufboot"],
        long_branch=config["authenticate"]["long_branch_factor"],
        max_iter=config["authenticate"]["max_iterations"],
        ufboot=config["authenticate"]["ufboot"],
        alrt=config["authenticate"]["alrt"],
        cross="--require-cross-platform" if config["authenticate"]["require_cross_platform"] else ""
    wildcard_constraints:
        locus="|".join(HTS_LOCI)
    threads: 16
    log:
        f"{OUT}/logs/authenticate/{{locus}}.filter.log"
    conda:
        "../envs/authenticate.yaml"
    shell:
        "python {SCRIPTS}/authenticate_filter.py --input {input.fasta} "
        "--outgroups {input.outgroups} --workdir {params.workdir} "
        "--threads {threads} --min-shalrt {params.min_shalrt} "
        "--min-ufboot {params.min_ufboot} --long-branch {params.long_branch} "
        "--max-iter {params.max_iter} --ufboot {params.ufboot} "
        "--alrt {params.alrt} {params.cross} "
        "--out-fasta {output.fasta} --out-removed {output.removed} "
        "--out-tree {output.tree} --out-aln {output.aln} "
        "--metrics {output.metrics} --label {wildcards.locus} > {log} 2>&1"


# One subset rule per locus, because the set of platforms (and so the output
# files) differs between loci.
for _l in HTS_LOCI:
    _plats = [p for p in platforms_of(_l) if p != "SAN"]
    _ill_table = f"{OUT}/illumina/{_l}/dada2/asv_table.tsv" if _l in ILL_MARKERS else []

    rule:
        # Split the authentic set by platform. A centroid belongs to every
        # platform/sample represented among the sequences it stands for.
        # Illumina sample membership comes from the DADA2 ASV table.
        name: f"subsets_{_l}"
        input:
            fasta=f"{OUT}/authenticate/{_l}/authentic.fasta",
            members=f"{OUT}/authenticate/{_l}/precluster_members.tsv",
            table=_ill_table
        output:
            subsets=[f"{OUT}/authenticate/{_l}/subset_{p}.fasta" for p in _plats],
            platforms=f"{OUT}/authenticate/{_l}/seq_platforms.tsv"
        params:
            outdir=f"{OUT}/authenticate/{_l}",
            table_arg=f"--asv-table {_ill_table}" if _ill_table else "",
            expected=" ".join(_plats)
        conda:
            "../envs/pytools.yaml"
        shell:
            "python {SCRIPTS}/make_subsets.py --fasta {input.fasta} "
            "--members {input.members} {params.table_arg} --expected {params.expected} "
            "--outdir {params.outdir} --platforms-out {output.platforms}"
