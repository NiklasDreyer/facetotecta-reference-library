# =============================================================================
# placement.smk: phylogenetic placement on a Facetotecta + outgroup reference tree.
# Used by the benchmark (non-target crustaceans) and by the public-data mining
# (candidate facetotectan ASVs).
#
#   reference: Sanger library + outgroups -> MAFFT -> IQ-TREE (GTR+G) ->
#              RAxML-NG model evaluation (EPA-ng needs a RAxML-NG model file)
#   queries:   MAFFT --addfragments --keeplength -> EPA-ng -> gappa assign
# =============================================================================

P = f"{OUT}/placement"


rule placement_reference:
    # Reference alignment, tree, model file and taxonomy for one marker.
    input:
        sanger=f"{OUT}/sanger/{{marker}}/{{marker}}.clean.fasta",
        outgroups=lambda w: os.path.join(config["authenticate"]["outgroups_dir"], f"{w.marker}.fasta")
    output:
        aln=f"{P}/{{marker}}/reference/ref.aln.fasta",
        tree=f"{P}/{{marker}}/reference/ref.raxml.bestTree",
        model=f"{P}/{{marker}}/reference/ref.raxml.bestModel",
        tax=f"{P}/{{marker}}/reference/ref_taxonomy.tsv"
    params:
        d=f"{P}/{{marker}}/reference"
    threads: 8
    log:
        f"{OUT}/logs/placement/{{marker}}.reference.log"
    conda:
        "../envs/placement.yaml"
    shell:
        "cat {input.sanger} {input.outgroups} > {params.d}/ref.fasta; "
        "mafft --auto --maxiterate 1000 --thread {threads} {params.d}/ref.fasta > {output.aln} 2> {log}; "
        "iqtree2 -s {output.aln} -m GTR+G -T {threads} --prefix {params.d}/ref -redo >> {log} 2>&1; "
        # RAxML-NG re-optimises the model on the IQ-TREE topology; EPA-ng needs
        # its model file, and its tree file has no support labels.
        "raxml-ng --evaluate --msa {output.aln} --tree {params.d}/ref.treefile --model GTR+G "
        "--threads {threads} --prefix {params.d}/ref --redo >> {log} 2>&1; "
        # Taxonomy for gappa: Sanger sequences = Facetotecta, the rest = Outgroup.
        "grep '^>' {output.aln} | sed 's/^>//' | awk '{{split($1,a,\"-\"); "
        "if (a[1]==\"SAN\") print $1\"\\tFacetotecta\"; else print $1\"\\tOutgroup\"}}' > {output.tax}"


rule place_queries:
    # Place query sequences ({qset}/query.fasta is written by another rule).
    input:
        query=f"{P}/{{marker}}/{{qset}}/query.fasta",
        aln=f"{P}/{{marker}}/reference/ref.aln.fasta",
        tree=f"{P}/{{marker}}/reference/ref.raxml.bestTree",
        model=f"{P}/{{marker}}/reference/ref.raxml.bestModel",
        tax=f"{P}/{{marker}}/reference/ref_taxonomy.tsv"
    output:
        jplace=f"{P}/{{marker}}/{{qset}}/epa_result.jplace",
        assign=f"{P}/{{marker}}/{{qset}}/assignments.tsv"
    params:
        d=f"{P}/{{marker}}/{{qset}}"
    threads: 8
    log:
        f"{OUT}/logs/placement/{{marker}}.{{qset}}.log"
    conda:
        "../envs/placement.yaml"
    shell:
        # The script adds the queries to the fixed reference alignment (MAFFT
        # --addfragments --keeplength), splits query and reference alignments,
        # runs EPA-ng and assigns each query with gappa examine assign.
        "python {SCRIPTS}/run_placement.py --query {input.query} --ref-aln {input.aln} "
        "--tree {input.tree} --model {input.model} --tax {input.tax} "
        "--workdir {params.d} --threads {threads} --output {output.assign} > {log} 2>&1"
