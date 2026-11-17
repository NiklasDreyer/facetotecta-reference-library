# =============================================================================
# delimitation.smk: MOTU delimitation (Section 2.10).
#
# Runs on every dataset in DELIM (common.smk): each Sanger locus, the
# concatenated Sanger matrix, all authentic HTS sequences per locus, and the
# per-platform HTS subsets. For each dataset:
#   - alignment (MAFFT) and ungapped sequences
#   - ML tree (IQ-TREE, ModelFinder, UFBoot -bnni, SH-aLRT), midpoint-rooted
#   - VSEARCH clustering at the locus threshold (manuscript_v1)
#   - PTP (single-rate, the model behind bPTP) and mPTP (multi-rate), ML
#     partitions plus MCMC support
#   - ASAP (optional, needs the ASAP binary)
#   - NEW: threshold sweep 90-99% (Fig. S8)
# Then NEW: congruence between methods and loci, and MOTU grading.
# Every method writes a two-column membership table: sequence <tab> MOTU.
# =============================================================================

D = f"{OUT}/delim"

wildcard_constraints:
    ds="|".join(DELIM_NAMES)


rule delim_prepare:
    # Make an alignment (MAFFT --auto --maxiterate 1000) and an ungapped copy.
    # The concatenated matrix is already aligned and is only copied.
    input:
        lambda w: DELIM[w.ds]["fasta"]
    output:
        aln=f"{D}/{{ds}}/aln.fasta",
        seqs=f"{D}/{{ds}}/seqs.fasta"
    params:
        aligned=lambda w: "yes" if DELIM[w.ds]["aligned"] else "no"
    threads: 4
    log:
        f"{OUT}/logs/delim/{{ds}}.prepare.log"
    conda:
        "../envs/align.yaml"
    shell:
        # An empty subset (e.g. no sequences from one platform) gives empty outputs.
        "if [ ! -s {input} ]; then touch {output.aln} {output.seqs}; exit 0; fi; "
        "if [ {params.aligned} = yes ]; then cp {input} {output.aln}; "
        "else mafft --auto --maxiterate 1000 --thread {threads} {input} > {output.aln} 2> {log}; fi; "
        # Ungapped sequences for VSEARCH (removes "-" and "?").
        "seqkit seq -g {output.aln} | sed '/^>/! s/[?]//g' > {output.seqs}"


rule delim_iqtree:
    # Single-locus (or concatenated, unpartitioned) ML tree used by PTP/mPTP.
    input:
        f"{D}/{{ds}}/aln.fasta"
    output:
        tree=f"{D}/{{ds}}/iqtree/tree.treefile"
    params:
        prefix=f"{D}/{{ds}}/iqtree/tree",
        ufboot=config["delimitation"]["ufboot"],
        alrt=config["delimitation"]["alrt"]
    threads: 8
    log:
        f"{OUT}/logs/delim/{{ds}}.iqtree.log"
    conda:
        "../envs/align.yaml"
    shell:
        # Fewer than 4 sequences: no tree can be inferred; write an empty tree
        # so tree-based methods report nothing for this dataset.
        "n=$(grep -c '^>' {input} || true); "
        "if [ \"$n\" -lt 4 ]; then mkdir -p $(dirname {output.tree}); echo '();' > {output.tree}; exit 0; fi; "
        "iqtree2 -s {input} -m MFP -B {params.ufboot} -alrt {params.alrt} -bnni "
        "-T {threads} --prefix {params.prefix} -redo > {log} 2>&1"


rule delim_root:
    # Midpoint-root the tree (PTP and mPTP need a rooted tree).
    input:
        f"{D}/{{ds}}/iqtree/tree.treefile"
    output:
        f"{D}/{{ds}}/iqtree/tree.rooted.nwk"
    conda:
        "../envs/mptp.yaml"
    shell:
        "python {SCRIPTS}/root_tree.py --input {input} --output {output}"


rule delim_ptp_ml:
    # Maximum-likelihood delimitation with mPTP: "--multi" = mPTP (one
    # coalescent rate per species), "--single" = PTP (one rate for all species).
    input:
        tree=f"{D}/{{ds}}/iqtree/tree.rooted.nwk",
        aln=f"{D}/{{ds}}/aln.fasta"
    output:
        raw=f"{D}/{{ds}}/{{model}}/{{model}}.txt",
        members=f"{D}/{{ds}}/membership/{{model}}.tsv"
    params:
        flag=lambda w: "--multi" if w.model == "mptp" else "--single",
        prefix=f"{D}/{{ds}}/{{model}}/{{model}}"
    wildcard_constraints:
        model="mptp|ptp"
    log:
        f"{OUT}/logs/delim/{{ds}}.{{model}}.log"
    conda:
        "../envs/mptp.yaml"
    shell:
        "bash {SCRIPTS}/run_mptp.sh ml {input.tree} {input.aln} {params.flag} "
        "{params.prefix} > {log} 2>&1; "
        "python {SCRIPTS}/parse_mptp.py --input {output.raw} --output {output.members} "
        "--motu-prefix {wildcards.model}"


rule delim_mptp_mcmc:
    # MCMC runs for delimitation support (5 runs of 1e8 generations in the
    # default config). Both models are run.
    input:
        tree=f"{D}/{{ds}}/iqtree/tree.rooted.nwk",
        aln=f"{D}/{{ds}}/aln.fasta"
    output:
        done=touch(f"{D}/{{ds}}/{{model}}/mcmc.done")
    params:
        flag=lambda w: "--multi" if w.model == "mptp" else "--single",
        prefix=f"{D}/{{ds}}/{{model}}/mcmc",
        gens=config["delimitation"]["mptp_mcmc_generations"],
        runs=config["delimitation"]["mptp_mcmc_runs"],
        sample=config["delimitation"]["mptp_mcmc_sample"]
    wildcard_constraints:
        model="mptp|ptp"
    log:
        f"{OUT}/logs/delim/{{ds}}.{{model}}.mcmc.log"
    conda:
        "../envs/mptp.yaml"
    shell:
        "bash {SCRIPTS}/run_mptp.sh mcmc {input.tree} {input.aln} {params.flag} "
        "{params.prefix} {params.gens} {params.runs} {params.sample} > {log} 2>&1"


rule delim_vsearch:
    # VSEARCH --cluster_fast at the locus threshold (manuscript_v1: threshold
    # from ASAP). Centroids and members go into the membership table.
    input:
        f"{D}/{{ds}}/seqs.fasta"
    output:
        uc=f"{D}/{{ds}}/vsearch/clusters.uc",
        members=f"{D}/{{ds}}/membership/vsearch.tsv"
    params:
        ident=lambda w: vsearch_id(w.ds)
    threads: 2
    log:
        f"{OUT}/logs/delim/{{ds}}.vsearch.log"
    conda:
        "../envs/vsearch.yaml"
    shell:
        "if [ ! -s {input} ]; then touch {output.uc}; else "
        "vsearch --cluster_fast {input} --id {params.ident} --uc {output.uc} "
        "--threads {threads} > {log} 2>&1; fi; "
        "python {SCRIPTS}/parse_uc.py --uc {output.uc} --output {output.members}"


rule delim_asap:
    # ASAP on the alignment (only when delimitation.asap.enabled is true).
    input:
        f"{D}/{{ds}}/aln.fasta"
    output:
        members=f"{D}/{{ds}}/membership/asap.tsv"
    params:
        outdir=f"{D}/{{ds}}/asap",
        cmd=config["delimitation"]["asap"]["command"],
        glob=config["delimitation"]["asap"]["partition_glob"]
    log:
        f"{OUT}/logs/delim/{{ds}}.asap.log"
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/run_asap.py --fasta {input} --outdir {params.outdir} "
        "--command '{params.cmd}' --partition-glob '{params.glob}' "
        "--output {output.members} > {log} 2>&1"


rule delim_sweep:
    # NEW: number of VSEARCH clusters at every identity from 90% to 99%.
    input:
        f"{D}/{{ds}}/seqs.fasta"
    output:
        tsv=f"{D}/{{ds}}/sweep/threshold_sweep.tsv",
        pdf=f"{D}/{{ds}}/sweep/threshold_sweep.pdf"
    params:
        s=config["delimitation"]["sweep_start"],
        e=config["delimitation"]["sweep_end"],
        step=config["delimitation"]["sweep_step"],
        chosen=lambda w: vsearch_id(w.ds)
    threads: 2
    log:
        f"{OUT}/logs/delim/{{ds}}.sweep.log"
    conda:
        "../envs/vsearch.yaml"
    shell:
        "python {SCRIPTS}/threshold_sweep.py --input {input} --start {params.s} "
        "--end {params.e} --step {params.step} --chosen {params.chosen} "
        "--threads {threads} --out-tsv {output.tsv} --out-pdf {output.pdf} > {log} 2>&1"


rule delim_counts:
    # Number of MOTUs and singletons per method for one dataset.
    input:
        expand(f"{D}/{{{{ds}}}}/membership/{{method}}.tsv", method=DELIM_METHODS)
    output:
        f"{D}/{{ds}}/motu_counts.metrics.tsv"
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/count_motus.py --dataset {wildcards.ds} "
        "--memberships {input} --output {output}"


rule delim_congruence:
    # NEW: agreement between methods (within each dataset) and between COX1 and
    # 18S MOTUs (per specimen), with match ratio and adjusted Rand index; MOTUs
    # of the primary dataset are graded A/B/C.
    input:
        expand(f"{D}/{{ds}}/membership/{{method}}.tsv", ds=DELIM_NAMES, method=DELIM_METHODS)
    output:
        pairs=f"{D}/congruence/congruence_pairs.tsv",
        grades=f"{D}/congruence/motu_grades.tsv",
        metrics=f"{D}/congruence/congruence.metrics.tsv"
    params:
        delim=D,
        datasets=" ".join(DELIM_NAMES),
        methods=" ".join(DELIM_METHODS),
        primary=config["delimitation"]["grading_primary"],
        nuclear=config["delimitation"]["grading_nuclear"]
    conda:
        "../envs/pytools.yaml"
    shell:
        "python {SCRIPTS}/delim_congruence.py --delim-dir {params.delim} "
        "--datasets {params.datasets} --methods {params.methods} "
        "--primary {params.primary} --nuclear {params.nuclear} "
        "--out-pairs {output.pairs} --out-grades {output.grades} --metrics {output.metrics}"
