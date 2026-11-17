# =============================================================================
# common.smk: shared settings, sample tables and helper functions.
# Included first by both workflow/Snakefile and workflow/Snakefile_mining.
# =============================================================================

import os
import pandas as pd

# Shortcut for the results folder and the scripts folder.
OUT = config["outdir"]
SCRIPTS = os.path.join(workflow.basedir, "scripts")

# Put the repository's bin/ (for the ASAP binary) on PATH for every job.
os.environ["PATH"] = os.path.abspath("bin") + os.pathsep + os.environ["PATH"]

# Wildcards for loci/markers are simple names (COX1, 18S, ...): no dots or
# slashes, so file patterns cannot match each other by accident.
wildcard_constraints:
    locus="[A-Za-z0-9]+",
    marker="[A-Za-z0-9]+",
    unit="[^/]+"

# -----------------------------------------------------------------------------
# Illumina libraries: one "unit" per sample x marker, e.g. "Y1__COX1".
# -----------------------------------------------------------------------------
UNITS = pd.read_csv(config["illumina_samples"], sep="\t", comment="#", dtype=str)
UNITS["unit"] = UNITS["sample"] + "__" + UNITS["marker"]
UNITS = UNITS.set_index("unit", drop=False)

# Markers with Illumina data (e.g. COX1 and 18S).
ILL_MARKERS = sorted(UNITS["marker"].unique())


def units_of(marker):
    """Return the unit IDs (sample__marker) for one marker."""
    return UNITS.loc[UNITS["marker"] == marker, "unit"].tolist()


def samples_of(marker):
    """Return the Illumina sample names (e.g. Y1, Y2) for one marker."""
    return UNITS.loc[UNITS["marker"] == marker, "sample"].tolist()


# -----------------------------------------------------------------------------
# Loci
# -----------------------------------------------------------------------------
SANGER_LOCI = config["sanger"]["loci"]
CODING = config["sanger"]["coding"]          # e.g. {"COX1": 5, "H3": 1}
HTS_LOCI = config["hts_loci"]                # loci with Illumina and/or PacBio data
PB_BLAST_LOCI = config["pacbio"]["blast_loci"]  # PacBio loci recovered by BLAST


def pacbio_has(locus):
    """True if PacBio reads are extracted for this locus."""
    return locus == "COX1" or locus in PB_BLAST_LOCI


def platforms_of(locus):
    """Platform tags present for a locus: SAN, PB, and ILL_<sample>."""
    tags = ["SAN"]
    if pacbio_has(locus):
        tags.append("PB")
    if locus in ILL_MARKERS:
        tags += [f"ILL_{s}" for s in samples_of(locus)]
    return tags


# -----------------------------------------------------------------------------
# Datasets used for species delimitation. Each has a FASTA (unaligned unless
# "aligned" is True) and a locus, which sets the VSEARCH threshold.
# -----------------------------------------------------------------------------
DELIM = {}
for _l in SANGER_LOCI:
    DELIM[f"sanger_{_l}"] = {"locus": _l, "aligned": False,
                             "fasta": f"{OUT}/sanger/{_l}/{_l}.clean.fasta"}
DELIM["sanger_concat"] = {"locus": "concat", "aligned": True,
                          "fasta": f"{OUT}/sanger/concat/concat.fasta"}
for _l in HTS_LOCI:
    # All authenticated sequences of a locus, from every platform.
    DELIM[f"hts_{_l}_all"] = {"locus": _l, "aligned": False,
                              "fasta": f"{OUT}/authenticate/{_l}/authentic.fasta"}
    # One subset per high-throughput platform (for the per-platform counts).
    for _p in platforms_of(_l):
        if _p != "SAN":
            DELIM[f"hts_{_l}_{_p}"] = {"locus": _l, "aligned": False,
                                       "fasta": f"{OUT}/authenticate/{_l}/subset_{_p}.fasta"}

DELIM_NAMES = list(DELIM.keys())

# Delimitation methods that produce a membership table for every dataset.
DELIM_METHODS = ["vsearch", "mptp", "ptp"]
if config["delimitation"]["asap"]["enabled"]:
    DELIM_METHODS.append("asap")


def vsearch_id(dataset):
    """VSEARCH identity threshold for a dataset, taken from its locus."""
    return config["delimitation"]["vsearch_id"][DELIM[dataset]["locus"]]


# -----------------------------------------------------------------------------
# Environments: one tiny rule per environment so install.sh can build them all
# with "snakemake --conda-create-envs-only build_envs".
# -----------------------------------------------------------------------------
ENV_NAMES = ["qc", "trim", "dada2", "blast", "align", "authenticate", "pytools",
             "vsearch", "transdecoder", "metamate", "mapping", "mptp", "inext",
             "placement", "sra", "cutadapt", "benchmark"]
if config["pacbio"]["deepconsensus"]["enabled"]:
    ENV_NAMES.append("deepconsensus")

for _env in ENV_NAMES:
    rule:
        name: f"env_{_env}"
        output: f"{OUT}/.envs/{_env}.ok"
        conda: f"../envs/{_env}.yaml"
        shell: "touch {output}"

rule build_envs:
    # Target used by install.sh: depends on one marker file per environment.
    input: expand(f"{OUT}/.envs/{{env}}.ok", env=ENV_NAMES)
