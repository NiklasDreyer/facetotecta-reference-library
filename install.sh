#!/usr/bin/env bash
# =============================================================================
# install.sh: installs everything the pipeline needs.
#
# Run once, on a LOGIN node (compute nodes often have no internet access):
#     bash install.sh
#
# Steps:
#   1. Find conda/mamba (or tell you how to install Miniforge).
#   2. Create the "yfacet-snakemake" environment (Snakemake 8 + SLURM plugin).
#   3. Pre-build every per-rule conda environment (workflow/envs/*.yaml), so
#      SLURM jobs never have to download anything.
#   4. Fetch the metaMATE default specifications file.
#   5. Print the steps for installing ASAP, which is not available on conda.
#   6. Export exact package versions of every environment to envs_lock/ (for the
#      Methods section and for reproducibility).
# =============================================================================

# Stop at the first error, on unset variables, and on failures inside pipes.
set -euo pipefail

# Work from the folder that contains this script (the repository root).
cd "$(dirname "$0")"

# -----------------------------------------------------------------------------
# 1. Find a conda front-end. mamba is much faster at solving environments.
# -----------------------------------------------------------------------------
if command -v mamba >/dev/null 2>&1; then
    CONDA=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA=conda
else
    # No conda found: explain how to get Miniforge (conda-forge's installer).
    echo "No conda or mamba found. Install Miniforge first, e.g.:"
    echo "  curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    echo "  bash Miniforge3-Linux-x86_64.sh -b -p \$HOME/miniforge3"
    echo "  source \$HOME/miniforge3/etc/profile.d/conda.sh"
    echo "Many clusters also provide it as a module (e.g. 'module load miniforge3')."
    exit 1
fi
echo "Using $CONDA"

# Strict channel priority avoids mixing incompatible builds.
conda config --set channel_priority strict

# -----------------------------------------------------------------------------
# 2. Create (or update) the environment that runs Snakemake.
# -----------------------------------------------------------------------------
if conda env list | grep -q "^yfacet-snakemake "; then
    $CONDA env update -n yfacet-snakemake -f workflow/envs/snakemake.yaml
else
    $CONDA env create -f workflow/envs/snakemake.yaml
fi

# Make "conda activate" work inside a non-interactive script.
eval "$(conda shell.bash hook)"
conda activate yfacet-snakemake

# -----------------------------------------------------------------------------
# 3. Build every per-rule environment now. The workflow has one small rule per
#    environment (rules/common.smk, "env_<name>") whose only job is to trigger
#    the build. --conda-create-envs-only builds without running anything.
# -----------------------------------------------------------------------------
snakemake --snakefile workflow/Snakefile \
          --software-deployment-method conda \
          --conda-prefix .snakemake/conda \
          --conda-create-envs-only \
          --cores 1 \
          build_envs

# -----------------------------------------------------------------------------
# 4. metaMATE ships a default "specifications" file with the threshold grid it
#    tests. Copy it into resources/metamate/ so the run is documented.
# -----------------------------------------------------------------------------
if [ ! -s resources/metamate/specifications.txt ]; then
    tmp=$(mktemp -d)
    git clone --depth 1 https://github.com/tjcreedy/metamate.git "$tmp/metamate"
    # The file name has varied between metaMATE versions, so look for it.
    spec=$(find "$tmp/metamate" -iname "specifications*.txt" | head -n 1)
    if [ -n "$spec" ]; then
        cp "$spec" resources/metamate/specifications.txt
        echo "metaMATE specifications copied from $spec"
    else
        echo "WARNING: no specifications file found in the metaMATE repository."
        echo "         Write resources/metamate/specifications.txt by hand (see metaMATE README)."
    fi
    rm -rf "$tmp"
fi

# -----------------------------------------------------------------------------
# 5. ASAP (Puillandre et al. 2021) is not on conda. Install the command-line
#    version yourself and put the binary in bin/ (added to PATH by the workflow).
# -----------------------------------------------------------------------------
mkdir -p bin
if [ ! -x bin/asap ]; then
    echo ""
    echo "ASAP is not installed. To enable the ASAP rules:"
    echo "  1. Download the command-line source from https://bioinfo.mnhn.fr/abi/public/asap/"
    echo "  2. Compile it (usually 'make') and copy the 'asap' binary to $(pwd)/bin/asap"
    echo "  3. Run 'bin/asap -h' and check the flags against delimitation.asap.command in config/config.yaml"
    echo "  4. Set delimitation.asap.enabled: true"
    echo "Until then, ASAP can still be run on the web server with the alignments in results/delim/*/aln.fasta."
fi

# -----------------------------------------------------------------------------
# 6. Record exact package versions of every environment that was built.
#    Snakemake stores each environment under .snakemake/conda/<hash>/ with a
#    copy of its YAML as .snakemake/conda/<hash>.yaml.
# -----------------------------------------------------------------------------
mkdir -p envs_lock
for y in .snakemake/conda/*.yaml; do
    [ -e "$y" ] || continue
    h=$(basename "$y" .yaml)
    # Use the first comment line of the YAML to label the lock file.
    label=$(grep -m1 '^#' "$y" | tr -cd '[:alnum:] ' | tr ' ' '_' | cut -c1-40)
    conda env export -p ".snakemake/conda/$h" --no-builds > "envs_lock/${h}_${label}.yaml" || true
done

echo ""
echo "Installation finished. Next steps:"
echo "  conda activate yfacet-snakemake"
echo "  edit config/config.yaml, config/*.tsv and profiles/slurm/config.yaml"
echo "  snakemake --profile profiles/slurm -n        # dry run"
echo "  snakemake --profile profiles/slurm           # full run"
