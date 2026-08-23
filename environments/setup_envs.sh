#!/usr/bin/env bash
# One-time setup for the pipeline's conda envs.
#
# main.nf's RUN_LUMPY, RUN_STATS, RUN_AMPLICONARCHITECT and ARCHIVE
# processes reference these envs by fixed path (params.*_env in nextflow.config,
# all under params.conda_envs = "$HOME/.conda/envs"). Build them once here.
# Rerun a single env's line below if you update that env's .yml.
set -euo pipefail

ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENVS="${HOME}/.conda/envs"

# Why specific miniconda?
# WEXAC's default `mamba`/`conda` module has a solver old enough to report
# false "nothing provides X" errors for packages that are perfectly real
# (mscorefonts, pomegranate, biopython all hit this) - load the modern one.
module load miniconda/26.1.1_environmentally

# Why ToS?
# Accepts Anaconda's ToS for the `defaults` channels - safe no-op once
# accepted. WEXAC's conda module has these channels active by default, so
# `conda env create` checks their ToS every run (unlike `conda create`,
# it has no --override-channels to skip this) and hangs on a prompt
# otherwise, even though none of our envs actually use `defaults`.
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Why mosek?
# ampliconsuite 1.6.0 unconditionally imports mosek at module load
# (bam_to_breakpoint.py), so RUN_AMPLICONARCHITECT fails with
# ModuleNotFoundError without it - even though the actual solver used is
# the free bundled Clarabel. The free/trial license at ~/mosek/mosek.lic
# is enough to satisfy the import. Revisit if a future release makes it lazy.

for name in lumpy-sv bw ampsuite; do
  target="${CONDA_ENVS}/${name}"
  if [ -d "${target}" ]; then
    echo "Skipping ${name}: ${target} already exists"
    continue
  fi
  echo "Creating ${name} at ${target}"
  conda env create -f "${ENV_DIR}/${name}.yml" -p "${target}"
done

echo "Done. Envs live under ${CONDA_ENVS}, matching params.conda_envs in nextflow.config."
