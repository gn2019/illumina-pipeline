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
# One-time-per-account, safe to run every time (no-op once accepted):
# accepts Anaconda's Terms of Service for the `defaults` channels
# (repo.anaconda.com/pkgs/main and /pkgs/r). WEXAC's conda module ships
# with these channels pre-configured as "active" for every account, so
# `conda env create` checks their ToS on every run regardless of what
# channels a given environment.yml lists - `conda env create`, unlike
# `conda create`, has no `--override-channels` flag to skip that check,
# so without this the loop below just hangs on an interactive prompt the
# first time any of these three envs gets created. None of lumpy-sv.yml,
# bw.yml or ampsuite.yml actually needs a package from `defaults`; this
# only silences the otherwise-hanging prompt.
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Why mosek?
# ampsuite.yml pulls in the `mosek` package (from the `mosek` channel) even
# though ampliconsuite's actual solver is the free bundled Clarabel. As of
# ampliconsuite 1.6.0, ampliconarchitectlib/bam_to_breakpoint.py does an
# unconditional `import mosek_solver` -> `import mosek` at module load time,
# so RUN_AMPLICONARCHITECT fails with ModuleNotFoundError before it ever
# gets to pick a solver if the mosek package isn't installed - regardless of
# which solver you actually intend to use. The existing free/trial license
# at ~/mosek/mosek.lic is enough to satisfy the import. Revisit this if a
# future ampliconsuite release makes that import lazy/optional.

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
