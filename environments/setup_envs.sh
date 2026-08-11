#!/usr/bin/env bash
# One-time setup for the pipeline's conda envs.
#
# Not a Nextflow process, by design: main.nf's RUN_LUMPY/RUN_STATS,
# RUN_AMPLICONARCHITECT and ARCHIVE processes reference these envs by
# fixed path (params.lumpy_env / params.ampsuite_env / params.bw_env in
# nextflow.config, all under params.conda_envs = "$HOME/.conda/envs").
# Build them once here so the paths already in nextflow.config resolve -
# no config changes needed. Rerun a single env's line below if you update
# that env's .yml.
set -euo pipefail

ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENVS="${HOME}/.conda/envs"

for name in lumpy-sv bw ampsuite; do
  target="${CONDA_ENVS}/${name}"
  if [ -d "${target}" ]; then
    echo "Skipping ${name}: ${target} already exists"
    continue
  fi
  echo "Creating ${name} at ${target}"
  conda env create -f "${ENV_DIR}/${name}.yml" -p "${target}"
done

# ampsuite (AmpliconSuite) needs a MOSEK license at runtime. This can't be
# fetched here - MOSEK only issues license files via their web form to a
# verified academic/institutional email, then emails you the .lic file. This
# just checks it's in the place MOSEK looks by default and points you at the
# request page if it's missing; nothing to configure in nextflow.config.
MOSEK_LIC="${HOME}/mosek/mosek.lic"
if [ -s "${MOSEK_LIC}" ]; then
  echo "Found MOSEK license at ${MOSEK_LIC}"
else
  echo
  echo "WARNING: no MOSEK license found at ${MOSEK_LIC}."
  echo "  RUN_AMPLICONARCHITECT (ampsuite env) will fail without one."
  echo "  Request a personal academic license (needs an academic email):"
  echo "    https://www.mosek.com/license/request/personal-academic/"
  echo "  Then save the emailed file to: ${MOSEK_LIC}"
fi

echo "Done. Envs live under ${CONDA_ENVS}, matching params.conda_envs in nextflow.config."
