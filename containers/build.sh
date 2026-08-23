#!/usr/bin/env bash
# One-time build step for the cgpwgs Singularity image.
#
# Build once, and point nextflow.config's params.cgpwgs_sif at it.
# Rerun this script only when you need to bump the cgpwgs version.
set -euo pipefail

IMAGE_TAG="${1:-2.1.1}"

EXPECTED_BASE="$HOME/illumina-pipeline"
OUT_FILE="${EXPECTED_BASE}/cgpwgs.sif"

singularity build \
  "${OUT_FILE}" \
  "docker://quay.io/wtsicgp/dockstore-cgpwgs:${IMAGE_TAG}"

echo "Built ${OUT_FILE} from dockstore-cgpwgs:${IMAGE_TAG}"
