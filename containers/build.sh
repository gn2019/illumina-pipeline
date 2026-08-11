#!/usr/bin/env bash
# One-time build step for the cgpwgs Singularity image.
#
# This is intentionally NOT a Nextflow process. Building a Singularity image
# from a docker:// source needs privileges (root or fakeroot) that most
# cluster/HPC compute nodes don't grant to pipeline jobs, and running it
# inside the DAG risks concurrent runs racing on the same cache path.
#
# Build once, store the resulting .sif somewhere shared/versioned, and point
# nextflow.config's params.cgpwgs_sif at it. Rerun this script only when you
# need to bump the cgpwgs version.
set -euo pipefail

IMAGE_TAG="${1:-2.1.1}"

EXPECTED_BASE="$HOME/illumina-pipeline"
OUT_FILE="${EXPECTED_BASE}/cgpwgs.sif"

singularity build \
  "${OUT_FILE}" \
  "docker://quay.io/wtsicgp/dockstore-cgpwgs:${IMAGE_TAG}"

echo "Built ${OUT_FILE} from dockstore-cgpwgs:${IMAGE_TAG}"
