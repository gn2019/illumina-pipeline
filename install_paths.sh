#!/bin/bash
# install_paths.sh - creates the directory layout nextflow.config expects
# under params.base, and copies this repo's scripts/ into place there.
#
# Usage: ./install_paths.sh [base_dir]
#   base_dir defaults to $HOME/illumina-pipeline, matching params.base in
#   nextflow.config.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${1:-$HOME/illumina-pipeline}"

# Mirrors the params {} block in nextflow.config - update here if that changes.
mkdir -p \
  "$BASE/scripts" \
  "$BASE/results" \
  "$BASE/refs/data_repo" \
  "$BASE/data/fastq" \
  "$BASE/data/genome/hg38"

if [[ -d "$REPO_DIR/scripts" && "$REPO_DIR/scripts" != "$BASE/scripts" ]]; then
  cp -rn "$REPO_DIR"/scripts/. "$BASE/scripts/"
  echo "Copied $REPO_DIR/scripts -> $BASE/scripts"
fi

echo "Directory layout ready under $BASE"
echo "Still needed by hand: genome fasta under data/genome/hg38/, fastq"
echo "files under data/fastq/<sample>/, and the ASCAT GC-correction file"
echo "under refs/hg38/.../ascat/ - see nextflow.config for exact filenames."
