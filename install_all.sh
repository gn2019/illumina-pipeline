#!/bin/bash
# install_all.sh - runs every installer in this repo: directory layout +
# scripts (install_paths.sh), conda environments, the cgpwgs container,
# and any tools/*/install.sh (e.g. the WEXAC/Seqera LSF head-job guard).
#
# Usage: ./install_all.sh [base_dir]
#   base_dir is passed through to install_paths.sh - see there for its
#   default.
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$DIR/install_paths.sh" "$@"
bash "$DIR/environments/setup_envs.sh"
bash "$DIR/containers/build.sh"
for installer in "$DIR"/tools/*/install.sh; do
  [[ -f "$installer" ]] && bash "$installer"
done

echo "All installers finished."
