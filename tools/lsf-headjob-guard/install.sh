#!/bin/bash
# install.sh - installs the LSF head-job preemption guard (see README.md).
# Usage: ./install.sh
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REAL_BSUB="$(type -ap bsub | grep -v "^$HOME/" | head -1)"
if [[ -z "$REAL_BSUB" ]]; then
  echo "ERROR: no real bsub found in PATH (run this where LSF is loaded)." >&2
  exit 1
fi

mkdir -p "$HOME/bin" "$HOME/.bsub_wrapper/scripts"
sed "s|__REAL_BSUB__|$REAL_BSUB|" "$DIR/bsub-wrapper.sh" > "$HOME/bin/bsub"
chmod +x "$HOME/bin/bsub"

echo "Installed ~/bin/bsub (real bsub: $REAL_BSUB)"
echo "Now restart tw-agent so it picks it up - see INSTALL.md."
