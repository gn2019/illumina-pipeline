#!/bin/bash
# install.sh - installs the LSF head-job preemption guard (see README.md)
# and its companion cancellation wrapper.
# Usage: ./install.sh
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REAL_BSUB="$(type -ap bsub | grep -v "^$HOME/" | head -1)"
if [[ -z "$REAL_BSUB" ]]; then
  echo "ERROR: no real bsub found in PATH (run this where LSF is loaded)." >&2
  exit 1
fi

REAL_BKILL="$(type -ap bkill | grep -v "^$HOME/" | head -1)"
if [[ -z "$REAL_BKILL" ]]; then
  echo "ERROR: no real bkill found in PATH (run this where LSF is loaded)." >&2
  exit 1
fi

mkdir -p "$HOME/bin" "$HOME/.bsub_wrapper/scripts"
touch "$HOME/.bsub_wrapper/intercepted_jobs"

sed "s|__REAL_BSUB__|$REAL_BSUB|" "$DIR/bsub-wrapper.sh" > "$HOME/bin/bsub"
chmod +x "$HOME/bin/bsub"

sed "s|__REAL_BKILL__|$REAL_BKILL|" "$DIR/bkill-wrapper.sh" > "$HOME/bin/bkill"
chmod +x "$HOME/bin/bkill"

echo "Installed ~/bin/bsub (real bsub: $REAL_BSUB)"
echo "Installed ~/bin/bkill (real bkill: $REAL_BKILL)"
echo "Now restart tw-agent so it picks them up - see INSTALL.md."
