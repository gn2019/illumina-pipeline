#!/bin/bash
# run_tw_agent.sh - runs tw-agent with ~/bin (the bsub wrapper) forced to
# the front of PATH, restarting it automatically if it exits.
#
# Two things this handles:
#   - tw-agent only picks up ~/bin/bsub if launched with it already ahead
#     of the real bsub in PATH (see ../README.md) - this sets that up.
#   - tw-agent has a known upstream bug (tower-agent#57, still open as of
#     v0.5.5) where it sometimes drops its connection and then fails to
#     reconnect on its own ("There is an active agent for this user and
#     connection ID"). This restarts it when that happens instead of
#     leaving it dead.
#
# Usage: ./run_tw_agent.sh <ENV_NAME> --access-token=<TOKEN>
# Run this (inside tmux/screen) instead of ./tw-agent directly, both for
# the first start and any later restart. Ctrl-C stops it for good.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
export PATH="$HOME/bin:$PATH"

module load Nextflow
trap 'echo stopping; exit' INT TERM

while true; do
  ./tw-agent "$@"
  echo "$(date '+%F %T') tw-agent exited ($?), restarting in 10s (Ctrl-C to stop)..." >&2
  sleep 10
done
