#!/bin/bash
# bkill-wrapper.sh
#
# Companion to bsub-wrapper.sh. install.sh fills in __REAL_BKILL__ with the
# path to the real bkill binary and drops the result at ~/bin/bkill, ahead
# of the real one in PATH.
#
# Why this exists:
#   bsub-wrapper.sh intercepts Seqera's Nextflow HEAD job submissions
#   (`#BSUB -J nf-workflow-*`) and runs them as a detached background
#   process instead of a real LSF job, printing a fake
#   "Job <PID> is submitted..." line so tw-agent/Seqera Platform treats
#   that PID as the LSF job ID for the run. It records `<pid> <jobname>`
#   in ~/.bsub_wrapper/intercepted_jobs.
#
#   When someone hits "Cancel" on that run in Seqera Platform, tw-agent
#   calls `bkill <that same fake job ID>`. The real bkill would just fail
#   (WEXAC LSF has no job with that ID) - so before this wrapper existed,
#   a head job could never be cancelled from Seqera; it kept running until
#   it finished on its own or someone killed the PID by hand on the host.
#
# What it does:
#   - If any argument looks like one of our fake job IDs (a bare number
#     matching a PID in the registry) or matches an intercepted job's name
#     (covers `bkill -J <jobname>`), it signals that process group
#     directly (SIGTERM, then SIGKILL if it's still alive after a couple
#     seconds), prints a bkill-style "is being terminated" line, and never
#     touches the real bkill.
#   - Anything else - i.e. every real pipeline task LSF job - passes
#     straight through to the real bkill, completely untouched. Normal
#     LSF job cancellation is unaffected.
#
# See ../README.md for why this exists.

REAL_BKILL="__REAL_BKILL__"
LOGFILE="$HOME/.bsub_wrapper/wrapper.log"
REGISTRY="$HOME/.bsub_wrapper/intercepted_jobs"
mkdir -p "$HOME/.bsub_wrapper"
touch "$REGISTRY"

ts() { date '+%Y-%m-%d %H:%M:%S'; }

# Look for a match against our registry ("<pid> <jobname>" per line): try
# each argument first as a fake job ID (bare PID), then as a job name (for
# `bkill -J <jobname>` style calls). Last match wins, same "don't overthink
# LSF's flag grammar" approach as bsub-wrapper.sh's stdin sniff.
match=""
for arg in "$@"; do
  case "$arg" in
    ''|*[!0-9]*) ;;  # not purely numeric - can't be a fake job ID
    *) line="$(grep "^${arg} " "$REGISTRY" 2>/dev/null | tail -1)"
       [ -n "$line" ] && match="$line" ;;
  esac
done
if [ -z "$match" ]; then
  for arg in "$@"; do
    line="$(grep " ${arg}\$" "$REGISTRY" 2>/dev/null | tail -1)"
    [ -n "$line" ] && match="$line"
  done
fi

if [ -n "$match" ]; then
  pid="$(printf '%s' "$match" | awk '{print $1}')"
  jobname="$(printf '%s' "$match" | cut -d' ' -f2-)"

  if kill -0 "$pid" 2>/dev/null; then
    echo "$(ts) INTERCEPTED kill of head job '$jobname' (fake job <$pid>)" >> "$LOGFILE"
    # setsid made $pid its own process group leader in bsub-wrapper.sh, so
    # "-$pid" signals the whole group (nextflow + whatever it spawned).
    # Fall back to a plain kill if the group signal fails for some reason.
    kill -TERM "-$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null
    sleep 2
    if kill -0 "$pid" 2>/dev/null; then
      kill -KILL "-$pid" 2>/dev/null || kill -KILL "$pid" 2>/dev/null
      echo "$(ts) -> PID $pid didn't stop on SIGTERM, sent SIGKILL" >> "$LOGFILE"
    else
      echo "$(ts) -> PID $pid stopped" >> "$LOGFILE"
    fi
  else
    echo "$(ts) INTERCEPTED kill of head job '$jobname' (fake job <$pid>) but it was already gone" >> "$LOGFILE"
  fi

  # Drop this entry so it isn't matched again by a future bkill/bsub. grep
  # -v exits 1 (not an error here) when that leaves nothing behind, e.g.
  # this was the only entry - so don't gate the mv on grep's exit status.
  grep -v "^${pid} " "$REGISTRY" > "$REGISTRY.tmp" 2>/dev/null
  mv "$REGISTRY.tmp" "$REGISTRY"

  # Mimic real bkill's normal confirmation line so anything parsing our
  # stdout (tw-agent / Seqera Platform) doesn't choke on unexpected output.
  echo "Job <$pid> is being terminated"
  exit 0
else
  echo "$(ts) PASSTHROUGH(bkill) args: $*" >> "$LOGFILE"
  exec "$REAL_BKILL" "$@"
fi
