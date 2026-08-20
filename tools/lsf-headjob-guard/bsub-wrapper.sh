#!/bin/bash
# bsub-wrapper.sh
#
# Template for the LSF head-job preemption guard. install.sh fills in
# __REAL_BSUB__ with the path to the real bsub binary and drops the result
# at ~/bin/bsub, ahead of the real one in PATH.
#
# What it does:
#   - Every bsub call whose piped-in script is tagged `#BSUB -J nf-workflow-*`
#     (Seqera Platform's own naming for a pipeline's Nextflow HEAD process)
#     is NOT handed to LSF. Instead it's written out and run directly as a
#     detached background process (setsid/nohup), so it can never be
#     preempted by LSF in the first place.
#   - Every other bsub call - i.e. every real pipeline task submitted by
#     Nextflow's own LSF executor (PREPROCESS, RUN_STATS, ...) - passes
#     straight through to the real bsub, completely untouched. Normal
#     LSF scheduling for actual pipeline work is unaffected.
#
# See ../README.md for why this exists. Cancellation of these detached
# head jobs from Seqera Platform is handled by the companion
# bkill-wrapper.sh, which reads the registry file this script writes to
# below.

REAL_BSUB="__REAL_BSUB__"
LOGFILE="$HOME/.bsub_wrapper/wrapper.log"
REGISTRY="$HOME/.bsub_wrapper/intercepted_jobs"
mkdir -p "$HOME/.bsub_wrapper/scripts"
touch "$REGISTRY"

ts() { date '+%Y-%m-%d %H:%M:%S'; }

if [ -t 0 ]; then
  # No script piped on stdin - not the pattern we care about, pass through.
  echo "$(ts) PASSTHROUGH(no-stdin) args: $*" >> "$LOGFILE"
  exec "$REAL_BSUB" "$@"
fi

script_content="$(cat)"

if printf '%s\n' "$script_content" | grep -q '^#BSUB -J nf-workflow-'; then
  jobname=$(printf '%s\n' "$script_content" | sed -n 's/^#BSUB -J \(.*\)/\1/p' | head -1)
  echo "$(ts) INTERCEPTED head job '$jobname'" >> "$LOGFILE"

  runfile="$HOME/.bsub_wrapper/scripts/${jobname}.$$.sh"
  printf '%s\n' "$script_content" > "$runfile"
  chmod +x "$runfile"

  setsid nohup bash "$runfile" </dev/null >/dev/null 2>&1 &
  disown
  pid=$!
  echo "$(ts) -> detached PID $pid, script $runfile" >> "$LOGFILE"

  # Record pid -> jobname so bkill-wrapper.sh can find and kill this
  # process when Seqera cancels the run by the fake job ID below (setsid
  # made $pid a process group leader, so bkill-wrapper can signal the
  # whole group).
  echo "$pid $jobname" >> "$REGISTRY"

  # Fake bsub's normal "submitted" line so anything parsing our stdout
  # (e.g. Seqera Platform / tw-agent) doesn't choke on unexpected output.
  echo "Job <$pid> is submitted to default queue <medium>."
  exit 0
else
  echo "$(ts) PASSTHROUGH(task) args: $*" >> "$LOGFILE"
  printf '%s\n' "$script_content" | "$REAL_BSUB" "$@"
fi
