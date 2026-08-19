# LSF head-job preemption guard

Works around a Seqera Platform + WEXAC LSF interaction that kills pipeline
runs outright, rather than just retrying a task.

## The problem

WEXAC's `medium`/`short` LSF queues preempt jobs often, under load. LSF's
`Re-runnable` preemption is normally transparent: a preempted job is killed
and automatically resubmitted under the same job ID, into the same work
directory, and just picks back up.

Seqera Platform submits the pipeline's **own Nextflow head process** as an
LSF job too (`bsub`'d, named `nf-workflow-<sessionID>`) - and that job is
just as preemptable as any other. The head job's launch script
(`nf-<id>.launcher.sh`) does a one-time `pre_run()` step near the top that
exports short-lived, single-use credentials (`TOWER_ACCESS_TOKEN`, and a
`NXF_SCM_FILE=https://api.cloud.seqera.io/ephemeral/<token>` URL) before
handing off to `nextflow run ...`. When LSF preempts and auto-resumes this
job, it re-executes the *entire script from scratch* - including
`pre_run()` - which re-exports credentials that were already spent by the
first attempt. The result is an immediate `403 Unable to access config
file ... ` and the whole Seqera-tracked run dies, even though nothing
about the actual pipeline failed.

Task-level jobs (`PREPROCESS`, `RUN_STATS`, etc.) don't have this problem -
Nextflow's own retry logic handles their preemption fine (see
`nxf_patch_130` in `nextflow.config` for a separate, unrelated fix for a
duplicate-execution issue on *those*). It's specifically the Platform's own
head job that this guard targets.

## What was tried and rejected

Catching Seqera's `nf-*.launcher.sh` the moment it appears, `bkill`-ing the
LSF job Platform just submitted for it, and running the launcher script
manually as a plain background process does fix the compute: the pipeline
runs fine, no duplicated work. But by the time this can be detected and
acted on, Nextflow has usually already made its first "workflow started"
call to Seqera's Tower API under that run's ID - so killing the LSF job
(even intentionally, even fast) reports that specific execution attempt as
dead to Seqera. A second process reporting under the same ID afterwards is
invisible in the Platform UI, even though it's doing the real work. This
is a race that can't be reliably won, and even when "won" it breaks
dashboard visibility.

## The actual fix

Seqera's `tw-agent` (Tower Agent, the process that connects this compute
environment to Seqera Cloud) submits the head job by piping a script into
`bsub` - the same script later found at `nf-<id>.launcher.sh`, tagged with
`#BSUB -J nf-workflow-<sessionID>`.

`bsub-wrapper.sh` is installed as `~/bin/bsub`, ahead of the real `bsub` in
`tw-agent`'s `PATH`. It inspects every script piped into it:

- If it's tagged `#BSUB -J nf-workflow-*` (a head job), the wrapper does
  **not** hand it to LSF. It writes the script out and runs it directly as
  a detached background process (`setsid`/`nohup`), immune to LSF
  preemption from the start - no race, no catching anything after the
  fact. It also prints a fake `Job <PID> is submitted...` line so nothing
  parsing `bsub`'s stdout trips over unexpected output.
- Anything else (every real pipeline task) is piped straight through to
  the real `bsub`, completely unmodified. Normal LSF scheduling for actual
  compute is unaffected.

`tw-agent` is long-running, though: it reads its environment once, at
startup, and every `bsub` it spawns for the rest of its life inherits that
same `PATH`. A plain `.bashrc`/`module load` fix doesn't reach an
already-running process, and (in practice) module loads run after
`.bashrc` tend to push `~/bin` back down the `PATH` anyway. So `tw-agent`
needs to be restarted with `~/bin` forced to the front of `PATH` at the
moment it starts - see `INSTALL.md` for the exact command.

## Install

See `INSTALL.md`.

## Verifying it's working

```bash
which bsub                          # should print ~/bin/bsub
tail -f ~/.bsub_wrapper/wrapper.log
```

Trigger a launch or a Resume in Seqera Platform. You should see a line
like:

```
2026-08-19 12:30:01 INTERCEPTED head job 'nf-workflow-<sessionID>'
```

Every real task submission logs a `PASSTHROUGH(task)` line instead, and
should keep showing up normally in `bjobs`.

## Known limitation

This relies on `tw-agent` resolving `bsub` via a plain `PATH` lookup. If
Seqera ever changes `tw-agent` to call an absolute path, or bundles its own
`bsub` resolution logic, this stops applying and the wrapper is simply
never hit (harmless - `bsub` for everything else keeps working normally
either way, it just stops intercepting head jobs). There's no way to
verify this ahead of time other than triggering a real run and checking
`wrapper.log`, as above.
