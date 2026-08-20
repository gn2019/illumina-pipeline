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

## Cancelling a head job from Seqera

Running the head job as a detached process instead of a real LSF job has a
side effect: cancelling. Seqera Platform's "Cancel" action calls
`bkill <jobID>`, using the fake job ID (the process's PID) that
`bsub-wrapper.sh` printed at submission time. Without anything else in
place, real `bkill` just fails on that ID (WEXAC LSF never heard of it),
so there'd be no way to stop the run from Seqera - only by killing the PID
by hand on the compute host.

`bkill-wrapper.sh` is installed the same way as `bsub-wrapper.sh`: dropped
at `~/bin/bkill`, ahead of the real `bkill` in `tw-agent`'s `PATH`.
`bsub-wrapper.sh` records `<pid> <jobname>` for every head job it
intercepts in `~/.bsub_wrapper/intercepted_jobs`; when `bkill-wrapper.sh`
sees a call whose job ID (or `-J` job name) matches an entry there, it
signals that process group directly (`SIGTERM`, then `SIGKILL` if it's
still around after a couple seconds) and reports success in `bkill`'s
usual format, instead of forwarding to LSF. Every other `bkill` call - i.e.
cancelling a real pipeline task - passes straight through untouched.

## Install

See `INSTALL.md`.

## Verifying it's working

```bash
which bsub                          # should print ~/bin/bsub
which bkill                         # should print ~/bin/bkill
tail -f ~/.bsub_wrapper/wrapper.log
```

Trigger a launch or a Resume in Seqera Platform. You should see a line
like:

```
2026-08-19 12:30:01 INTERCEPTED head job 'nf-workflow-<sessionID>'
```

Every real task submission logs a `PASSTHROUGH(task)` line instead, and
should keep showing up normally in `bjobs`.

Then hit Cancel on that run in Seqera Platform and watch for:

```
2026-08-19 12:31:12 INTERCEPTED kill of head job 'nf-workflow-<sessionID>' (fake job <PID>)
2026-08-19 12:31:14 -> PID <PID> stopped
```

and confirm the process is actually gone with `ps -fp <PID>`. Cancelling a
real pipeline task instead logs a `PASSTHROUGH(bkill)` line and behaves
exactly like plain `bkill`.

## Known limitation

Both wrappers rely on `tw-agent` resolving `bsub`/`bkill` via a plain
`PATH` lookup. If Seqera ever changes `tw-agent` to call an absolute path,
or bundles its own resolution logic for either binary, this stops applying
and the relevant wrapper is simply never hit (harmless - `bsub`/`bkill` for
everything else keeps working normally either way, it just stops
intercepting head jobs). There's no way to verify this ahead of time other
than triggering a real run and checking `wrapper.log`, as above.

`bkill-wrapper.sh` matches fake job IDs by treating any bare numeric
argument as a candidate PID, and any other argument as a candidate job
name - it doesn't parse `bkill`'s full flag grammar. This mirrors how
`bsub-wrapper.sh` only looks for the `#BSUB -J nf-workflow-` tag rather
than fully parsing LSF directives, and is deliberately best-effort rather
than a complete `bkill` reimplementation.
