# Installing the LSF head-job preemption guard

Run everything below as the WEXAC user that runs `tw-agent` (i.e. the one
Seqera Platform's compute environment is configured to use).

## 1. Run the installer

```bash
cd tools/lsf-headjob-guard
./install.sh
```

This writes `~/bin/bsub` (the wrapper, with the real `bsub` path baked in).
Safe to re-run any time.

## 2. Restart tw-agent

`tw-agent` only picks up the wrapper if it's *launched* with `~/bin` ahead
of the real `bsub` in `PATH` - a running instance keeps whatever
environment it already had, so it needs a restart with that PATH set
explicitly.

Find the running instance and note its exact arguments:

```bash
pgrep -fa tw-agent
```

You'll see something like:

```
123456 ./tw-agent WEXAC --access-token=eyJ0aWQiOiAx...
```

Kill it, then restart it from the same directory with `~/bin` forced to
the front of `PATH`, using the exact same `<ENV_NAME>` and
`--access-token=...` from the `pgrep -fa` output (don't regenerate a
token):

```bash
kill 123456
cd /path/to/your/tw-agent/dir
PATH="$HOME/bin:$PATH" ./tw-agent WEXAC --access-token=eyJ0aWQiOiAx...
```

If `tw-agent` normally runs inside `tmux`/`screen` for persistence, run
this restart inside that same session so it survives you disconnecting,
exactly as before.

## 3. Verify

```bash
which bsub                          # -> ~/bin/bsub
tail -f ~/.bsub_wrapper/wrapper.log
```

Trigger a launch or Resume in Seqera Platform and watch for:

```
INTERCEPTED head job 'nf-workflow-<sessionID>'
```

If instead you only ever see `PASSTHROUGH(...)` lines even after a
launch/resume, `tw-agent` isn't resolving `bsub` through the wrapper - see
the "Known limitation" section in `README.md`.

## Installing for a different user

Everything above is per-user (each user's own `~/bin`, own `tw-agent`
instance, own `--access-token`). To set this up for someone else, they
need to run steps 1-3 themselves under their own account - there's nothing
here that needs root or any shared/system-level change.
