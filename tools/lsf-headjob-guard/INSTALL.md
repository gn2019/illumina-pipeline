# Installing the LSF head-job preemption guard

Run everything below as the WEXAC user that runs `tw-agent` (i.e. the one
Seqera Platform's compute environment is configured to use).

## 1. Run the installer

```bash
cd tools/lsf-headjob-guard
./install.sh
```

This writes `~/bin/bsub` and `~/bin/bkill` (the wrappers, with the real
`bsub`/`bkill` paths baked in) and creates
`~/.bsub_wrapper/intercepted_jobs`, the registry the two wrappers share.
Safe to re-run any time.

## 2. Restart tw-agent

`tw-agent` only picks up the wrappers if it's *launched* with `~/bin`
ahead of the real `bsub`/`bkill` in `PATH` - a running instance keeps
whatever environment it already had, so it needs a restart with that PATH
set explicitly.

Find the running instance and note its exact arguments:

```bash
pgrep -fa tw-agent
```

You'll see something like:

```
123456 ./tw-agent WEXAC --access-token=eyJ0aWQiOiAx...
```

Kill it, then restart it with `run_tw_agent.sh` instead of `./tw-agent`
directly, using the exact same `<ENV_NAME>` and `--access-token=...` from
the `pgrep -fa` output (don't regenerate a token):

```bash
kill 123456
cp tools/lsf-headjob-guard/run_tw_agent.sh /path/to/your/tw-agent/dir/
cd /path/to/your/tw-agent/dir
./run_tw_agent.sh WEXAC --access-token=eyJ0aWQiOiAx...
```

`run_tw_agent.sh` sets `PATH` correctly and also restarts `tw-agent`
automatically if it exits - `tw-agent` has a known upstream bug
([tower-agent#57](https://github.com/seqeralabs/tower-agent/issues/57),
still open as of v0.5.5) where it sometimes drops its connection and then
fails to reconnect with `There is an active agent for this user and
connection ID`, needing a restart to recover. It waits 10s between
restarts so it doesn't tight-loop if `tw-agent` ever exits immediately for
some other reason; `Ctrl-C` stops it for good.

If `tw-agent` normally runs inside `tmux`/`screen` for persistence, run
`run_tw_agent.sh` inside that same session so it survives you
disconnecting, exactly as before.

## 3. Verify

```bash
which bsub                          # -> ~/bin/bsub
which bkill                         # -> ~/bin/bkill
tail -f ~/.bsub_wrapper/wrapper.log
```

Trigger a launch or Resume in Seqera Platform and watch for:

```
INTERCEPTED head job 'nf-workflow-<sessionID>'
```

Then hit Cancel on that run in Seqera Platform and watch for:

```
INTERCEPTED kill of head job 'nf-workflow-<sessionID>' (fake job <PID>)
-> PID <PID> stopped
```

If instead you only ever see `PASSTHROUGH(...)` lines even after a
launch/resume/cancel, `tw-agent` isn't resolving `bsub`/`bkill` through the
wrappers - see the "Known limitation" section in `README.md`.

## Installing for a different user

Everything above is per-user (each user's own `~/bin`, own `tw-agent`
instance, own `--access-token`). To set this up for someone else, they
need to run steps 1-3 themselves under their own account - there's nothing
here that needs root or any shared/system-level change.
