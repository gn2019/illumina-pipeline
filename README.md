# Illumina Pipeline

Builds bam files from fastq files, then runs genome analysis tools.

## One-time setup

Run everything at once:

```
bash install_all.sh
```

...or the individual pieces, if you only need some of them:

```
bash install_paths.sh             # creates results/refs/data/scripts dirs under
                                   # $HOME/illumina-pipeline, copies scripts/ into place
bash environments/setup_envs.sh   # creates lumpy-sv, bw, ampsuite under $HOME/.conda/envs
bash containers/build.sh          # builds cgpwgs.sif into the repo root
```

`environments/*.yml` are the source definitions; `setup_envs.sh` and
`build.sh` just land the built envs/image at the paths `nextflow.config`
already points at, so no config changes are needed when re-running them.

`install_paths.sh` only creates directories and copies `scripts/` - the
genome fasta, fastq files, and ASCAT GC-correction reference still need to
be placed by hand (it prints exactly where when it runs).

## Running via Seqera Platform on WEXAC

If launches/resumes through Seqera Platform die with a `403 Unable to
access config file` error shortly after LSF shows a preemption, see
`tools/lsf-headjob-guard/README.md` - it's a one-time per-user setup that
stops LSF from ever being able to preempt the pipeline's own Nextflow head
process in the first place.
