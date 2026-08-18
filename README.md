# Illumina Pipeline

Builds bam files from fastq files, then runs genome analysis tools.

## One-time environment setup

The pipeline references its conda envs and the cgpwgs container by fixed path
(see the "Environments & Containers" section of `nextflow.config`) rather than
building them as part of the DAG. Set them up once before running:

```
bash environments/setup_envs.sh   # creates lumpy-sv, bw, ampsuite under $HOME/.conda/envs
bash containers/build.sh          # builds cgpwgs.sif into the repo root
```

`environments/*.yml` are the source definitions; `setup_envs.sh` and
`build.sh` just land the built envs/image at the paths `nextflow.config`
already points at, so no config changes are needed when re-running them.
