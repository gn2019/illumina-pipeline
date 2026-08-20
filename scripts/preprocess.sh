#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 100G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --job-name pre
#SBATCH --output slurm-preprocess.txt

set -x
set -e

# conda activate ofer
#- python 3.7
#- trim-galore 0.6.10 bioconda/noarch
#- fastqc      0.12.1 bioconda/noarch
#- bwa         0.7.17 bioconda/linux-64
#- picard      3.1.1 bioconda/noarch
#- gatk        4.4.0.0 bioconda/noarch

date
hostname

sample=$1
results=$2
fastq1=$3
fastq2=$4
genome=$5

name1=`echo "$fastq1" | awk -F"[/\.]" '{print $(NF-2)}'`
name2=`echo "$fastq2" | awk -F"[/\.]" '{print $(NF-2)}'`

# Per-lane run id (e.g. "ERR3085889_1" -> "ERR3085889"). Derived from the
# fastq filename itself rather than a caller-supplied argument, so it's
# always correct and always unique per read-pair/lane - a sample with
# multiple lanes (like SAMEA104693125 -> ERR3085889 + ERR3085891) gets one
# distinct run id per lane instead of everything colliding on a single
# hardcoded name.
run="${name1%_1}"

# Nextflow gives every task attempt (including retries and non-cached
# -resume launches) a brand-new, EMPTY work directory - that's fundamental
# to how it isolates and caches tasks. A checkpoint file living only inside
# that ephemeral directory can never survive a retry, so it never actually
# skips anything. To make checkpointing meaningful, all real work below
# happens in a persistent, per-sample-per-lane directory on shared storage
# (under params.results, passed in as $results) that outlives any single
# task attempt.
#
# NOTE: this assumes a given (sample, run) pair is never processed by two
# concurrently-running task attempts at once (e.g. two overlapping -resume
# launches on the same run) - concurrent writers to the same persist dir
# would race. Normal retry-after-failure and -resume-after-relaunch usage
# is safe since only one attempt is ever actively running at a time.
persist="${results}/${sample}/preprocess/${run}"
mkdir -p "${persist}/qc_raw" "${persist}/qc_trimmed" "${persist}/trimmed_reports"

CKPT="${persist}/.checkpoints"
touch "$CKPT"

has_ckpt () { grep -qx "$1" "$CKPT"; }
add_ckpt () { echo "$1" >> "$CKPT"; }

##get quality of fastqs
if ! has_ckpt fastqc_raw; then
  date
  fastqc -o "${persist}/qc_raw" ${fastq1} ${fastq2}
  add_ckpt fastqc_raw
fi

##trim reads
if ! has_ckpt trim; then
  date
  trim_galore -q 20 --phred33 --illumina --paired -j 8 --length 20 --output_dir "${persist}/trimmed_reports" ${fastq1} ${fastq2}
  rm -rf "${persist}/trimmed_reports"/*_trimmed.fq.gz
  mv "${persist}/trimmed_reports"/*.fq.gz "${persist}/"
  mv "${persist}/${name1}_val_1.fq.gz" "${persist}/${run}_1.fq.gz"
  mv "${persist}/${name2}_val_2.fq.gz" "${persist}/${run}_2.fq.gz"
  add_ckpt trim
fi

##mapping
if ! has_ckpt map; then
  date
  bwa mem -q -t 8 -M \
  -R $(echo "@RG\tID:${run}\tSM:${sample}\tLB:lib_${sample}\tPL:ILLUMINA") \
  ${genome} \
  "${persist}/${run}_1.fq.gz" \
  "${persist}/${run}_2.fq.gz" \
  | samtools sort -@8 -O BAM -o "${persist}/${run}_sorted.bam" -

  ##index bam file
  samtools index "${persist}/${run}_sorted.bam"

  add_ckpt map
fi

##mark duplicates
if ! has_ckpt dedup; then
  date
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
	I="${persist}/${run}_sorted.bam" \
	O="${persist}/${run}_sorted.dedup.bam" \
	M="${persist}/${run}_sorted.dedup.metrix"

  ##index bam file
  samtools index "${persist}/${run}_sorted.dedup.bam"

  add_ckpt dedup
fi

##get quality of trimmed reads
if ! has_ckpt fastqc_trimmed; then
  date
  fastqc -o "${persist}/qc_trimmed" "${persist}/${run}_1.fq.gz" "${persist}/${run}_2.fq.gz"
  add_ckpt fastqc_trimmed
fi

##create coverage overview in bigwig format
if ! has_ckpt coverage; then
  date
  bamCoverage --bam "${persist}/${run}_sorted.dedup.bam" -o "${persist}/${run}_sorted.dedup.bw"
  add_ckpt coverage
fi

## remove large pre-dedup sorted bam once dedup has succeeded
rm -f "${persist}/${run}_sorted.bam"*

# --- Hand the result back to THIS Nextflow task's own (ephemeral) work dir.
# main.nf's PREPROCESS process declares its output as
# path("${meta.id}/*.dedup.bam") - i.e. it expects a directory named after
# the sample, not "noERX" or anything else. Symlinking (cheap, no data copy)
# the persisted dedup bam in here is what makes that output glob
# actually match, and what publishDir then copies out.
mkdir -p "${sample}"
ln -f "${persist}/${run}_sorted.dedup.bam" "${sample}/${run}_sorted.dedup.bam"
[ -f "${persist}/${run}_sorted.dedup.bam.bai" ] && ln -f "${persist}/${run}_sorted.dedup.bam.bai" "${sample}/${run}_sorted.dedup.bam.bai"

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
