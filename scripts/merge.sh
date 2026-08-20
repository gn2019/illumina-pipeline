#!/bin/bash
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32000]"

set -x
set -e

sample=$1          # e.g. PD29429h_DM3_40
bam_list=$2        # text file with paths to *.dedup.bam

merged_bam=${sample}.bam

echo "Merging BAMs for ${sample}"
samtools merge -@ "${LSB_DJOB_NUMPROC:-1}" -b "${bam_list}" "${merged_bam}"

echo "Indexing merged BAM"
samtools index "${merged_bam}"

echo "Done: ${merged_bam}"
