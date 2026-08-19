#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 100G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --job-name pre
#SBATCH --output slurm-cnv.txt

set -x
set -e

# conda activate lumpy-sv
#- python 2.7
#- lumpy 0.2.13
#- samtools 1.10
#- survivor 1.0.7
#- vcftools 0.1.16

date
hostname

sample=$1
bam=$2
genome=$3
exclude=$4

threads=16

mkdir -p lumpy

CKPT=.sv_checkpoints
touch "$CKPT"

has_ckpt () { grep -qx "$1" "$CKPT"; }
add_ckpt () { echo "$1" >> "$CKPT"; }


############### Preprocessing

if ! has_ckpt preprocess; then
  date

  # Resolve extractSplitReads_BwaMem
  if command -v extractSplitReads_BwaMem >/dev/null 2>&1; then
    EXTRACT_SPLIT_READS=$(command -v extractSplitReads_BwaMem)
  else
    EXTRACT_SPLIT_READS=$(ls -d "$CONDA_PREFIX"/share/lumpy-sv-*/scripts/extractSplitReads_BwaMem 2>/dev/null | sort -V | tail -1)
  fi
  [[ -z "$EXTRACT_SPLIT_READS" ]] && { echo "ERROR: extractSplitReads_BwaMem not found in PATH or under \$CONDA_PREFIX/share/lumpy-sv-*/scripts" >&2; exit 1; }

  # Extract the discordant paired-end alignments.
  samtools view -b -F 1294 $bam > lumpy/discordants.unsorted.bam

  # Extract the split-read alignments
  samtools view -h $bam | "$EXTRACT_SPLIT_READS" -i stdin | samtools view -Sb - > lumpy/splitters.unsorted.bam

  # Sort both alignments
  samtools sort -@ $threads -o lumpy/discordants.bam lumpy/discordants.unsorted.bam
  samtools sort -@ $threads -o lumpy/splitters.bam lumpy/splitters.unsorted.bam

  # Index
  samtools index -@ $threads lumpy/discordants.bam
  samtools index -@ $threads lumpy/splitters.bam
  rm lumpy/discordants.unsorted.bam lumpy/splitters.unsorted.bam

  add_ckpt preprocess
fi

############### Run Lumpy

if ! has_ckpt lumpy; then
  date

  lumpyexpress \
    -B $bam \
    -S lumpy/splitters.bam \
    -D lumpy/discordants.bam \
    -x $exclude \
    -R $genome \
    -o $bam.vcf

  mv $bam.vcf lumpy/lumpy_sv.vcf

  add_ckpt lumpy
fi

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
