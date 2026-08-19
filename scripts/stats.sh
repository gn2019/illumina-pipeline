#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 100G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --job-name stats
#SBATCH --time 02:00:00

set -x
set -e

# conda activate samtools
# samtools 1.16.1
# mosdepth 0.3.3
# picard 2.27.5


sample=$1
bam=$2
outdir=stats

mkdir -p $outdir

mosdepth -t 10 --no-per-base --fast-mode -Q 5 $outdir/stats $bam

cat $outdir/stats.mosdepth.summary.txt | head -24 | tail -n+2 | awk '{sum+=$4} END {print "mean_coverage =",sum/NR}' >> $outdir/stats.mosdepth.summary.txt  

samtools stats -d -F 2304 -@10 $bam > $outdir/samtools.stats


