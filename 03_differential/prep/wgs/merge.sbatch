#!/bin/bash
#SBATCH --job-name=wgs_merge
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=1-2
#SBATCH --output=logs/merge_%a.log
# Merge one line's chunks, remove duplicates, index. Chunks are left in place for checking.
set -euo pipefail
module load miniforge/25.11.0-0; conda activate fastq2pheno
cd "$SLURM_SUBMIT_DIR"
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" runs.tsv | cut -f2)
n_expected=$(awk -v n="$name" '$2 == n' chunks.tsv | wc -l)
n_found=$(ls chunks/${name}_*.bam | grep -vc tmp)
[ "$n_found" -eq "$n_expected" ] || { echo "only $n_found of $n_expected chunks"; exit 1; }
samtools merge -@ 16 -o - chunks/${name}_*.bam | samtools markdup -r -@ 8 -s - ${name}.bam
samtools index -@ 8 ${name}.bam
samtools flagstat -@ 8 ${name}.bam > ${name}.flagstat.txt
cat ${name}.flagstat.txt
