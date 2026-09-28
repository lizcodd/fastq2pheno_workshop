#!/bin/bash
#SBATCH --job-name=wgs_align
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=logs/align_%a.log
# One chunk of 25M read pairs per array task; chunks.tsv lists run, name, first spot, last spot.
set -euo pipefail
module load miniforge/25.11.0-0; conda activate fastq2pheno
cd "$SLURM_SUBMIT_DIR"
read -r run name first last < <(sed -n "${SLURM_ARRAY_TASK_ID}p" chunks.tsv)
ref=$(grep '^fasta:' ../../../02_pipeline/params.yaml | awk '{print $2}')
out=chunks/${name}_$(printf %03d $SLURM_ARRAY_TASK_ID)
mkdir -p chunks
fastq-dump --split-spot --stdout -N "$first" -X "$last" sra/"$run" \
  | bwa mem -p -t 16 -R "@RG\tID:${name}\tSM:${name}" "$ref" - \
  | samtools fixmate -m -@ 4 - - \
  | samtools sort -@ 8 -m 1500M -o ${out}.tmp.bam -
mv ${out}.tmp.bam ${out}.bam
echo "chunk $SLURM_ARRAY_TASK_ID ($name spots $first-$last) done"
