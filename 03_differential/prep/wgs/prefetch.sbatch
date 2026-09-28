#!/bin/bash
#SBATCH --job-name=wgs_prefetch
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --array=1-2
#SBATCH --output=logs/prefetch_%a.log
set -euo pipefail
module load miniforge/25.11.0-0; conda activate fastq2pheno
cd "$SLURM_SUBMIT_DIR"
run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" runs.tsv | cut -f1)
prefetch --max-size u --output-directory sra "$run"
vdb-validate sra/"$run"
