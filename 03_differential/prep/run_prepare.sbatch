#!/bin/bash
#SBATCH --job-name=prep_differential
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/prepare_%j.log
# Run prepare_data.ipynb (writes ../data/). Submit from this folder: sbatch run_prepare.sbatch
set -euo pipefail
module load miniforge/25.11.0-0
conda activate fastq2pheno
export R_LIBS_USER=/dev/null
cd "$SLURM_SUBMIT_DIR"
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 prepare_data.ipynb
