#!/bin/bash
#SBATCH --job-name=s2p_prepare
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/prepare_%j.log
set -euo pipefail
module load miniforge/25.11.0-0
conda activate scprinter                 # made by seq2print/build_env.sh
cd "$SLURM_SUBMIT_DIR"
python prepare.py
