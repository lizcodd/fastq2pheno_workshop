#!/bin/bash
#SBATCH --job-name=s2p_tfbs
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-2
#SBATCH --output=logs/tfbs_%a.log
# seq2PRINT TF binding scores for each line's model (about 4 GPU-hours each). Submit after train.sbatch:  sbatch tfbs.sbatch
set -euo pipefail
module load miniforge/25.11.0-0
conda activate scprinter                 # made by seq2print/build_env.sh
cd "$SLURM_SUBMIT_DIR"
python tfbs.py $(echo COLO320DM COLO320HSR | cut -d' ' -f$SLURM_ARRAY_TASK_ID)
ls -la ../results/seq2print/*TFBS*
