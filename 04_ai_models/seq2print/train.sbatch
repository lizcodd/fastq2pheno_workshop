#!/bin/bash
#SBATCH --job-name=s2p_train
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00         # seq2PRINT allows 300 epochs (early stop after 5 without improvement), then attribution
#SBATCH --array=1-2
#SBATCH --output=logs/train_%a.log
# One seq2PRINT model per line (both replicates pooled, fold 0). wandb off; --run_id fixes the model file name:
#   ../results/seq2print/model/COLO320_<line>_fold0-<line>_fold0.pt  (savename comes from the config)
# Submit from this folder after configs.sbatch:  sbatch train.sbatch
set -euo pipefail
module load miniforge/25.11.0-0
conda activate scprinter                 # made by seq2print/build_env.sh
cd "$SLURM_SUBMIT_DIR"
W=$(readlink -f ../results/seq2print)
line=$(echo COLO320DM COLO320HSR | cut -d' ' -f$SLURM_ARRAY_TASK_ID)
mkdir -p $W/temp $W/model
seq2print_train --config $W/configs/COLO320_${line}_fold0.JSON --temp_dir $W/temp --model_dir $W/model \
    --data_dir $W --run_id ${line}_fold0
echo "$line done"
