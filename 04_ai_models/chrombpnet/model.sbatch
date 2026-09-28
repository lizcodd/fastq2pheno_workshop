#!/bin/bash
#SBATCH --job-name=cbp_model
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-2
#SBATCH --output=logs/model_%a.log
# One chromBPNet model per line (both replicates pooled, fold 0): train, evaluate, then contribution scores on all
# its peaks (counts + profile). About 12 GPU-hours per model on an L40S. The steps are separate, and each is skipped
# when its output exists, so a job stopped by a GPU time limit continues where it left off when resubmitted.
# Without time limits, `chrombpnet pipeline` does training, evaluation and profile motifs in one command.
# Output: ../results/chrombpnet/models/<line>_fold0/ and ../results/chrombpnet/contribs/<line>_fold0/
# Submit from this folder after prepare.sbatch:  sbatch model.sbatch
set -euo pipefail
module load apptainer/1.5.2
cd "$SLURM_SUBMIT_DIR"
SIF=${SIF:-chrombpnet_latest.sif}
cb() { apptainer exec --nv "$SIF" "$@"; }

line=$(echo COLO320DM COLO320HSR | cut -d' ' -f$SLURM_ARRAY_TASK_ID)
id=${line}_fold0
I=../results/chrombpnet/inputs; m=../results/chrombpnet/models/$id; c=../results/chrombpnet/contribs/$id

if [ ! -f $m/models/chrombpnet_nobias.h5 ]; then
    [ -d "$m" ] && { echo "$m exists but is incomplete; move it aside"; exit 1; }
    cb chrombpnet train -ibam $I/$line.main.bam -d ATAC -g $I/genome.fa -c $I/chrom.sizes \
        -p $I/${line}_peaks.narrowPeak -n $I/${line}_nonpeaks_negatives.bed \
        -fl $I/folds/fold_0.json -b $I/bias_fold_0.h5 -o $m --tmpdir "${TMPDIR:-/tmp}"
    echo "$id trained"
fi
# Evaluation on the test chromosomes (qc_direct.py runs chromBPNet's QC step on the training folder)
if [ ! -f $m/evaluation/chrombpnet_metrics.json ]; then
    cb python qc_direct.py -bw $m/auxiliary/data_unstranded.bw -cm $m/models/chrombpnet.h5 -cmb $m/models/chrombpnet_nobias.h5 \
        -g $I/genome.fa -c $I/chrom.sizes -d ATAC -p $m/auxiliary/filtered.peaks.bed -n $m/auxiliary/filtered.nonpeaks.bed \
        -fl $I/folds/fold_0.json -o $m
    echo "$id evaluated"
fi
# Contribution scores of the bias-free model on all of the line's peaks
if [ ! -f $c/$id.profile_scores.bw ]; then
    mkdir -p $c
    cb chrombpnet contribs_bw -m $m/models/chrombpnet_nobias.h5 -r $m/auxiliary/filtered.peaks.bed \
        -g $I/genome.fa -c $I/chrom.sizes -pc counts profile -op $c/$id
    echo "$id contribution scores done"
fi
