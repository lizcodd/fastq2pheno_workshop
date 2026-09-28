#!/bin/bash
#SBATCH --job-name=cbp_modisco
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/modisco_%a.log
# TF-MoDISco on each model's contribution scores: task 1 = DM counts, 2 = DM profile, 3 = HSR counts, 4 = HSR profile.
# CPU only, 3-6 h each. Output: ../results/chrombpnet/contribs/<line>_fold0/modisco_<counts|profile>{.h5,_report/}
# Submit from this folder after model.sbatch:  sbatch modisco.sbatch
set -euo pipefail
module load apptainer/1.5.2
cd "$SLURM_SUBMIT_DIR"
SIF=${SIF:-chrombpnet_latest.sif}
cb() { apptainer exec "$SIF" "$@"; }

line=$(echo COLO320DM COLO320DM COLO320HSR COLO320HSR | cut -d' ' -f$SLURM_ARRAY_TASK_ID)
k=$([ $((SLURM_ARRAY_TASK_ID % 2)) -eq 1 ] && echo counts || echo profile)
id=${line}_fold0; d=../results/chrombpnet/contribs/$id
MEME=$(cb python -c "from chrombpnet.data import DefaultDataFile, get_default_data_path; print(get_default_data_path(DefaultDataFile.motifs_meme))")
[ -f $d/modisco_$k.h5 ] || cb modisco motifs -i $d/$id.${k}_scores.h5 -n 500000 -w 500 -o $d/modisco_$k.h5
cb modisco report -i $d/modisco_$k.h5 -o $d/modisco_${k}_report/ -m "$MEME"
echo "$id $k done"
