#!/bin/bash
#SBATCH --job-name=cbp_prep
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/prep_%j.log
# Build 04_sequence_models/data/ and chrombpnet_data.tar.gz from ../results/ (the trained chromBPNet models) and Step 3's data.
# Submit from this folder:  sbatch run_prep.sbatch      (CPU only; about 3 hours)
set -euo pipefail
module load miniforge/25.11.0-0 apptainer/1.5.2
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs work/contribs
SIF=${SIF:-../chrombpnet/chrombpnet_latest.sif}
# the container sees the workshop folder, plus the genome's real location if it is a link to elsewhere
G=../../02_pipeline/results/genome/genome.fa
cb() { apptainer exec -B "$(readlink -f ../..)" -B "$(dirname "$(readlink -f $G)")" "$SIF" "$@"; }

# Step 3's DESeq2 calls, to pick featured peaks and for the DM vs HSR comparison
conda activate fastq2pheno
Rscript calls.R

# featured peaks, motif-planting backgrounds
cb python prepare_bundle.py select

# contribution scores of each line's model on the featured peaks
grep -P '^chr([0-9]+|X)\t' $G.fai | cut -f1,2 > work/chrom.sizes
while read -r line model; do
    out=work/contribs/chrombpnet_${line}
    [ -f $out.counts_scores.bw ] || cb chrombpnet contribs_bw -m "$model" -r work/featured.narrowPeak -g $G \
        -c work/chrom.sizes -pc counts -op $out
done < work/chrombpnet_models.txt

# every peak, and every motif site in peaks with the site (or a random 6 bp) scrambled; both models predict them all
[ -f work/peak_seqs.txt ] || cb python peaks.py
[ -f work/ablate_seqs.txt ] || cb python ablation.py
ls work/pred_chrombpnet_COLO320{DM,HSR}_{peak,ablate}.npy >/dev/null 2>&1 || cb python predict.py

cb python prepare_bundle.py build
