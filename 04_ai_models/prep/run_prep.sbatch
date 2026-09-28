#!/bin/bash
#SBATCH --job-name=ai_prep
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/prep_%j.log
# Build 04_ai_models/data/ and ai_models_data.tar.gz from ../results/ (the trained models) and Step 3's data.
# Submit from this folder:  sbatch run_prep.sbatch
set -euo pipefail
module load miniforge/25.11.0-0 apptainer/1.5.2
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs work/contribs

# DESeq2 calls (as in Step 3), to pick the featured peaks
conda activate fastq2pheno
Rscript calls.R

conda activate scprinter                 # made by seq2print/build_env.sh
python prepare_bundle.py select
# chromBPNet contribution scores of each line's model on the featured peaks (CPU is enough for ~750 peaks)
SIF=${SIF:-../chrombpnet/chrombpnet_latest.sif}
G=../../02_pipeline/results/genome/genome.fa
grep -P '^chr([0-9]+|X)\t' $G.fai | cut -f1,2 > work/chrom.sizes
while read -r line k model; do
    out=work/contribs/chrombpnet_${line}_${k}
    [ -f $out.counts_scores.bw ] || apptainer exec "$SIF" chrombpnet contribs_bw -m "$model" \
        -r work/featured.narrowPeak -g $G -c work/chrom.sizes -pc counts -op $out
done < work/chrombpnet_models.txt
python prepare_bundle.py build
