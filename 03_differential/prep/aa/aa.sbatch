#!/bin/bash
#SBATCH --job-name=aa_colo320
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-2
#SBATCH --output=logs/aa_%a.log

# AmpliconSuite-pipeline (CNVkit seeds -> AmpliconArchitect -> AmpliconClassifier) on the Wu 2019 full-depth WGS.
# AA downsamples to 10x by default, which is what it is tuned for.
set -euo pipefail
module load miniforge/25.11.0-0
conda activate ampliconsuite_ac2       # AmpliconSuite-pipeline (bioconda "ampliconsuite", with a free MOSEK academic license)
export AA_DATA_REPO=$SLURM_SUBMIT_DIR/aa_data_repo   # GRCh38 only: the shared repo's GRCh37 folder is incomplete
export MOSEKLM_LICENSE_FILE=$HOME/mosek/mosek.lic
cd "$SLURM_SUBMIT_DIR"

names=(COLO320DM COLO320HSR); s=${names[$((SLURM_ARRAY_TASK_ID - 1))]}
AmpliconSuite-pipeline.py -s "$s" -t 16 --bam ../wgs/${s}.bam --ref GRCh38 \
    --run_AA --run_AC -o "$s"
