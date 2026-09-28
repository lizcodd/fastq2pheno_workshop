#!/bin/bash
#SBATCH --job-name=atacseq
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=atacseq_%j.log

# Run nf-core/atacseq on the COLO320 samples. Submit from this folder:  sbatch run.sbatch
# This job only runs Nextflow; Nextflow submits every pipeline step as its own SLURM job.
set -euo pipefail

module load miniforge/25.11.0-0 apptainer/1.5.2
conda activate fastq2pheno       # from environment.yml

export NXF_SINGULARITY_CACHEDIR=$PWD/containers   # container images, kept outside work/
export APPTAINER_CACHEDIR=$PWD/containers/tmp      # apptainer's download cache, out of the small home directory
export NXF_SYNTAX_PARSER=v1                        # nf-core/atacseq 2.1.2 needs Nextflow's older parser

# local.config (not in git; see local.config.example): the reference genome paths, and optionally a lab QOS
nextflow run nf-core/atacseq -r 2.1.2 -profile engaging -params-file params.yaml -c local.config -resume
