#!/bin/bash
#SBATCH --job-name=download_fastq
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=download_%j.log

# Download the COLO320 ATAC-seq runs from NCBI SRA and convert them to gzipped FASTQ.
# Submit from this folder:  sbatch download.sbatch
set -euo pipefail

module load miniforge/25.11.0-0
conda activate fastq2pheno       # from environment.yml

cd "$SLURM_SUBMIT_DIR/fastq"

# Run accessions and expected read pairs, from samples.tsv (columns 1 and 5)
while read -r run expected; do
    echo "== $run"
    prefetch --max-size u --output-directory sra "$run" < /dev/null     # download sra/<run>/<run>.sra
    fasterq-dump --split-files --threads 8 --temp tmp sra/"$run" < /dev/null # write <run>_1.fastq, <run>_2.fastq
    pigz -p 8 "${run}_1.fastq" "${run}_2.fastq"                          # compress to .fastq.gz

    # No archive checksums for FASTQs made this way, so check that every read arrived instead
    n=$(( $(pigz -dc "${run}_1.fastq.gz" | wc -l) / 4 ))
    if [ "$n" -eq "$expected" ]; then
        echo "$run: $n read pairs, OK"
    else
        echo "$run: $n read pairs, expected $expected: FAILED"
        exit 1
    fi
done < <(tail -n +2 ../samples.tsv | cut -f1,5)

echo "All runs downloaded. The sra/ folder can now be deleted."
