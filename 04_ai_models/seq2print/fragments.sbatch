#!/bin/bash
#SBATCH --job-name=s2p_frags
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/fragments_%a.log
# Fragment file per replicate (10x style: chrom, start, end, sample; Tn5-shifted +4/-5) from the pipeline's filtered BAMs,
# in ../results/seq2print/fragments/. Submit from this folder:  sbatch fragments.sbatch
set -euo pipefail
module load miniforge/25.11.0-0
conda activate fastq2pheno
cd "$SLURM_SUBMIT_DIR"
F=../results/seq2print/fragments; mkdir -p $F logs
s=$(echo COLO320DM_REP1 COLO320DM_REP2 COLO320HSR_REP1 COLO320HSR_REP2 | cut -d' ' -f$SLURM_ARRAY_TASK_ID)
bam=../../02_pipeline/results/bwa/merged_library/$s.mLb.clN.sorted.bam
tmp=${TMPDIR:-/tmp}/$s.$SLURM_JOB_ID; mkdir -p $tmp
samtools sort -n -@ 8 -m 2G -T $tmp/ns -o $tmp/byname.bam $bam
bedtools bamtobed -bedpe -i $tmp/byname.bam 2> $tmp/bedpe.warn \
  | awk -v OFS='\t' -v s=$s '$1 == $4 && $1 ~ /^chr([0-9]+|X)$/ {
        st = ($2 < $5 ? $2 : $5) + 4; en = ($3 > $6 ? $3 : $6) - 5
        if (en > st) print $1, st, en, s }' \
  | sort -k1,1 -k2,2n -S 8G -T $tmp --parallel=8 | bgzip -@ 8 > $F/$s.frags.tsv.gz
tabix -p bed $F/$s.frags.tsv.gz
echo "$s: $(zcat $F/$s.frags.tsv.gz | wc -l) fragments"
rm -r $tmp
