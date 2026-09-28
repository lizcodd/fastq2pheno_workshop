#!/bin/bash
#SBATCH --job-name=cbp_prepare
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/prepare_%j.log
# Prepare chromBPNet inputs for each line from the nf-core/atacseq results (Step 2), in ../results/chrombpnet/inputs.
# Submit from this folder:  sbatch prepare.sbatch
set -euo pipefail
module load miniforge/25.11.0-0 apptainer/1.5.2
conda activate fastq2pheno
cd "$SLURM_SUBMIT_DIR"

RES=../../02_pipeline/results
I=../results/chrombpnet/inputs
SIF=${SIF:-chrombpnet_latest.sif}                 # apptainer pull chrombpnet_latest.sif docker://kundajelab/chrombpnet:latest
cb() { apptainer exec "$SIF" "$@"; }
mkdir -p $I logs

# Train/validation/test chromosome split (chromBPNet's Zenodo record 7443683); fold 0 is used
[ -f $I/folds/fold_0.json ] || { curl -sL -o $I/folds.zip "https://zenodo.org/records/7443683/files/folds.zip?download=1"
                                  unzip -o -q $I/folds.zip -d $I; }
# Tn5 bias model: ENCODE's K562 ATAC bias model, fold 0 (annotation ENCSR467RSV, file ENCFF984RAF, 727 MB)
if [ ! -f $I/bias_fold_0.h5 ]; then
    curl -sL -o $I/ENCFF984RAF.tar.gz \
        https://encode-public.s3.amazonaws.com/2024/09/06/356dd29e-20d8-4647-a484-2b62c7689fe3/ENCFF984RAF.tar.gz
    echo "0da121e548f1fc977e6ea206ef5b7a52  $I/ENCFF984RAF.tar.gz" | md5sum -c -
    tar -xzf $I/ENCFF984RAF.tar.gz -C $I ./fold_0/model.bias_scaled.fold_0.ENCSR868FGK.h5
    mv $I/fold_0/model.bias_scaled.fold_0.ENCSR868FGK.h5 $I/bias_fold_0.h5 && rmdir $I/fold_0 && rm $I/ENCFF984RAF.tar.gz
fi

# Genome: a link to the pipeline's FASTA plus a fresh index next to it (chromBPNet rebuilds an index it thinks is stale)
ln -sfn "$(readlink -f $RES/genome/genome.fa)" $I/genome.fa
samtools faidx $I/genome.fa
grep -P '^chr([0-9]+|X)\t' $I/genome.fa.fai | cut -f1,2 > $I/chrom.sizes      # chr1-22 and X only

# Regions to leave out: the ENCODE blacklist, and the amplicons (amplicons.bed: AmpliconArchitect intervals of both
# lines, padded; chromBPNet predicts reads from sequence and can't know a region has 100 copies)
cat <(grep -P '^chr([0-9]+|X)\t' ../../02_pipeline/hg38-blacklist.v3.bed | cut -f1-3) <(cut -f1-3 amplicons.bed) \
    | sort -k1,1 -k2,2n > $I/exclude.bed
bedtools slop -i $I/exclude.bed -g $I/chrom.sizes -b 1057 > $I/exclude_ext.bed

for s in COLO320DM COLO320HSR; do
    # Reads (both replicates merged) on chr1-22 and X only: chromBPNet can't handle the small unplaced contigs
    samtools view -@ 8 -b -o $I/$s.main.bam $RES/bwa/merged_replicate/$s.mRp.clN.sorted.bam $(cut -f1 $I/chrom.sizes)
    samtools index $I/$s.main.bam
    # Peaks: the line's MACS2 peaks whose 2,114 bp model window stays clear of excluded regions and chromosome ends
    grep -P '^chr([0-9]+|X)\t' $RES/bwa/merged_replicate/macs2/narrow_peak/$s.mRp.clN_peaks.narrowPeak \
      | bedtools intersect -v -a - -b $I/exclude_ext.bed \
      | awk -v OFS='\t' -v sizes=$I/chrom.sizes 'BEGIN { while ((getline < sizes) > 0) L[$1] = $2 }
                                                  { c = $2 + $10; if (c - 1057 > 0 && c + 1057 < L[$1]) print }' \
      | sort -k1,1 -k2,2n > $I/${s}_peaks.narrowPeak
    # Background regions matched to the peaks' GC content (fold 0)
    cb chrombpnet prep nonpeaks -g $I/genome.fa -p $I/${s}_peaks.narrowPeak -c $I/chrom.sizes \
        -fl $I/folds/fold_0.json -br $I/exclude.bed -o $I/${s}_nonpeaks
    echo "$s: $(wc -l < $I/${s}_peaks.narrowPeak) peaks, $(wc -l < $I/${s}_nonpeaks_negatives.bed) background regions"
done
