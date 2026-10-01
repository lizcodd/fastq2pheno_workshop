#!/bin/bash
#SBATCH --job-name=cbp_prepare
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/prepare_%j.log
# Prepare chromBPNet inputs from the nf-core/atacseq results (Step 2), in ../results/chrombpnet/inputs: each line's
# reads, and the peaks and background regions both lines' models train on.
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

# Reads (both replicates merged) of each line on chr1-22 and X only: chromBPNet can't handle the small unplaced contigs
for s in COLO320DM COLO320HSR; do
    samtools view -@ 8 -b -o $I/$s.main.bam $RES/bwa/merged_replicate/$s.mRp.clN.sorted.bam $(cut -f1 $I/chrom.sizes)
    samtools index $I/$s.main.bam
done

# Peaks: both lines' models train on the same regions, the consensus peaks of all four samples (Step 2), so the two
# models differ only in the reads they learn from. Kept when the 2,114 bp model window stays clear of excluded regions
# and chromosome ends; consensus peaks have no summit, so the model window is centered on the middle of the peak.
grep -P '^chr([0-9]+|X)\t' $RES/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.bed \
  | bedtools intersect -v -a - -b $I/exclude_ext.bed \
  | awk -v OFS='\t' -v sizes=$I/chrom.sizes 'BEGIN { while ((getline < sizes) > 0) L[$1] = $2 }
      { s = int(($3 - $2) / 2); c = $2 + s; if (c - 1057 > 0 && c + 1057 < L[$1]) print $1, $2, $3, $4, 0, ".", 0, 0, 0, s }' \
  | sort -k1,1 -k2,2n > $I/consensus_peaks.narrowPeak
# Background regions matched to the peaks' GC content (fold 0)
cb chrombpnet prep nonpeaks -g $I/genome.fa -p $I/consensus_peaks.narrowPeak -c $I/chrom.sizes \
    -fl $I/folds/fold_0.json -br $I/exclude.bed -o $I/consensus_nonpeaks
echo "$(wc -l < $I/consensus_peaks.narrowPeak) peaks, $(wc -l < $I/consensus_nonpeaks_negatives.bed) background regions"
