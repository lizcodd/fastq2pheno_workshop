# Step 4: What makes DNA accessible? Sequence models with chromBPNet

Step 3 found the peaks that differ between DM and HSR, and HOMER found the AP-1 motif enriched in the peaks more
open in HSR. Enrichment can't say whether those sites cause the opening, or which sites, in which peaks.
[chromBPNet](https://github.com/kundajelab/chrombpnet) learns how sequence sets accessibility: trained on one line,
it reads DNA sequence and predicts that line's ATAC-seq signal. With one model per line, both can be given the same
sequence, with a motif removed or added, and their answers compared, down to single bases.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lizcodd/fastq2pheno_workshop/blob/main/04_sequence_models/chrombpnet.ipynb)

```
04_sequence_models/
  chrombpnet.ipynb   the Colab notebook: explore the trained models (no GPU needed)
  chrombpnet/        scripts to train, evaluate and interpret the models
  prep/              builds the data bundle the notebook downloads
  results/           what the scripts write (not in git; about 16 GB)
  data/              the notebook's data bundle (not in git; packed as chrombpnet_data.tar.gz)
```

Training takes hours of GPU time, so the notebook starts from the trained models and runs them live where that's
quick. The sections below show how to train them.

## What chromBPNet does

chromBPNet reads **2,114 bp of sequence** and predicts, for the central 1,000 bp, how many reads the region gets
(**counts**) and where exactly Tn5 cuts (**profile**, per base). Tn5's own sequence preference is handled by a
separate, pretrained **bias model** (here ENCODE's, from K562 ATAC-seq); the main model learns everything else on top
of it, and after training the main model alone (**bias-free**) is used for interpretation.

- **Contribution score**: how much each base pushes the prediction up or down. Drawn as letters scaled by their
  contribution, a bound motif stands out as a block of tall letters.
- **TF-MoDISco**: collects the stretches with high contribution across all peaks (**seqlets**), clusters them into
  motifs, then matches each to databases of known motifs.

**Our models.** One per line, trained on **both replicates pooled**. Both are trained on the **same regions**, the
consensus peaks of all four samples plus GC-matched background regions, so the two models differ only in the reads
they learn from. (chromBPNet's usual recipe trains each model on its own sample's peaks. For comparing two
conditions that is worse: each model then rarely sees the regions open only in the other line, and here it missed
much of the DM vs HSR difference.) Each model trains on most chromosomes, stops training based on two (validation), and is tested
on chromosomes it never saw; this split is a **fold**. The models use chromBPNet's **fold 0** (test chr1, 3, 6;
validation chr8, 20). ENCODE averages five folds, which costs five times the GPU time.

---

## Train the models

**Needs**: one L40S GPU per model (about 12 GPU-hours, longer than `mit_normal_gpu`'s 6-hour limit, so the GPU
scripts need a SLURM QOS that allows longer jobs, added when submitting), and the chromBPNet container:

```bash
cd 04_sequence_models/chrombpnet
module load apptainer/1.5.2
apptainer pull chrombpnet_latest.sif docker://kundajelab/chrombpnet:latest     # about 5 GB
```

**Run**:

```bash
sbatch prepare.sbatch     # inputs (CPU)
sbatch model.sbatch       # one GPU job per line: train, evaluate, contribution scores
sbatch modisco.sbatch     # TF-MoDISco motifs, after model.sbatch (CPU, 3-6 h)
```

- `prepare.sbatch` downloads the fold definitions ([Zenodo](https://zenodo.org/records/7443683)) and ENCODE's
  fold-0 Tn5 bias model for K562 ATAC-seq (file ENCFF984RAF, checked against its MD5), then makes the
  inputs: each line's merged-replicate BAM on chr1–22 and X, the consensus peaks (Step 2), and GC-matched
  background regions.
- Peaks near the ENCODE blacklist or an amplicon are left out. The model predicts reads from sequence and
  can't know that a region has 100 copies, so amplified regions would distort it; `amplicons.bed` lists the
  amplicons of both lines.
- `model.sbatch` trains the model on top of the bias model, evaluates it on the test chromosomes, then computes
  per-base contribution scores on every peak, for **counts** (which bases make a region accessible) and
  **profile** (which bases shape the cut pattern; TFs with strong footprints, like CTCF, stand out).
- These steps run separately, and each is skipped when its output exists, to guard against GPU time limits: a
  job stopped at its limit continues where it left off when resubmitted. Where time limits aren't an issue, the
  single command `chrombpnet pipeline` trains, evaluates and finds motifs (profile contributions on 30,000 peaks)
  in one job.

**Outputs**:

```
results/chrombpnet/models/<line>_fold0/
  models/chrombpnet_nobias.h5      <- the model with Tn5 bias removed: use this for biology
  evaluation/overall_report.html   <- the QC report
results/chrombpnet/contribs/<line>_fold0/
  *.counts_scores.bw, *.profile_scores.bw           <- contribution tracks (open in IGV)
  modisco_counts_report/, modisco_profile_report/    <- motifs
```

**Reading the QC report**:

1. **Tn5 check** (`bw_shift_qc.png`): the motif should be Tn5's symmetric insertion preference; if not, reads
   were shifted wrongly.
2. **Bias model alone** (`bias_metrics.json`): should predict reads per peak barely at all; chromBPNet warns
   below −0.3.
3. **Test accuracy** (`chrombpnet_metrics.json`): counts Pearson r of 0.6–0.8 is typical for bulk ATAC; the
   normalized profile JSD runs from 0 (no better than flat) to 1 (perfect) and depends on depth.
4. **Bias removed?** (`chrombpnet_nobias.tn5_*.footprint.png`): the bias-free model should barely react to
   Tn5's preferred sequence.

**TF-MoDISco vs HOMER.** HOMER (Step 3) finds motifs over-represented in one peak set compared with another.
TF-MoDISco describes one model: the motifs it relies on to predict accessibility, including motifs common
everywhere (like CTCF), variants and fixed-spacing motif pairs, and it scores every individual site. Its top motifs
describe what drives accessibility in that line overall, not what differs between lines. Check the logo, not only
the database name.

---

## Comparing DM and HSR

Each model describes one line. Their TF-MoDISco lists come from the same peaks, but a motif missing from one list
may just have fallen below TF-MoDISco's reporting threshold. The notebook compares the lines by giving **the same
sequence** to both models:

- **Motif sites in peaks**: scramble the site, predict again with both models, and compare the change. One peak is
  done live; every AP-1 and CTCF site near a peak's center is precomputed (`prep/ablation.py`).
- **Planted motifs**: insert a motif into background sequences and compare how much each model's prediction rises
  (live in the notebook).

---

## The notebook's data

```bash
cd 04_sequence_models/prep
sbatch run_prep.sbatch     # CPU, about 3 hours
```

`run_prep.sbatch` recomputes Step 3's differential calls (`calls.R`), picks 750 featured peaks (the 250 most
significant up in DM and up in HSR, and 250 strong unchanged ones, all off the amplicons with equal copy number,
each centered on its coverage summit), computes both models' contribution scores on them, cuts out every peak's sequence (`peaks.py`), scrambles each motif site near a peak's
center plus a random 6 bp control (`ablation.py`), predicts all of these with both models (`predict.py`; minutes on
a GPU, hours on a CPU), and packs everything with the two models into
`data/` and `chrombpnet_data.tar.gz` (`prepare_bundle.py`). Steps whose output already exists are skipped.

## The notebook

1. What chromBPNet is: the network and the Tn5 bias model
2. Did the models learn? (training curves, test accuracy, bias checks)
3. What each model relies on: its TF-MoDISco motifs
4. One peak, base by base (ranking tables, contribution scores)
5. Edit the DNA: scramble a site and plant a motif (live), then every AP-1 site at once
6. Takeaways
