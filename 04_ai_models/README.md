# Step 4: Sequence models of accessibility (chromBPNet and seq2PRINT)

> **In progress.** The models are trained; their evaluation, motif results and the DM vs HSR comparison are still
> being added (see the TODOs below).

Step 3 asked which peaks differ and which known motifs they contain. Here two deep-learning models ask **which
DNA sequences make a region accessible in each line**. Both read DNA sequence and predict the ATAC signal; the
bases that drive each prediction spell out the TF motifs the model learned from the data.

```
04_ai_models/
  ai_models.ipynb   the Colab notebook: explore the trained models (no GPU needed)
  chrombpnet/       scripts to train, evaluate and interpret chromBPNet
  seq2print/        scripts to train and interpret seq2PRINT (scPrinter)
  prep/             builds the data bundle the notebook downloads
  results/          what the scripts write (not in git; about 26 GB)
  data/             the notebook's data bundle (not in git; packed as ai_models_data.tar.gz)
```

Training takes days of GPU time, so the notebook starts from trained models. The sections below show how to
train them.

## The two tools

| | chromBPNet | seq2PRINT |
|---|---|---|
| Reads | 2,114 bp of sequence | ~1,800 bp of sequence |
| Predicts | reads per region (**counts**) and where Tn5 cuts inside it (**profile**) | **footprints** at many widths (TF- to nucleosome-sized), plus reads per region |
| Tn5 bias | a pretrained **bias model**, plugged in during training | removed before training, with a genome-wide Tn5 track |
| Outputs | QC report, contribution scores, TF-MoDISco motifs | contribution scores, **TF binding scores** (0–1 per base) |
| Source | [kundajelab/chrombpnet](https://github.com/kundajelab/chrombpnet) | [buenrostrolab/scPrinter](https://github.com/buenrostrolab/scPrinter) |

- **Footprint**: where a protein sits on DNA, Tn5 can't cut, leaving a dip in the cut pattern (~20 bp for a TF,
  ~150 bp for a nucleosome).
- **Contribution score**: how much each base pushes a model's prediction up or down. Drawn as letters scaled by
  their contribution, a bound motif stands out as a block of tall letters.
- **TF binding score**: a second, pretrained model's probability that some TF is bound at each position.

**Models here**: one per line and tool, trained on **both replicates pooled**. Each model trains on most
chromosomes, stops training based on two (validation), and is tested on chromosomes it never saw; this split
is a **fold**. All models use chromBPNet's **fold 0** (test chr1, 3, 6; validation chr8, 20). ENCODE averages
five folds, which costs five times the GPU time.

---

## chromBPNet

**Needs**: one L40S GPU per model (about 12 GPU-hours, longer than `mit_normal_gpu`'s 6-hour limit, so the GPU
scripts need a SLURM QOS that allows longer jobs, added when submitting), and the chromBPNet container:

```bash
cd 04_ai_models/chrombpnet
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
  fold-0 Tn5 bias model for K562 ATAC-seq (file ENCFF984RAF, checked against its MD5), then makes each line's
  inputs: the merged-replicate BAM on chr1–22 and X, its MACS2 peaks, and GC-matched background regions.
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

**TF-MoDISco vs HOMER.** HOMER (Step 3) compares two peak sets against a motif database. TF-MoDISco describes
one model: it clusters the stretches of sequence the model found important (**seqlets**) into motifs, then
matches them to databases, so it also finds variants and fixed-spacing motif pairs. Its top motifs describe
what drives accessibility in that line overall, not what differs between lines. Check the logo, not only the
database name.

TODO: QC numbers and motif results (evaluation and TF-MoDISco still running).

---

## seq2PRINT (scPrinter)

**PRINT** measures footprints at many widths at once, corrected for Tn5's sequence preference. **seq2PRINT**
trains a sequence model to predict them, which removes much of the noise of low read counts, and a second model
(trained by the authors on ChIP-seq) turns them into TF binding scores. scPrinter was built for single cells;
on bulk data each sample is one "big cell".

**Install.** scPrinter's dependencies have moved on since its last release (v1.2.0), so a plain install fails.
`build_env.sh` pins the versions that work together (tangermeme 0.4.4, snapATAC2 2.8.0 with numpy 1.26, pyBigWig
from pip, MACS2, plus two undeclared dependencies) and keeps scPrinter's downloads (about 16 GB) in
`results/scprinter_data/`:

```bash
cd 04_ai_models/seq2print
bash build_env.sh          # creates the env "scprinter"; about 30 minutes
```

**Run**:

```bash
sbatch fragments.sbatch    # fragment files, one per replicate (CPU)
sbatch prepare.sbatch      # import them and call training peaks (CPU; after ../chrombpnet/prepare.sbatch)
sbatch configs.sbatch      # one training configuration per line, replicates pooled, fold 0
sbatch train.sbatch        # one GPU job per line, 7-9 h each
sbatch tfbs.sbatch         # TF binding scores, after train.sbatch (GPU, about 4 h each)
```

- **Fragment files**: one line per DNA fragment (chromosome, start, end, sample), with ends at the Tn5 cut
  sites (+4 / −5 bp from the read ends).
- **Peaks**: scPrinter's `seq2PRINT` peak setting, minus the blacklist and amplicons. scPrinter imports files
  with worker processes that re-run the script, so script code must sit under `if __name__ == "__main__":`.
- **Training** runs up to 300 epochs (stopping after 5 without improvement), then computes contribution scores.
  A job killed at its time limit loses those, and a restarted job starts again from the first epoch.

**Outputs**:

```
results/seq2print/
  model/COLO320_<line>_fold0-<line>_fold0.pt               <- the model
  model/<model>.pt_COLO320_<line>/attr.count.shap_hypo_0_.0.85.bigwig   <- contribution scores
  <line>_TFBS.bigwig                                        <- TF binding scores
```

---

## Comparing DM and HSR

Each model describes one line, and TF-MoDISco motifs are clustered separately per model, so they can't be
compared directly. Two comparisons on a common scale:

- **Motif insertion**: insert a motif into background sequences and compare how much each model's predicted
  accessibility rises (`chrombpnet footprints`, scPrinter's `delta_effects_seq2print`).
- **The same genomic sites** scored with the DM and HSR models.

TODO: scripts and results.

---

## The notebook's data

```bash
cd 04_ai_models/prep
sbatch run_prep.sbatch     # CPU
```

`run_prep.sbatch` recomputes Step 3's differential calls (`calls.R`), picks 750 featured peaks (the 250 most
significant up in DM and up in HSR, off the amplicons, and 250 strong unchanged ones), computes chromBPNet
contribution scores on them, and packs sequences, tracks, both tools' models, QC numbers and motifs into
`data/` and `ai_models_data.tar.gz`. Until the evaluation and TF-MoDISco finish, the bundle is built without
them and the notebook skips those sections.
