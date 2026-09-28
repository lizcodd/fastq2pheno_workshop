# fastq2pheno: bulk ATAC-seq from FASTQ to biology

A worked example of a bulk ATAC-seq analysis on MIT's Engaging cluster: download the reads, run a standard
Nextflow pipeline, check QC, find differential peaks, then use deep-learning models (chromBPNet and
seq2PRINT) to ask which DNA sequences drive accessibility. Step 4 is in progress.

```
fastq2pheno_workshop/
  environment.yml   conda env with the tools used outside the pipeline
  01_data/          Step 1: sample table, download script
  02_pipeline/      Step 2: nf-core/atacseq
  03_differential/  Step 3: differential peaks (notebook)
  04_ai_models/     Step 4: sequence models of accessibility (notebook; in progress)
```

---

## Step 1: Get the data

### The dataset

**Wu et al. 2019, *Nature* 575:699, "Circular ecDNA promotes accessible chromatin and high oncogene
expression"** (GEO [GSE131956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131956)). COLO320 is a
colorectal cancer line with a large amplification of the *MYC* oncogene, in two versions from the same patient:

- **COLO320DM**: the amplicon is on **extrachromosomal DNA (ecDNA)**, small DNA circles outside the chromosomes;
- **COLO320HSR**: the same amplicon is **integrated into a chromosome**.

The paper reports that ecDNA is more accessible than chromosomal DNA. The workshop uses four of its ATAC-seq
runs, two replicates per line (paired-end, HiSeq 4000; also in `01_data/samples.tsv`):

| Run | Sample | Replicate | Read pairs |
|---|---|---|---|
| SRR8236757 | COLO320DM | 1 | 62.3 M |
| SRR8236758 | COLO320DM | 2 | 66.2 M |
| SRR9140502 | COLO320HSR | 1 | 37.5 M |
| SRR9140503 | COLO320HSR | 2 | 70.3 M |

**Limits of this dataset:**

1. **Line and batch are the same thing**: the DM runs were made and deposited together, the HSR runs later.
   Any DM vs HSR difference may partly be batch.
2. **Copy number differs** between the lines, and more copies of DNA give more reads.
3. **Two replicates per line**: enough to explore, not for strong claims.

### Download

Published reads live in **SRA** (NCBI) and its mirrors ENA and DDBJ; GEO links each sample to its runs
(`SRR...`). Download with the SRA Toolkit: `prefetch` fetches a run, `fasterq-dump` converts it to FASTQ, and
`pigz` compresses it. These tools and everything else outside the pipeline are in `environment.yml`:

```bash
module load miniforge/25.11.0-0
conda env create -f environment.yml     # once; creates an env called fastq2pheno
conda activate fastq2pheno
```

Conda envs hold tens of thousands of files, and the home directory allows only 1 million. Put envs in project
space by listing a project folder first under `envs_dirs` in `~/.condarc`:

```yaml
envs_dirs:
  - /orcd/data/<lab>/<you>/conda_envs
  - ~/.conda/envs
```

Then run the download as a batch job (next section):

```bash
cd 01_data
sbatch download.sbatch            # prints the job ID
tail download_<jobid>.log         # each run should end with "OK"
```

For each run, the script does:

```bash
prefetch --max-size u --output-directory sra SRR8236757
fasterq-dump --split-files --threads 8 --temp tmp sra/SRR8236757     # SRR8236757_1.fastq, _2.fastq
pigz -p 8 SRR8236757_1.fastq SRR8236757_2.fastq
```

**Check that every download is complete.** A truncated FASTQ can fail silently much later. The script counts
the reads in each file and compares them with the archive's count (`read_pairs` in `samples.tsv`), stopping
with `FAILED` if they differ. Once every run says OK, delete `sra/`.

**Your own data** from the MIT BioMicroCenter arrives in a shared folder, with a copy command like
`rsync -av /path/to/core/data /path/to/destination`. Run it in a batch job, into project space
(`/orcd/data/<lab>/...`), not home. Rerunning `rsync` only copies what is missing, so it is safe to repeat.

### Long jobs: submit them with sbatch

Downloads, copies and pipelines outlast a laptop's connection. Put the commands in a script with `#SBATCH`
lines requesting resources, and submit it; SLURM runs it on a compute node and writes a log:

```bash
#!/bin/bash
#SBATCH --job-name=copy_run
#SBATCH --partition=mit_normal     # up to 12 hours
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=06:00:00
#SBATCH --output=copy_%j.log       # %j becomes the job ID

rsync -av /path/to/core/data /path/to/destination
```

```bash
sbatch copy.sbatch       # submit
squeue --me              # queued and running jobs
tail -f copy_<jobid>.log # watch the log (Ctrl-C stops watching, not the job)
scancel <jobid>          # cancel
sacct -j <jobid>         # afterwards: COMPLETED or FAILED?
```

### Look inside the FASTQs

In `01_data/fastq/`, each run is a pair of files: `_1` has read 1 and `_2` read 2, from the two ends of each
DNA fragment, in the same order. Read them with `zcat`:

```bash
zcat SRR8236757_1.fastq.gz | head -4
```

```
@SRR8236757.1 1 length=51
NTATACAATAAGTTATCTAGAGCAACCATCAAAAAGCTATACAAAGGGATA
+SRR8236757.1 1 length=51
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
```

Every read is 4 lines: name, sequence (`N` = no base call), `+`, and one quality character per base (Phred score
= ASCII code − 33; `J` = 41, about 1 error in 10,000; `#` = 2, no confidence).

```bash
echo $(( $(zcat SRR8236757_1.fastq.gz | wc -l) / 4 ))                                  # number of reads
zcat SRR8236757_1.fastq.gz | head -400000 | awk 'NR % 4 == 2 {print length($0)}' \
    | sort -n | uniq -c | sort -nr | head -3                                            # read lengths
zcat SRR8236757_1.fastq.gz | head -400000 | grep -c CTGTCTCTTATACACATCT                  # Tn5 adapter left?
```

These reads are up to 51 bases, many shorter, and contain no Tn5 adapter: the authors trimmed them before
uploading. Reads straight from a sequencing core are untrimmed; the pipeline trims them either way.

---

## Step 2: Process the reads with nf-core/atacseq

**Nextflow** runs a chain of tools for you: it starts each step when its inputs are ready, runs samples in
parallel, submits each step as its own SLURM job, and with `-resume` restarts a failed run where it stopped.
**nf-core** is a collection of peer-reviewed Nextflow pipelines whose tools come in containers, so only
Nextflow needs installing. For bulk ATAC-seq: [nf-core/atacseq](https://nf-co.re/atacseq/2.1.2).

For each sample it trims adapters, aligns (BWA-MEM), marks duplicates, filters (mitochondrial reads, ENCODE's
blacklist, duplicates, multi-mapping reads), makes coverage tracks (bigWig) and calls peaks (MACS2). It then
builds consensus peaks across samples with a table of reads per peak per sample (the input to Step 3), and a
MultiQC report. It does all this once per replicate and once with each line's replicates merged. Its
"differential" DESeq2 step only normalizes counts for a PCA; the test is Step 3.

### Inputs

**Nextflow** is in the `fastq2pheno` env, and **apptainer** (runs the containers) is a module. Input files in
`02_pipeline/`:

`samplesheet.csv`: one row per FASTQ pair; rows with the same `sample` are replicates.

```
sample,fastq_1,fastq_2,replicate
COLO320DM,../01_data/fastq/SRR8236757_1.fastq.gz,../01_data/fastq/SRR8236757_2.fastq.gz,1
...
```

`params.yaml`: the settings.

```yaml
input:       samplesheet.csv
outdir:      results
blacklist:   hg38-blacklist.v3.bed   # regions with artifactual signal in every experiment
mito_name:   chrM                    # the chromosome to filter out
macs_gsize:  2.7e9                   # mappable genome size, for MACS2
narrow_peak: true                    # ATAC peaks are narrow (the default is broad)
```

`local.config`: the reference genome on your cluster, kept out of git. Copy `local.config.example` and fill in
the paths: a GRCh38 FASTA with UCSC chromosome names (`chr1`, `chrM`) and no alternate haplotypes, its BWA index,
and a gene annotation (GTF). Without reference files on hand, set `genome = 'GRCh38'` and `read_length = 50`
in its `params` block instead, and the pipeline downloads them from AWS iGenomes.

`run.sbatch`: the job that runs Nextflow.

```bash
module load miniforge/25.11.0-0 apptainer/1.5.2
conda activate fastq2pheno

export NXF_SINGULARITY_CACHEDIR=$PWD/containers   # container images, kept outside work/
export APPTAINER_CACHEDIR=$PWD/containers/tmp      # apptainer's download cache, out of the small home directory
export NXF_SYNTAX_PARSER=v1                        # nf-core/atacseq 2.1.2 needs Nextflow's older parser

nextflow run nf-core/atacseq -r 2.1.2 -profile engaging -params-file params.yaml -c local.config -resume
```

`-r 2.1.2` pins the pipeline version; `-profile engaging` runs each step as a SLURM job with apptainer. Without
`NXF_SYNTAX_PARSER=v1`, current Nextflow stops at once with `Unexpected input: '('`, because this 2023 pipeline's
config defines a function.

### Run it

```bash
cd 02_pipeline
sbatch run.sbatch
tail -f atacseq_<jobid>.log      # ends with "Pipeline completed successfully"
```

It takes about 2 hours for the four samples (187 steps). If it fails, fix the problem and submit again:
`-resume` skips finished steps. To test a setup first, `-profile test,engaging` runs a tiny built-in dataset.
Once the run is finished, delete `02_pipeline/work/` (about 100 GB of intermediate files).

### What comes out

```
results/
  multiqc/narrow_peak/multiqc_report.html     <- the QC report
  bwa/merged_library/                         <- per replicate ("mLb")
    COLO320DM_REP1.mLb.clN.sorted.bam         <- filtered alignments
    bigwig/                                   <- coverage tracks
    macs2/narrow_peak/consensus/
      consensus_peaks.mLb.clN.featureCounts.txt  <- reads per peak per sample (Step 3)
  bwa/merged_replicate/                       <- the same, replicates merged ("mRp")
```

### Look at the tracks in IGV

The pipeline's IGV session only opens on the server. `make_igv_local.sh` packages the coverage tracks and peaks
with a session that uses IGV's built-in hg38, to download and open in [IGV](https://igv.org/doc/desktop/) on
your computer (File → Open Session):

```bash
cd 02_pipeline
bash make_igv_local.sh          # writes igv_local/ and igv_local.tar (1.1 GB)
```

Coverage is scaled to reads per million, not corrected for copy number, so the *MYC* amplicon towers over
everything. Step 3 saves the differential peaks as BED files to drag into this session.

### Check quality with MultiQC

Open `results/multiqc/narrow_peak/multiqc_report.html` in a browser. Section names start with a level:
**LIB** (one FASTQ pair), **MERGED LIB** (one replicate) and **MERGED REP** (replicates merged). Here each
replicate is one FASTQ pair, so LIB and MERGED LIB repeat each other. Most of what matters is in MERGED LIB,
filtered:

1. **SAMtools / Picard**: reads aligned, duplicates, and what filtering removed.
2. **Picard insert size**: a nucleosome ladder, with a peak under 100 bp (open DNA) and smaller ones at ~200
   and ~400 bp (one and two nucleosomes).
3. **deepTools**: signal concentration (fingerprint) and reads peaking at transcription start sites.
4. **MACS2 peak count and FRiP** (fraction of reads in peaks): the headline quality measure.
5. **DESeq2 PCA**: replicates of the same line should sit together.

---

## Step 3: Find differential peaks

`03_differential/differential_atac.ipynb` starts from the pipeline's reads per peak: DESeq2, checks for
technical biases, copy number, whether ecDNA is more accessible per copy, motifs (HOMER), and BED files for IGV.
It runs in Colab (Runtime → Change runtime type → R) and downloads its data (12 MB, built by
`03_differential/prep/prepare_data.ipynb`); on Engaging, use the `R (fastq2pheno)` kernel. About 5 minutes.

Things that change read counts besides biology, and how to check them:

| Bias | How to check |
|---|---|
| **Batch** | Know how samples were made; put every condition in every batch. |
| **Copy number** | Copy number (WGS, or reads between peaks) and fold changes along the genome. |
| **Depth, signal-to-noise** | MultiQC (reads, duplicates, FRiP); fold change vs peak strength (MA plot). |
| **GC content** | Fold change vs peak GC; if it trends, check reads between peaks. Correct (e.g. `cqn`) only if technical. |
| **Mappability** | Compare only libraries with the same read length; be careful in repeats. |

Before reading any list of top peaks, plot fold changes along the genome and against GC and peak strength.

---

## Step 4: Sequence models of accessibility

> **In progress.** The models are trained; their evaluation, motif results and the DM vs HSR comparison are still
> being added.

chromBPNet and seq2PRINT: setup, runs and outputs are in [04_ai_models/README.md](04_ai_models/README.md);
the notebook is `04_ai_models/ai_models.ipynb`.
