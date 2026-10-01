# Every consensus peak on chr1-22 and X with Step 3's call: 2,114 bp around its center, for predict.py.
#   work/peaks.tsv, work/peak_seqs.txt
# Runs in the chromBPNet container (CPU).
import os
import pandas as pd, pyfaidx

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))     # the workshop folder
W = f"{ROOT}/04_sequence_models/prep/work"
GENOME = f"{ROOT}/02_pipeline/results/genome/genome.fa"
HALF = 1057


def peaks(fa):
    calls = pd.read_csv(f"{ROOT}/04_sequence_models/prep/work/calls.tsv", sep="\t")
    gc = pd.read_csv(f"{ROOT}/03_differential/data/peaks.tsv.gz", sep="\t", usecols=["peak", "gc", "feature"])
    p = calls.merge(gc, on="peak")
    p["center"] = (p.start + p.end) // 2
    p = p[p.chr.isin([f"chr{c}" for c in list(range(1, 23)) + ["X"]])].reset_index(drop=True)
    with open(f"{W}/peak_seqs.txt", "w") as f:
        for c, x in zip(p.chr, p.center):
            f.write(str(fa[c][x - HALF:x + HALF]).upper() + "\n")
    p.to_csv(f"{W}/peaks.tsv", sep="\t", index=False)
    print("peaks:", len(p))


if __name__ == "__main__":
    os.makedirs(W, exist_ok=True)
    peaks(pyfaidx.Fasta(GENOME))
