# seq2PRINT TF binding scores for one line (following scPrinter's PBMC_bulkATAC_tutorial).
# Regions: pipeline consensus peaks, 500 bp around each peak center, chr1-22/X, > 200 reads across the 4 samples.
# Usage: python tfbs.py COLO320DM.  Output: ../results/seq2print/<line>_TFBS.bigwig
import os, sys, json
import pandas as pd
import scprinter as scp

SAMPLES = ["COLO320DM_REP1", "COLO320DM_REP2", "COLO320HSR_REP1", "COLO320HSR_REP2"]
WORK = os.path.abspath("../results/seq2print")


def make_regions(path):
    fc = pd.read_csv("../../02_pipeline/results/bwa/merged_library/macs2/narrow_peak/consensus/"
                     "consensus_peaks.mLb.clN.featureCounts.txt", sep="\t", comment="#")
    main = [f"chr{c}" for c in list(range(1, 23)) + ["X"]]
    fc = fc[fc["Chr"].isin(main) & (fc.iloc[:, 6:10].sum(axis=1) > 200)]
    center = ((fc["Start"] + fc["End"]) // 2).astype(int)
    regions = pd.DataFrame({"chr": fc["Chr"], "start": center - 250, "end": center + 250})
    regions.sort_values(["chr", "start"]).to_csv(path, sep="\t", header=False, index=False)
    print(len(regions), "regions written to", path)


def main():
    sample = sys.argv[1]            # a line: COLO320DM or COLO320HSR
    regions = f"{WORK}/tfbs_regions.bed"
    if not os.path.exists(regions):
        make_regions(regions)
    fold = 0
    config = json.load(open(f"{WORK}/configs/COLO320_{sample}_fold{fold}.JSON"))
    scp.tl.seq_tfbs_seq2print(
        seq_attr_count=None, seq_attr_footprint=None, genome=scp.genome.hg38,
        region_path=regions, gpus=[0], model_type="seq2print",
        model_path=f"{WORK}/model/COLO320_{sample}_fold{fold}-{sample}_fold{fold}.pt",
        lora_config=config, group_names=[sample], verbose=True, launch=True,
        return_adata=False, overwrite_seqattr=True, post_normalize=True,
        save_key=f"COLO320_{sample}", save_path=WORK)


if __name__ == "__main__":
    main()
