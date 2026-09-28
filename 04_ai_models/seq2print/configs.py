# One seq2PRINT training config per line, both replicates pooled, fold 0 (the same input as the chromBPNet models),
# following scPrinter's PBMC_bulkATAC_tutorial. Paths stay absolute (no path_swap): the training script joins them
# onto --data_dir with os.path.join, which leaves absolute paths unchanged.
import os
import scprinter as scp

WORK = os.path.abspath("../results/seq2print")
LINES = {"COLO320DM": ["COLO320DM_REP1", "COLO320DM_REP2"], "COLO320HSR": ["COLO320HSR_REP1", "COLO320HSR_REP2"]}


def main():
    os.makedirs(f"{WORK}/configs", exist_ok=True)
    printer = scp.load_printer(f"{WORK}/colo320_scprinter.h5ad", scp.genome.hg38)
    for line, reps in LINES.items():
        scp.tl.seq_model_config(
            printer, region_path=f"{WORK}/seq2print_peaks.bed", cell_grouping=[reps], group_names=[line],
            genome=printer.genome, fold=0, overwrite_bigwig=False, model_name=f"COLO320_{line}",
            additional_config={"tags": ["COLO320", line, "pooled", "fold0"]},
            config_save_path=f"{WORK}/configs/COLO320_{line}_fold0.JSON")
        print("config written for", line, reps)
    printer.close()


if __name__ == "__main__":
    main()
