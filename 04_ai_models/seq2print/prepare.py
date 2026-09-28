# Import the four replicate fragment files into scPrinter and call the seq2PRINT training peaks
# (following scPrinter's PBMC_bulkATAC_tutorial), then drop peaks in amplicons / blacklist.
import os, subprocess
import pandas as pd
import scprinter as scp

WORK = "../results/seq2print"


# scPrinter imports multiple fragment files with 'spawn' worker processes, which re-import this script:
# everything must run under the __main__ guard.
def main():
    samples = ["COLO320DM_REP1", "COLO320DM_REP2", "COLO320HSR_REP1", "COLO320HSR_REP2"]
    frags = [f"{WORK}/fragments/{s}.frags.tsv.gz" for s in samples]
    os.makedirs(WORK, exist_ok=True)

    if os.path.exists(f"{WORK}/colo320_scprinter.h5ad"):          # already imported: reuse it
        printer = scp.load_printer(f"{WORK}/colo320_scprinter.h5ad", scp.genome.hg38)
    else:
        printer = scp.pp.import_fragments(
            path_to_frags=frags, barcodes=[None] * len(frags), sample_names=samples,
            savename=f"{WORK}/colo320_scprinter.h5ad", genome=scp.genome.hg38,
            min_num_fragments=1000, min_tsse=7, sorted_by_barcode=False, low_memory=False)
    print("samples:", list(printer.insertion_file.obs_names))
    printer.insertion_file.obs_names = samples

    # seq2PRINT preset: the peak set the authors recommend for training
    scp.pp.call_peaks(printer=printer, frag_file=frags, cell_grouping=[None], group_names=["all"],
                      preset="seq2PRINT", overwrite=False)
    peaks = pd.DataFrame(printer.uns["peak_calling"]["all_cleaned"][:])
    peaks.to_csv(f"{WORK}/seq2print_peaks_all.bed", sep="\t", header=False, index=False)
    printer.close()

    # drop peaks near amplicons or the blacklist (same exclusion as chromBPNet, including its 1,057 bp margin)
    excl = "../results/chrombpnet/inputs/exclude_ext.bed"      # made by ../chrombpnet/prepare.sbatch
    out = subprocess.run(f"bedtools intersect -v -a {WORK}/seq2print_peaks_all.bed -b {excl} | sort -k1,1 -k2,2n",
                         shell=True, capture_output=True, text=True, check=True).stdout
    open(f"{WORK}/seq2print_peaks.bed", "w").write(out)
    print("peaks:", len(peaks), "called;", out.count("\n"), "kept after excluding amplicons/blacklist")


if __name__ == "__main__":
    main()
