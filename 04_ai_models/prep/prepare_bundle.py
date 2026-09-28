# Build the data bundle for ai_models.ipynb (runs on Engaging, in the scprinter env; CPU only).
# Three stages, run in order by run_prep.sbatch:
#   python prepare_bundle.py select   picks the featured peaks -> work/featured.narrowPeak
#   (chrombpnet contribs_bw for every chromBPNet model on those peaks -> work/contribs/, in the container)
#   python prepare_bundle.py build    everything else, and the archive
#
# Writes 04_ai_models/data/ and 04_ai_models/ai_models_data.tar.gz:
#   regions.npz            featured peaks: 2,114 bp sequence + per-base tracks over the central 1,000 bp
#   motifs_chrombpnet.json TF-MoDISco motifs of each line's chromBPNet model: logo matrix, seqlets, matches
#   qc_chrombpnet.tsv      test-chromosome accuracy of each chromBPNet model
#   models/                chromBPNet (.h5) and seq2PRINT (TorchScript .pt) models for live prediction
#
# One model per tool and line: both replicates pooled, fold 0 (see ../README.md). Evaluation metrics and TF-MoDISco
# results are optional: while missing, qc_chrombpnet.tsv / motifs_chrombpnet.json are left out and the notebook skips them.
import glob, json, os, re, shutil, tarfile
import numpy as np, pandas as pd, pyBigWig, pyfaidx, h5py, torch

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))     # the workshop folder
RESULTS = f"{ROOT}/04_ai_models/results"
OUT = f"{ROOT}/04_ai_models/data"
WORK = f"{ROOT}/04_ai_models/prep/work"
LINES = ["DM", "HSR"]

CBP, S2P = f"{RESULTS}/chrombpnet", f"{RESULTS}/seq2print"
MODELS = {
    "chrombpnet": {line: [dict(model=f"{CBP}/models/COLO320{line}_fold0/models/chrombpnet_nobias.h5",
                               metrics=f"{CBP}/models/COLO320{line}_fold0/evaluation/chrombpnet_metrics.json",
                               modisco=f"{CBP}/contribs/COLO320{line}_fold0",
                               name=f"COLO320{line}_fold0")]
                   for line in LINES},
    "seq2print": {line: [dict(model=f"{S2P}/model/COLO320_COLO320{line}_fold0-COLO320{line}_fold0.pt",
                              attr=f"{S2P}/model/COLO320_COLO320{line}_fold0-COLO320{line}_fold0.pt_COLO320_COLO320{line}/attr.count.shap_hypo_0_.0.85.bigwig",
                              tfbs=f"{S2P}/COLO320{line}_TFBS.bigwig", name=f"COLO320{line}_fold0")]
                  for line in LINES},
}
COVERAGE = {line: f"{ROOT}/02_pipeline/results/bwa/merged_replicate/bigwig/COLO320{line}.mRp.clN.bigWig" for line in LINES}
GENOME = f"{ROOT}/02_pipeline/results/genome/genome.fa"
CALLS = f"{WORK}/calls.tsv"                    # DESeq2 calls as in Step 3, from calls.R
TFBS_REGIONS = f"{S2P}/tfbs_regions.bed"       # the peaks seq2PRINT scored (center +- 250 bp)
N_PER_CLASS = 250
HALF_IN, HALF_TRACK = 1057, 500               # 2,114 bp model input; 1,000 bp of tracks


def select_regions():
    calls = pd.read_csv(CALLS, sep="\t")
    peaks = pd.read_csv(f"{ROOT}/03_differential/data/peaks.tsv.gz", sep="\t",
                        usecols=["peak", "gene", "COLO320DM_REP1", "COLO320DM_REP2", "COLO320HSR_REP1", "COLO320HSR_REP2"])
    calls = calls.merge(peaks, on="peak")
    calls["mean_count"] = calls[[c for c in calls if c.startswith("COLO320")]].mean(axis=1)
    calls["center"] = (calls.start + calls.end) // 2
    tfbs = pd.read_csv(TFBS_REGIONS, sep="\t", header=None, names=["chr", "s", "e"])
    scored = set(zip(tfbs.chr, tfbs.s + 250))
    calls = calls[[(c, x) in scored for c, x in zip(calls.chr, calls.center)]]
    calls = calls[calls.region_class.str.startswith("off-amplicon")]
    up_dm = calls[calls.call == "up in DM"].nsmallest(N_PER_CLASS, "padj")
    up_hsr = calls[calls.call == "up in HSR"].nsmallest(N_PER_CLASS, "padj")
    same = calls[(calls.call == "not significant") & (calls.mean_count >= calls.mean_count.quantile(0.75))]
    same = same.sample(N_PER_CLASS, random_state=0)
    sel = pd.concat([up_dm, up_hsr, same]).reset_index(drop=True)
    print(f"{len(calls)} candidate peaks; selected {len(sel)}:", sel.call.value_counts().to_dict())
    os.makedirs(WORK, exist_ok=True)
    sel.to_csv(f"{WORK}/featured.tsv", sep="\t", index=False)
    # chromBPNet region format: 10-column narrowPeak, column 10 = summit offset from the start
    pd.DataFrame({"chr": sel.chr, "s": sel.center - 500, "e": sel.center + 500, "name": sel.peak, "score": 0,
                  "strand": ".", "a": 0, "b": 0, "c": 0, "summit": 500}).to_csv(
        f"{WORK}/featured.narrowPeak", sep="\t", header=False, index=False)
    # the chromBPNet models to interpret, one per line: "<line> <k> <model path>"
    with open(f"{WORK}/chrombpnet_models.txt", "w") as f:
        for line in LINES:
            for k, m in enumerate(MODELS["chrombpnet"][line]):
                f.write(f"{line} {k} {m['model']}\n")
    return sel


def track(bw_paths, chrom, start, end):
    """Mean over bigwigs of per-base values (missing = 0)."""
    vals = []
    for p in bw_paths:
        bw = pyBigWig.open(p)
        v = np.array(bw.values(chrom, start, end), dtype=np.float32)
        bw.close()
        vals.append(np.nan_to_num(v))
    return np.mean(vals, axis=0)


def build_regions(sel):
    fa = pyfaidx.Fasta(GENOME)
    names, arrays = [], []
    for line in LINES:
        names += [f"coverage_{line}", f"chrombpnet_contrib_{line}", f"seq2print_contrib_{line}", f"seq2print_tfbs_{line}"]
    seqs = []
    out = np.zeros((len(sel), len(names), 2 * HALF_TRACK), dtype=np.float16)
    for i, r in sel.iterrows():
        c = int(r.center)
        seqs.append(str(fa[r.chr][c - HALF_IN:c + HALF_IN]).upper())
        s, e = c - HALF_TRACK, c + HALF_TRACK
        k = 0
        for line in LINES:
            out[i, k] = track([COVERAGE[line]], r.chr, s, e)
            out[i, k + 1] = track(sorted(glob.glob(f"{WORK}/contribs/chrombpnet_{line}_*.counts_scores.bw")), r.chr, s, e)
            out[i, k + 2] = track([m["attr"] for m in MODELS["seq2print"][line]], r.chr, s, e)
            out[i, k + 3] = track([m["tfbs"] for m in MODELS["seq2print"][line]], r.chr, s, e)
            k += 4
    np.savez_compressed(f"{OUT}/regions.npz", tracks=out, track_names=np.array(names), seq=np.array(seqs),
                        chrom=sel.chr.values, center=sel.center.values, peak=sel.peak.values,
                        gene=sel.gene.fillna("").values, call=sel.call.values,
                        log2FC=sel.log2FC.values, padj=sel.padj.values)
    print("regions.npz:", out.shape)


def parse_report(html):
    rows = re.findall(r"<tr>(.*?)</tr>", open(html).read(), flags=re.S)
    table = {}
    for row in rows:
        cells = re.findall(r"<td>(.*?)</td>", row, flags=re.S)
        if not cells:
            continue
        matches = [(cells[i], cells[i + 1]) for i in (4, 7, 10) if i + 1 < len(cells) and cells[i] not in ("", "NaN", "nan")]
        table[cells[0]] = matches
    return table


def trim(cwm, frac=0.3):
    score = np.abs(cwm).sum(axis=1)
    keep = np.where(score >= frac * score.max())[0]
    return cwm[max(keep.min() - 2, 0):keep.max() + 3]


def build_motifs():
    out = {}
    for line in LINES:
        m = MODELS["chrombpnet"][line][0]
        if not os.path.exists(f"{m['modisco']}/modisco_counts_report/motifs.html"):
            print(f"TF-MoDISco for {m['name']} not finished: motifs_chrombpnet.json left out")
            return
        matches = parse_report(f"{m['modisco']}/modisco_counts_report/motifs.html")
        pats = []
        with h5py.File(f"{m['modisco']}/modisco_counts.h5") as f:
            for sign in ("pos_patterns", "neg_patterns"):
                if sign not in f:
                    continue
                for p in sorted(f[sign], key=lambda x: int(x.split("_")[1])):
                    g = f[sign][p]
                    pats.append(dict(pattern=f"{sign}.{p}", n_seqlets=int(g["seqlets"]["n_seqlets"][()][0]),
                                     cwm=trim(g["contrib_scores"][()]).round(4).tolist(),
                                     matches=matches.get(f"{sign}.{p}", [])))
        out[line] = pats
    json.dump(out, open(f"{OUT}/motifs_chrombpnet.json", "w"))
    print("motifs:", {k: len(v) for k, v in out.items()})


def build_qc():
    rows = []
    for line in LINES:
        for m in MODELS["chrombpnet"][line]:
            if not os.path.exists(m["metrics"]):
                print(f"evaluation of {m['name']} not finished: qc_chrombpnet.tsv left out")
                return
            q = json.load(open(m["metrics"]))
            rows.append(dict(line=line, model=m["name"], counts_pearson_r=q["counts_metrics"]["peaks"]["pearsonr"],
                             counts_spearman_r=q["counts_metrics"]["peaks"]["spearmanr"],
                             profile_median_jsd=q["profile_metrics"]["peaks"]["median_jsd"],
                             profile_median_norm_jsd=q["profile_metrics"]["peaks"]["median_norm_jsd"]))
    pd.DataFrame(rows).round(3).to_csv(f"{OUT}/qc_chrombpnet.tsv", sep="\t", index=False)


def build_models():
    os.makedirs(f"{OUT}/models", exist_ok=True)
    for line in LINES:
        for k, m in enumerate(MODELS["chrombpnet"][line]):
            shutil.copy(m["model"], f"{OUT}/models/chrombpnet_{line}_{k}.h5")
        for k, m in enumerate(MODELS["seq2print"][line]):
            net = torch.load(m["model"], map_location="cpu", weights_only=False).eval()
            x = torch.zeros(1, 4, net.dna_len); x[:, 0] = 1
            with torch.no_grad():
                torch.jit.trace(net, x, check_trace=False).save(f"{OUT}/models/seq2print_{line}_{k}.pt")
    print("models:", sorted(os.listdir(f"{OUT}/models")))


if __name__ == "__main__":
    import sys
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "select":
        shutil.rmtree(OUT); os.makedirs(OUT)          # start clean, so no file from older models survives
        select_regions()
        raise SystemExit
    build_regions(pd.read_csv(f"{WORK}/featured.tsv", sep="\t"))
    build_motifs()
    build_qc()
    build_models()
    with tarfile.open(f"{ROOT}/04_ai_models/ai_models_data.tar.gz", "w:gz") as t:
        t.add(OUT, arcname="data")
    print("archive MB:", round(os.path.getsize(f"{ROOT}/04_ai_models/ai_models_data.tar.gz") / 1e6, 1))
