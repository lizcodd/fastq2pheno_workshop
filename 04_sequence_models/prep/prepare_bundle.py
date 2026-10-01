# Build the data bundle for chrombpnet.ipynb. Runs in the chromBPNet container (CPU is enough); run_prep.sbatch runs
# the stages in order:
#   python prepare_bundle.py select   featured peaks, motif-planting backgrounds
#   chrombpnet contribs_bw            contribution scores of both lines' models on the featured peaks
#   python peaks.py                   every peak's sequence
#   python ablation.py                motif sites in peaks, each with the site scrambled and a scrambled control
#   python predict.py                 both models' predictions for every peak and every motif-removal pair
#   python prepare_bundle.py build    everything into data/ and chrombpnet_data.tar.gz
#
# data/ (what the notebook reads):
#   regions.npz         750 featured peaks: 2,114 bp sequence; observed coverage, observed Tn5 cuts and both models'
#                       contribution scores over the central 1,000 bp
#   peaks.tsv.gz        every consensus peak off the amplicons with equal copy number: Step 3's call and fold change, observed reads, both
#                       models' predicted reads
#   ablation.tsv.gz     motif sites in peaks: predicted change when the site (or, as a control, a random 6 bp) is
#                       scrambled, in each model
#   motifs.json         each model's TF-MoDISco motifs: logo matrix, number of seqlets, best database matches
#   qc.json             each model's evaluation: test-chromosome accuracy, bias checks, training curve
#   backgrounds.txt     16 genomic background sequences (GC-matched non-peak regions) for planting motifs
#   diagram.png         chromBPNet's architecture
#   models/             chromBPNet_<line>.h5, the bias-free model of each line
import glob, json, os, re, shutil, sys, tarfile
import numpy as np, pandas as pd, pyBigWig, pyfaidx, h5py

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))     # the workshop folder
RES = f"{ROOT}/04_sequence_models/results/chrombpnet"
OUT = f"{ROOT}/04_sequence_models/data"
WORK = f"{ROOT}/04_sequence_models/prep/work"
LINES = ["DM", "HSR"]
MODEL = {l: f"{RES}/models/COLO320{l}_fold0" for l in LINES}               # one model per line: replicates pooled, fold 0
COVERAGE = {l: f"{ROOT}/02_pipeline/results/bwa/merged_replicate/bigwig/COLO320{l}.mRp.clN.bigWig" for l in LINES}
GENOME = f"{ROOT}/02_pipeline/results/genome/genome.fa"
CALLS = f"{WORK}/calls.tsv"                   # Step 3's DESeq2 calls, from calls.R
N_PER_CLASS = 250
HALF_IN, HALF_TRACK = 1057, 500               # 2,114 bp model input; 1,000 bp of tracks
TEST_CHROMS = ["chr1", "chr3", "chr6"]        # fold 0: the models never saw these


def all_peaks():
    calls = pd.read_csv(CALLS, sep="\t")
    info = pd.read_csv(f"{ROOT}/03_differential/data/peaks.tsv.gz", sep="\t",
                       usecols=["peak", "gene", "gc", "COLO320DM_REP1", "COLO320DM_REP2", "COLO320HSR_REP1", "COLO320HSR_REP2"])
    p = calls.merge(info, on="peak")
    p = p[(p.region_class == "off-amplicon, CN equal") & p.chr.isin([f"chr{c}" for c in list(range(1, 23)) + ["X"]])]
    p["center"] = (p.start + p.end) // 2
    p["reads_DM"] = p.COLO320DM_REP1 + p.COLO320DM_REP2
    p["reads_HSR"] = p.COLO320HSR_REP1 + p.COLO320HSR_REP2
    return p.reset_index(drop=True)


def summit(chrom, start, end):
    # consensus peaks have no summit: use the base with the most coverage (both lines, 10 bp average) inside the peak
    cov = sum(track(COVERAGE[l], chrom, start, end) for l in LINES)
    return start + int(np.argmax(np.convolve(cov, np.ones(10) / 10, mode="same")))


def select():
    shutil.rmtree(OUT, ignore_errors=True); os.makedirs(OUT, exist_ok=True); os.makedirs(WORK, exist_ok=True)
    fa = pyfaidx.Fasta(GENOME)
    p = all_peaks()
    p["seq"] = [str(fa[c][x - HALF_IN:x + HALF_IN]).upper() for c, x in zip(p.chr, p.center)]
    # featured peaks: the most significant up in DM and up in HSR, and strong unchanged ones, centered on their summit
    up_dm = p[p.call == "up in DM"].nsmallest(N_PER_CLASS, "padj")
    up_hsr = p[p.call == "up in HSR"].nsmallest(N_PER_CLASS, "padj")
    same = p[(p.call == "not significant") & (p.reads_DM + p.reads_HSR >= (p.reads_DM + p.reads_HSR).quantile(0.75))]
    sel = pd.concat([up_dm, up_hsr, same.sample(N_PER_CLASS, random_state=0)]).reset_index(drop=True)
    sel["center"] = [summit(c, s, e) for c, s, e in zip(sel.chr, sel.start, sel.end)]
    sel.drop(columns="seq").to_csv(f"{WORK}/featured.tsv", sep="\t", index=False)
    pd.DataFrame({"chr": sel.chr, "s": sel.center - 500, "e": sel.center + 500, "name": sel.peak, "score": 0,
                  "strand": ".", "a": 0, "b": 0, "c": 0, "summit": 500}).to_csv(
        f"{WORK}/featured.narrowPeak", sep="\t", header=False, index=False)       # chromBPNet's region format
    print("featured peaks:", sel.call.value_counts().to_dict())
    # genomic background sequences for planting motifs
    bg = pd.read_csv(f"{MODEL['DM']}/auxiliary/filtered.nonpeaks.bed", sep="\t", header=None)
    out = []
    for c, s in zip(*bg.sample(200, random_state=0)[[0, 1]].values.T):
        x = str(fa[c][s:s + 2 * HALF_IN]).upper()
        if len(x) == 2 * HALF_IN and set(x) <= set("ACGT"):
            out.append(x)
    open(f"{OUT}/backgrounds.txt", "w").write("\n".join(out[:16]) + "\n")
    with open(f"{WORK}/chrombpnet_models.txt", "w") as f:
        for l in LINES:
            f.write(f"{l} {MODEL[l]}/models/chrombpnet_nobias.h5\n")


def track(path, chrom, start, end):
    bw = pyBigWig.open(path)
    v = np.nan_to_num(np.array(bw.values(chrom, start, end), dtype=np.float32))
    bw.close()
    return v


def build_regions():
    sel = pd.read_csv(f"{WORK}/featured.tsv", sep="\t")
    fa = pyfaidx.Fasta(GENOME)
    names = [f"{kind}_{l}" for l in LINES for kind in ("coverage", "cuts", "contrib")]
    out = np.zeros((len(sel), len(names), 2 * HALF_TRACK), dtype=np.float16)
    for i, r in sel.iterrows():
        s, e = r.center - HALF_TRACK, r.center + HALF_TRACK
        for k, l in enumerate(LINES):
            out[i, 3 * k] = track(COVERAGE[l], r.chr, s, e)                                     # reads covering each base
            out[i, 3 * k + 1] = track(f"{MODEL[l]}/auxiliary/data_unstranded.bw", r.chr, s, e)  # Tn5 cuts at each base
            out[i, 3 * k + 2] = track(f"{WORK}/contribs/chrombpnet_{l}.counts_scores.bw", r.chr, s, e)
    seqs = [str(fa[c][x - HALF_IN:x + HALF_IN]).upper() for c, x in zip(sel.chr, sel.center)]
    np.savez_compressed(f"{OUT}/regions.npz", tracks=out, track_names=np.array(names), seq=np.array(seqs),
                        peak=sel.peak.values, chrom=sel.chr.values, center=sel.center.values, gene=sel.gene.fillna("").values,
                        call=sel.call.values, log2FC=sel.log2FC.values, padj=sel.padj.values, gc=sel.gc.values)
    print("regions:", out.shape)


def build_tables():
    order = pd.read_csv(f"{WORK}/peaks.tsv", sep="\t", usecols=["peak"])       # the rows predict.py predicted
    for l in LINES:
        order[f"pred_{l}"] = np.load(f"{WORK}/pred_chrombpnet_COLO320{l}_peak.npy")   # predicted log reads (natural log)
    p = all_peaks().merge(order, on="peak")                                          # off the amplicons
    p["chromosomes"] = np.where(p.chr.isin(TEST_CHROMS), "test", "train")
    p[["peak", "chr", "center", "gene", "gc", "call", "log2FC", "padj", "reads_DM", "reads_HSR", "pred_DM", "pred_HSR",
       "chromosomes"]].to_csv(f"{OUT}/peaks.tsv.gz", sep="\t", index=False)
    a = pd.read_csv(f"{WORK}/ablate.tsv", sep="\t")
    a["peak"] = order.peak.values[a.peak_idx]
    for l in LINES:
        x = np.load(f"{WORK}/pred_chrombpnet_COLO320{l}_ablate.npy")
        a[f"change_{l}"] = (x[1::2] - x[0::2]) / np.log(2)                # scrambled minus original, log2
    a = a[a.peak.isin(set(p.peak))][["peak", "motif", "kind", "start", "change_DM", "change_HSR"]]
    a.to_csv(f"{OUT}/ablation.tsv.gz", sep="\t", index=False)
    print("peaks:", len(p), " ablation pairs:", len(a))


def report_matches(html):
    out = {}
    for row in re.findall(r"<tr>(.*?)</tr>", open(html).read(), flags=re.S):
        cells = [re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", row, flags=re.S)]
        if cells and "pattern" in cells[0]:
            out[cells[0]] = [(cells[i], float(cells[i + 1])) for i in (4, 7, 10)
                             if i + 1 < len(cells) and cells[i] not in ("", "nan", "NaN")]
    return out


def trim(cwm, frac=0.3):
    score = np.abs(cwm).sum(axis=1)
    keep = np.where(score >= frac * score.max())[0]
    return cwm[max(keep.min() - 2, 0):keep.max() + 3]


def build_motifs():
    out = {}
    for l in LINES:
        d = f"{RES}/contribs/COLO320{l}_fold0"
        matches = report_matches(f"{d}/modisco_counts_report/motifs.html")
        pats = []
        with h5py.File(f"{d}/modisco_counts.h5") as f:
            for sign in ("pos_patterns", "neg_patterns"):
                for name in sorted(f.get(sign, {}), key=lambda x: int(x.split("_")[1])):
                    g = f[sign][name]
                    pats.append(dict(pattern=f"{sign}.{name}", n_seqlets=int(g["seqlets"]["n_seqlets"][()][0]),
                                     cwm=trim(g["contrib_scores"][()]).round(4).tolist(),
                                     matches=matches.get(f"{sign}.{name}", [])))
        out[l] = pats
    json.dump(out, open(f"{OUT}/motifs.json", "w"))
    print("motifs:", {k: len(v) for k, v in out.items()})


def build_qc():
    qc = {}
    for l in LINES:
        e = f"{MODEL[l]}/evaluation"
        qc[l] = dict(test=json.load(open(f"{e}/chrombpnet_metrics.json")), bias_only=json.load(open(f"{e}/bias_metrics.json")),
                     bias_response=float(open(f"{e}/chrombpnet_nobias_max_bias_response.txt").read().split("_")[1]),  # "corrected_<max>_<per Tn5 motif>"
                     training=pd.read_csv(f"{MODEL[l]}/logs/chrombpnet.log").to_dict(orient="list"),
                     profile_example=profile_example(l))
    json.dump(qc, open(f"{OUT}/qc.json", "w"))


def profile_example(l):
    # the full model's (bias included) prediction, from chromBPNet's evaluation, at a typical strong test-chromosome
    # peak (inside a consensus peak with equal copy number), next to those cuts. Typical: among the top 10% of peaks by
    # observed cuts, the one where the model's predicted total is as close as usual (median predicted/observed ratio);
    # picking by observed cuts alone would favor peaks the model underpredicts.
    cons = all_peaks()
    with h5py.File(f"{MODEL[l]}/evaluation/chrombpnet_predictions.h5") as f:
        chrom = f["coords/coords_chrom"][()].astype(str); center = f["coords/coords_center"][()]
        is_peak = f["coords/coords_peak"][()] == 1
        for c in set(chrom):
            q = cons[cons.chr == c].sort_values("start")
            k = np.searchsorted(q.start.values, center, side="right") - 1
            inside = (chrom == c) & (k >= 0) & (center < np.r_[q.end.values, [0]][np.clip(k, 0, None)])
            is_peak &= (chrom != c) | inside
        bw = pyBigWig.open(f"{MODEL[l]}/auxiliary/data_unstranded.bw")
        obs = np.array([bw.stats(c, x - HALF_TRACK, x + HALF_TRACK, type="sum")[0] or 0 if p else -1
                        for c, x, p in zip(chrom, center, is_peak)])
        bw.close()
        top = np.where(obs >= np.quantile(obs[obs >= 0], 0.9))[0]
        ratio = np.exp(f["predictions/logcounts"][()][top]) / obs[top]
        k = int(top[np.argsort(ratio)[len(top) // 2]])
        prof, logcounts = f["predictions/profs"][k], float(f["predictions/logcounts"][k])
    # the bias-free model at the same peak: what the notebook's live predictions show
    import tensorflow as tf
    seq = str(pyfaidx.Fasta(GENOME)[chrom[k]][center[k] - HALF_IN:center[k] + HALF_IN]).upper()
    x = np.zeros((1, 2 * HALF_IN, 4), np.float32)
    for j, base in enumerate(seq):
        if base in "ACGT":
            x[0, j, "ACGT".index(base)] = 1
    logits, lc = tf.keras.models.load_model(f"{MODEL[l]}/models/chrombpnet_nobias.h5", compile=False).predict(x, verbose=0)
    nb = np.exp(logits[0] - logits[0].max())
    return dict(chrom=chrom[k], center=int(center[k]), profile_sum=float(prof.sum()),
                observed=track(f"{MODEL[l]}/auxiliary/data_unstranded.bw", chrom[k], center[k] - HALF_TRACK,
                               center[k] + HALF_TRACK).round(2).tolist(),
                predicted=(prof / prof.sum() * np.exp(logcounts)).round(3).tolist(),
                predicted_bias_free=(nb / nb.sum() * np.exp(lc[0, 0])).round(3).tolist())


def build_diagram():
    # chromBPNet's architecture (chrombpnet/training/models/chrombpnet_with_bias_model.py), with our models' sizes:
    # main model 512 filters and 8 dilated layers; ENCODE's bias model 128 filters and 4 dilated layers
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    fig = plt.figure(figsize=(16, 8.0))
    ax = fig.add_axes([0, 0.35, 1, 0.65]); ax.axis("off"); ax.set_xlim(0, 16); ax.set_ylim(0, 4.95)
    RED, BLUE, GRAY, YELLOW = "#f4c7c3", "#c9daf8", "#f0f0f0", "#fff2cc"
    R = fig.canvas.get_renderer()
    M, GAP = 0.18, 0.12                                   # margin inside each box; space between title and text

    def size(t):                                          # a text's width and height in data units
        bb = t.get_window_extent(R).transformed(ax.transData.inverted())
        return bb.width, bb.height

    def box(x, yc, title, body, color, fs=12, min_w=0):
        # a box just big enough for its text, left edge at x, centered on yc; returns its edges
        t1 = ax.text(0, 0, title, ha="center", va="bottom", fontsize=fs + 1, weight="bold")
        t2 = ax.text(0, 0, body, ha="center", va="top", fontsize=fs, linespacing=1.35)
        (w1, h1), (w2, h2) = size(t1), size(t2)
        w, h = max(w1 + 2 * M, w2 + 2 * M, min_w), h1 + GAP + h2 + 2 * M
        top = yc + h / 2
        t1.set_position((x + w / 2, top - M - h1)); t2.set_position((x + w / 2, top - M - h1 - GAP))
        ax.add_patch(FancyBboxPatch((x, yc - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.1", fc=color, ec="gray"))
        return dict(l=x, r=x + w, b=yc - h / 2, t=top, y=yc)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="->", color="gray", lw=1.6, shrinkA=0, shrinkB=0))

    TOP, BOTTOM = 3.7, 0.85                               # centers of the main-model row and the bias-model row
    dna = box(0.1, (TOP + BOTTOM) / 2, "DNA sequence", "2,114 bp, written as\n2,114 rows of A/C/G/T\n(one-hot)", GRAY)
    l1 = box(dna["r"] + 0.6, TOP, "Layer 1", "512 filters, each\n21 bp wide, slide\nalong the DNA:\nmotif scanners", RED)
    l29 = box(l1["r"] + 0.5, TOP, "Layers 2–9", "each combines the layer\nbelow at positions 2, 4,\n8 … 256 bp apart:\nmotif combinations\nand spacing", RED)
    hx = l29["r"] + 0.6
    prof = box(hx, TOP + 0.6, "Profile head", "a score for each\nof the 1,000 bp", RED, fs=11, min_w=2.3)
    cnt = box(hx, TOP - 0.6, "Counts head", "average over the\nwindow: log reads", RED, fs=11, min_w=2.3)
    ax.text(l1["l"], l29["t"] + 0.12, "Main model (trained on this line's reads)", fontsize=14, weight="bold", color="firebrick")
    bias = box(l1["l"], BOTTOM, "Same design, smaller", "128 filters, 4 dilated layers: each output base sees only ~150 bp of DNA,\n"
               "enough for Tn5's own sequence preference but not for TF motif combinations.\n"
               "Outputs its own per-base scores and log reads.", BLUE)
    ax.text(bias["l"], bias["t"] + 0.12, "Tn5 bias model (pretrained by ENCODE, frozen)", fontsize=14, weight="bold", color="steelblue")
    comb = box(max(prof["r"], cnt["r"]) + 0.6, TOP - 0.55, "Combined",
               "profile: add the two models'\nscores at each base, then\nturn them into shares of\nthe cuts (softmax)\n\n"
               "counts: predicted reads\n= main model's reads\n+ bias model's reads", YELLOW)
    arrow(dna["r"], dna["y"] + 0.3, l1["l"], l1["y"]); arrow(dna["r"], dna["y"] - 0.3, bias["l"], bias["y"])
    arrow(l1["r"], TOP, l29["l"], TOP)
    arrow(l29["r"], TOP, prof["l"], prof["y"]); arrow(l29["r"], TOP, cnt["l"], cnt["y"])
    arrow(prof["r"], prof["y"], comb["l"], min(prof["y"], comb["t"] - 0.3))
    arrow(cnt["r"], cnt["y"], comb["l"], cnt["y"])
    arrow(bias["r"], bias["y"], comb["l"], comb["b"] + 0.3)

    # dilation sketch: which input positions one output position sees, layer by layer
    bx = fig.add_axes([0.04, 0.02, 0.58, 0.25]); bx.axis("off")
    n = 64; rows = [("layer 1", 0)] + [(f"layer {k + 2}", 2 ** (k + 1)) for k in range(4)]
    # draw top-down: the output position at the top row, the positions it depends on below
    centers = [{32}]
    for r in range(len(rows) - 1, 0, -1):
        d = rows[r][1]
        centers.append({s + o for s in centers[-1] for o in (-d, 0, d)})
    for row, (label, d) in enumerate(reversed(rows)):
        y = row
        on = centers[row]
        bx.scatter(range(n), [y] * n, s=10, color="lightgray")
        bx.scatter(sorted(p for p in on if 0 <= p < n), [y] * len([p for p in on if 0 <= p < n]), s=28, color="firebrick")
        bx.text(-2, y, f"{label}: gap {d}" if d else label, ha="right", va="center", fontsize=11)
        if row < len(rows) - 1:
            dd = rows[len(rows) - 1 - row][1]
            for s in on:
                for o in (-dd, 0, dd):
                    if 0 <= s + o < n:
                        bx.plot([s, s + o], [y, y + 1], color="firebrick", lw=0.6, alpha=0.5)
    bx.set_ylim(len(rows) - 0.5, -0.8); bx.set_xlim(-14, n)
    bx.set_title("Dilated layers 2–5 of 2–9: each looks at points twice as far apart as the layer before",
                 fontsize=12, loc="left")
    fig.text(0.66, 0.24, "Training: compare the combined prediction\nwith the observed cuts and reads, and\n"
             "adjust only the main model.\n\nAfter training: drop the bias model.\nThe main model alone is the bias-free\n"
             "model used in this notebook.", fontsize=13, va="top", linespacing=1.4)
    fig.savefig(f"{OUT}/diagram.png", dpi=110, bbox_inches="tight")


def build_models():
    os.makedirs(f"{OUT}/models", exist_ok=True)
    for l in LINES:
        shutil.copy(f"{MODEL[l]}/models/chrombpnet_nobias.h5", f"{OUT}/models/chrombpnet_{l}.h5")


if __name__ == "__main__":
    if sys.argv[1] == "select":
        select()
    else:
        build_regions(); build_tables(); build_motifs(); build_qc(); build_diagram(); build_models()
        with tarfile.open(f"{ROOT}/04_sequence_models/chrombpnet_data.tar.gz", "w:gz") as t:
            t.add(OUT, arcname="data", filter=lambda x: None if os.path.basename(x.name).startswith(".") else x)
        print("archive MB:", round(os.path.getsize(f"{ROOT}/04_sequence_models/chrombpnet_data.tar.gz") / 1e6, 1))
