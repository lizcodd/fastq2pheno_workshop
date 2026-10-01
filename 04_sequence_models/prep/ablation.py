# In-peak motif ablation: for each consensus peak with a motif site in its central 300 bp, the 2,114 bp sequence with
# that site scrambled, plus the same peak with a random 6-bp stretch near the center scrambled (control).
# Writes work/ablate.tsv and work/ablate_seqs.txt (original, then mutated, per row). Run after peaks.py.
import re, numpy as np, pandas as pd
rng = np.random.default_rng(0)
MOTIFS = {"E-box, MYC/MAX-type flanks (CCACGTGG)": r"[CG]CACGTG[CG]", "E-box, USF/TFE-type flanks (TCACGTGA)": r"TCACGTGA",
          "E-box, any (CACGTG)": r"CACGTG", "AP-1 (TGASTCA)": r"TGA[CG]TCA", "SIX1-like (GTAATATGA)": r"GTAAT[AC]TGA|TCA[GT]ATTAC",
          "CTCF core (CCASYAGRKGG)": r"CCA[CG][CT]AG[AG][GT]GG|CC[AC][CT]CT[AG][GC]TGG"}
peaks = pd.read_csv("work/peaks.tsv", sep="\t")
seqs = [l.strip() for l in open("work/peak_seqs.txt")]
rows, out = [], []
def scramble(s):
    # a shuffle of the same letters; a single repeated letter (e.g. AAAAAA) can't be shuffled, so use random letters
    if len(set(s)) == 1:
        return "".join(rng.choice([b for b in "ACGT" if b != s[0]], len(s)))
    while True:
        t = "".join(rng.permutation(list(s)))
        if t != s:
            return t
for i, s in enumerate(seqs):
    center = s[907:1207]
    for name, pat in MOTIFS.items():
        m = re.search(pat, center)
        if not m:
            continue
        a, b = 907 + m.start(), 907 + m.end()
        rows.append(dict(peak_idx=i, motif=name, kind="motif", start=a)); out += [s, s[:a] + scramble(s[a:b]) + s[b:]]
        c = 907 + int(rng.integers(0, 294))
        rows.append(dict(peak_idx=i, motif=name, kind="control", start=c)); out += [s, s[:c] + scramble(s[c:c + 6]) + s[c + 6:]]
pd.DataFrame(rows).to_csv("work/ablate.tsv", sep="\t", index=False)
open("work/ablate_seqs.txt", "w").write("\n".join(out) + "\n")
print(pd.DataFrame(rows).query("kind == 'motif'").motif.value_counts())
