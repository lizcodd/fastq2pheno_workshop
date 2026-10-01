# Predicted log counts of both lines' chromBPNet models for work/{peak,ablate}_seqs.txt
# -> work/pred_chrombpnet_<model>_<set>.npy. Runs in the chromBPNet container (TF 2.8); minutes on a GPU, hours on a CPU.
import os, numpy as np, tensorflow as tf
E = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results", "chrombpnet", "models"))
MODELS = ["COLO320DM", "COLO320HSR"]
LUT = np.zeros(256, dtype=np.int64) + 4; LUT[[ord(c) for c in "ACGT"]] = np.arange(4)
EYE = np.vstack([np.eye(4, dtype=np.float32), np.zeros((1, 4), np.float32)])

def batches(path, n=2048):
    buf = []
    for line in open(path):
        buf.append(line.strip())
        if len(buf) == n:
            yield buf; buf = []
    if buf:
        yield buf

def onehot(seqs):
    idx = LUT[np.frombuffer("".join(seqs).encode(), dtype=np.uint8)].reshape(len(seqs), -1)
    return EYE[idx]

for name in MODELS:
    m = tf.keras.models.load_model(f"{E}/{name}_fold0/models/chrombpnet_nobias.h5", compile=False)
    for s in ["peak", "ablate"]:
        out = [m.predict(onehot(b), batch_size=256, verbose=0)[1][:, 0] for b in batches(f"work/{s}_seqs.txt")]
        np.save(f"work/pred_chrombpnet_{name}_{s}.npy", np.concatenate(out).astype(np.float32))
        print(name, s, "done", flush=True)
