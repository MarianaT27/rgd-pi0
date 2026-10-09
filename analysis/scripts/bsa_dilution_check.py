import numpy as np

TARGETS = ["LD2", "CxC", "Cu", "Sn"]

def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}

d = {t: load(f"results_full/bsa_{t}.csv") for t in TARGETS}

# bins present (good fit) in all four targets
good = []
for t in TARGETS:
    ok = np.isfinite(d[t]["A_LU"]) & (d[t]["sA_LU"] > 0) & (np.abs(d[t]["A_LU"]) < 0.5)
    good.append(set(d[t]["bin4d"][ok].astype(int)))
common = set.intersection(*good)
print(f"bins common to all targets: {len(common)}")

def wmean(v, s):
    w = 1.0 / s**2
    return np.sum(w * v) / np.sum(w), 1.0 / np.sqrt(np.sum(w))

print("\n i_z  target   A_raw             S/(S+B)   A_LU              A_raw/A_raw(D)   A_LU/A_LU(D)")
for iz in range(5):
    ref = {}
    for t in TARGETS:
        b = d[t]
        sel = np.array([(int(k) in common) and (int(z) == iz) for k, z in zip(b["bin4d"], b["i_z"])])
        ar, ear = wmean(b["A_raw"][sel], b["sA_raw"][sel])
        al, eal = wmean(b["A_LU"][sel], b["sA_LU"][sel])
        w = b["n_events"][sel]
        sf = np.sum(w * b["signal_fraction"][sel]) / np.sum(w)
        if t == "LD2":
            ref = {"ar": ar, "al": al}
        print(f" {iz:3d}  {t:5s}   {ar:+.4f} +- {ear:.4f}   {sf:.3f}     {al:+.4f} +- {eal:.4f}   "
              f"{ar / ref['ar']:6.2f}           {al / ref['al']:6.2f}")
    print()
