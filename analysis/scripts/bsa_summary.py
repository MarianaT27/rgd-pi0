import sys
import numpy as np

def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}

def wmean(a, s):
    w = 1.0 / s**2
    m = np.sum(w * a) / np.sum(w)
    e = 1.0 / np.sqrt(np.sum(w))
    chi2 = np.sum(w * (a - m)**2) / max(len(a) - 1, 1)
    return m, e, chi2

for t in sys.argv[1:]:
    d = load(f"bsa3_{t}.csv")
    ok = np.isfinite(d["A_LU"]) & np.isfinite(d["sA_LU"]) & (d["sA_LU"] > 0) & (np.abs(d["A_LU"]) < 0.5)
    a, s = d["A_LU"][ok], d["sA_LU"][ok]
    m, e, c = wmean(a, s)
    print(f"\n{t}: {ok.sum()} bins used of {len(ok)}  |  A_LU = {m:+.4f} +- {e:.4f}  (chi2/ndf = {c:.2f})")
    print(f"  bins > 3 sigma from zero: {np.sum(np.abs(a / s) > 3)}")
    print("  i_z   n_bins   A_LU      error")
    for iz in np.unique(d["i_z"][ok]).astype(int):
        k = d["i_z"][ok] == iz
        mz, ez, _ = wmean(a[k], s[k])
        print(f"  {iz:3d}   {k.sum():6d}   {mz:+.4f}   {ez:.4f}")