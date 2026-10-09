"""
bsa.py : beam-spin asymmetry A_LU (sin phi_h) for one sample (diag3 or full; f10 has no BSA csv).

Usage (from the repo root):
    python3 analysis/scripts/bsa.py full

BEFORE running this script, these must already exist for the chosen sample:
  1. Stage A slims, make_grid (config/binning/grid_*.json) and Stage B
     (binned ROOT files, e.g. analysis/full/binned_full/<T>.root)
  2. pi0.bsa for LD2, CxC, Cu and Sn, which writes the BSA csv files:
       PYTHONPATH=python python3 -m pi0.bsa --file <binned>/<T>.root --config config \
           --out <results>/bsa_<T>.csv --allow-unpublishable
     (A_LU is already divided by the beam polarization in config/cuts.json and by S/(S+B))
  3. pi0.ratio for CxC (ratio_CxC.csv), only used to place each z bin at its <z> on the x axis
     (run_full_chain.sh does steps 1 to 3 for the full sample)

Reads:   full : analysis/full/results_full/bsa_{LD2,CxC,Cu,Sn}.csv
         diag3: analysis/diag3/bsa3_{LD2,CxC,Cu,Sn}.csv
         and <results>/ratio_CxC.csv of the SAME sample for <z>
Writes:  analysis/<sample>/bsa_vs_z_<sample>.pdf/.png   A_LU vs <z> (top) and A_LU^A / A_LU^D (bottom)
Prints:  1. summary: A_LU averaged over all bins and per z bin, for each target
         2. dilution check: A_raw, S/(S+B) and A_LU per z bin, on the bins common to all 4 targets

Bins with |A_LU| >= 0.5 or no valid error are treated as failed fits and left out.
Averages are inverse-variance weighted.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLES = {
    "diag3": {"dir": "analysis/diag3", "results": "results3",     "bsa": "analysis/diag3/bsa3_{t}.csv",
              "label": "3% sample"},
    "f10":   {"dir": "analysis/f10",   "results": "results_f10",  "bsa": None,
              "label": "10% sample"},
    "full":  {"dir": "analysis/full",  "results": "results_full", "bsa": "analysis/full/results_full/bsa_{t}.csv",
              "label": "full sample"},
}
TARGETS = ["LD2", "CxC", "Cu", "Sn"]
LABELS = {"LD2": "D", "CxC": "C", "Cu": "Cu", "Sn": "Sn"}
COLORS = {"LD2": "black", "CxC": "tab:blue", "Cu": "tab:orange", "Sn": "tab:red"}

if len(sys.argv) != 2 or sys.argv[1] not in SAMPLES:
    sys.exit("usage: python3 analysis/scripts/bsa.py <diag3|full>")
S = sys.argv[1]
if SAMPLES[S]["bsa"] is None:
    sys.exit(f"no BSA csv files exist for sample {S} (pi0.bsa was not run on it)")
OUT = SAMPLES[S]["dir"]
RES = f"{OUT}/{SAMPLES[S]['results']}"


def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


def wmean(v, s):
    w = 1.0 / s**2
    m = np.sum(w * v) / np.sum(w)
    chi2 = np.sum(w * (v - m)**2) / max(len(v) - 1, 1)
    return m, 1.0 / np.sqrt(np.sum(w)), chi2


d = {t: load(SAMPLES[S]["bsa"].format(t=t)) for t in TARGETS}
good = {t: np.isfinite(d[t]["A_LU"]) & np.isfinite(d[t]["sA_LU"]) & (d[t]["sA_LU"] > 0) & (np.abs(d[t]["A_LU"]) < 0.5)
        for t in TARGETS}
zbins = sorted(set(d["LD2"]["i_z"].astype(int)))

# <z> of each z bin from this sample's ratio_CxC.csv (Y_A-weighted)
r = load(f"{RES}/ratio_CxC.csv")
rk = np.isfinite(r["Y_A"]) & (r["Y_A"] > 0)
zc = np.array([np.sum(r["Y_A"][rk & (r["i_z"] == iz)] * r["z_mean"][rk & (r["i_z"] == iz)])
               / np.sum(r["Y_A"][rk & (r["i_z"] == iz)]) for iz in zbins])

# ---------------------------------------------------------------- 1. summary
print(f"=== A_LU summary, sample {S}")
res = {}
for t in TARGETS:
    ok = good[t]
    a, s, iz = d[t]["A_LU"][ok], d[t]["sA_LU"][ok], d[t]["i_z"][ok].astype(int)
    m, e, c = wmean(a, s)
    print(f"\n{LABELS[t]} ({t}): {ok.sum()} bins used of {len(ok)}  |  A_LU = {m:+.4f} +- {e:.4f}  (chi2/ndf = {c:.2f})")
    print(f"  bins > 3 sigma from zero: {np.sum(np.abs(a / s) > 3)}")
    print("  i_z    <z>   n_bins   A_LU      error")
    ms, es = [], []
    for k, z in zip(zbins, zc):
        mz, ez, _ = wmean(a[iz == k], s[iz == k])
        ms.append(mz); es.append(ez)
        print(f"  {k:3d}  {z:.3f}  {np.sum(iz == k):6d}   {mz:+.4f}   {ez:.4f}")
    res[t] = (np.array(ms), np.array(es))

# ---------------------------------------------------------------- plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 7), sharex=True,
                               gridspec_kw={"height_ratios": [3, 2], "hspace": 0.05})
for i, t in enumerate(TARGETS):
    m, e = res[t]
    ax1.errorbar(zc + (i - 1.5) * 0.006, m, yerr=e, fmt="o", color=COLORS[t], label=LABELS[t], capsize=2)
ax1.axhline(0, color="gray", lw=0.8, ls="--")
ax1.set_ylabel(r"$A_{LU}^{\sin\phi_h}$")
ax1.legend(frameon=False, ncol=4)
ax1.text(0.02, 0.95, f"DIAGNOSTIC ({SAMPLES[S]['label']}, RG-A photon ID)", transform=ax1.transAxes,
         va="top", fontsize=8, color="gray")
mD, eD = res["LD2"]
for i, t in enumerate(["CxC", "Cu", "Sn"]):
    m, e = res[t]
    q = m / mD
    eq = np.abs(q) * np.sqrt((e / m)**2 + (eD / mD)**2)
    ax2.errorbar(zc + (i - 1) * 0.006, q, yerr=eq, fmt="o", color=COLORS[t], capsize=2)
ax2.axhline(1, color="gray", lw=0.8, ls="--")
ax2.set_ylabel(r"$A_{LU}^{A} / A_{LU}^{D}$")
ax2.set_xlabel(r"$\langle z \rangle$")
ax2.set_ylim(-1.5, 3.5)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/bsa_vs_z_{S}.{ext}", dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- 2. dilution check
common = set.intersection(*[set(d[t]["bin4d"][good[t]].astype(int)) for t in TARGETS])
print(f"\n=== dilution check, {len(common)} 4D bins common to all targets")
print(" i_z  target   A_raw             S/(S+B)   A_LU              A_raw/A_raw(D)   A_LU/A_LU(D)")
for k in zbins:
    ref = {}
    for t in TARGETS:
        b = d[t]
        sel = np.array([(int(n) in common) and (int(z) == k) for n, z in zip(b["bin4d"], b["i_z"])])
        ar, ear, _ = wmean(b["A_raw"][sel], b["sA_raw"][sel])
        al, eal, _ = wmean(b["A_LU"][sel], b["sA_LU"][sel])
        w = b["n_events"][sel]
        sf = np.sum(w * b["signal_fraction"][sel]) / np.sum(w)
        if t == "LD2":
            ref = {"ar": ar, "al": al}
        print(f" {k:3d}  {t:5s}   {ar:+.4f} +- {ear:.4f}   {sf:.3f}     {al:+.4f} +- {eal:.4f}   "
              f"{ar / ref['ar']:6.2f}           {al / ref['al']:6.2f}")
    print()

print(f"wrote {OUT}/bsa_vs_z_{S}.pdf/.png")
