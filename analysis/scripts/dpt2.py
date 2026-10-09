"""
dpt2.py : pT broadening Delta<pT2> for one sample (diag3, f10 or full).

Usage (from the repo root):
    python3 analysis/scripts/dpt2.py full

BEFORE running this script, these must already exist for the chosen sample:
  1. Stage A slims, make_grid (config/binning/grid_*.json) and Stage B
     (binned ROOT files, e.g. analysis/full/binned_full/<T>.root)
  2. pi0.broadening for CxC, Cu and Sn, which writes the input csv files:
       PYTHONPATH=python python3 -m pi0.broadening --target <binned>/<T>.root --ld2 <binned>/LD2.root \
           --config config --out <results>/broadening_<T>.csv --allow-unpublishable
     (run_full_chain.sh does steps 1 and 2 for the full sample)

Reads:   analysis/<sample>/<results>/broadening_{CxC,Cu,Sn}.csv
Writes:  analysis/<sample>/dpt2_vs_A13_<sample>.pdf/.png   Delta<pT2> vs A^(1/3), one colour per z bin,
                                                            with a straight line through the origin
Prints:  Delta<pT2> per target and z bin, and the fitted slope.

Only 3D bins that are good (finite delta, positive errors) in all three targets are used,
so every target is averaged over the same bins. Weights are 1/spt2_D^2 (the LD2 error),
which are the same for all targets in a given bin.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLES = {
    "diag3": {"dir": "analysis/diag3", "results": "results3",     "label": "3% sample"},
    "f10":   {"dir": "analysis/f10",   "results": "results_f10",  "label": "10% sample"},
    "full":  {"dir": "analysis/full",  "results": "results_full", "label": "full sample"},
}
TARGETS = {"CxC": 12.011, "Cu": 63.546, "Sn": 118.71}
LABELS = {"CxC": "C", "Cu": "Cu", "Sn": "Sn"}

if len(sys.argv) != 2 or sys.argv[1] not in SAMPLES:
    sys.exit("usage: python3 analysis/scripts/dpt2.py <diag3|f10|full>")
S = sys.argv[1]
OUT = SAMPLES[S]["dir"]
RES = f"{OUT}/{SAMPLES[S]['results']}"


def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


d = {t: load(f"{RES}/broadening_{t}.csv") for t in TARGETS}

good = {}
for t in TARGETS:
    ok = np.isfinite(d[t]["delta"]) & (d[t]["sdelta"] > 0) & (d[t]["spt2_D"] > 0)
    good[t] = set(d[t]["bin3d"][ok].astype(int))
common = set.intersection(*good.values())

zbins = sorted(set(d["CxC"]["i_z"].astype(int)))
x = {t: A ** (1 / 3) for t, A in TARGETS.items()}
xs = np.array([x[t] for t in TARGETS])

print(f"=== Delta<pT2> [GeV2], sample {S} ({RES}), {len(common)} common 3D bins")
print("i_z  n_bins   " + "   ".join(f"{LABELS[t]:>16s}" for t in TARGETS) + "   slope [GeV2]")
fig, ax = plt.subplots(figsize=(6.5, 4.8))
colors = plt.cm.viridis(np.linspace(0, 0.9, len(zbins)))
for iz, col in zip(zbins, colors):
    vals, errs = [], []
    for t in TARGETS:
        b = d[t]
        sel = np.array([(int(k) in common) and (int(z) == iz) for k, z in zip(b["bin3d"], b["i_z"])])
        nb = sel.sum()
        w = 1.0 / b["spt2_D"][sel] ** 2
        vals.append(np.sum(w * b["delta"][sel]) / np.sum(w))
        errs.append(np.sqrt(np.sum(w ** 2 * b["sdelta"][sel] ** 2)) / np.sum(w))
    vals, errs = np.array(vals), np.array(errs)
    slope = np.sum(xs * vals / errs ** 2) / np.sum(xs ** 2 / errs ** 2)
    print(f"{iz:3d}  {nb:6d}   " + "   ".join(f"{v:.5f} +- {e:.5f}" for v, e in zip(vals, errs)) + f"   {slope:.5f}")
    ax.errorbar(xs, vals, yerr=errs, fmt="o", color=col, capsize=2, label=f"z bin {iz}")
    xx = np.linspace(0, xs.max() * 1.05, 50)
    ax.plot(xx, slope * xx, color=col, lw=0.8, ls="--")

for t in TARGETS:
    ax.annotate(LABELS[t], (x[t], 1), xycoords=("data", "axes fraction"), xytext=(0, 4),
                textcoords="offset points", ha="center", fontsize=9)
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel(r"$A^{1/3}$")
ax.set_ylabel(r"$\Delta\langle p_T^2\rangle$ [GeV$^2$]")
ax.set_xlim(0, xs.max() * 1.1)
ax.legend(frameon=False, fontsize=9)
ax.text(0.02, 0.97, f"DIAGNOSTIC ({SAMPLES[S]['label']}, RG-A photon ID); common bins, D-error weights",
        transform=ax.transAxes, va="top", fontsize=8, color="gray")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/dpt2_vs_A13_{S}.{ext}", dpi=150)
print(f"\nwrote {OUT}/dpt2_vs_A13_{S}.pdf/.png")
