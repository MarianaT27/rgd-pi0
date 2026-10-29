"""
dpt2.py : pT broadening Delta<pT2> for one sample, one panel per z bin.

Usage (from the repo root):
    python3 analysis/scripts/dpt2.py binning_v2

Reads:   analysis/<sample>/<results>/broadening_{CxC,Cu,Sn}.csv
         config/binning/grid_B_z_pt2.json for the z bin ranges (must be the grid this sample used)
Writes:  analysis/<sample>/dpt2_vs_A13_<sample>.pdf/.png   Delta<pT2> vs A^(1/3), one panel per z bin,
                                                            targets in the repo colours, with a straight
                                                            line through the origin
Prints:  Delta<pT2> per target and z bin, and the fitted slope.

Only 3D bins that are good (finite delta, positive errors) in all three targets are used,
so every target is averaged over the same bins. Weights are 1/spt2_D^2 (the LD2 error),
which are the same for all targets in a given bin.
"""
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLES = {
    "diag3": {"dir": "analysis/diag3", "results": "results3",     "label": "3% sample"},
    "f10":   {"dir": "analysis/f10",   "results": "results_f10",  "label": "10% sample"},
    "full_z10": {"dir": "analysis/full_z10", "results": "results_full", "label": "full sample, 10 z bins"},
    "full_z10_v2": {"dir": "analysis/full_z10_v2", "results": "results_full", "label": "full sample, 10 z bins, pool from xB 0.06"},
    "full":  {"dir": "analysis/full",  "results": "results_full", "label": "full sample"},
    "ndis":  {"dir": "analysis/ndis",  "results": "results_full", "label": "full sample, inclusive N_DIS"},
    "binning_v2": {"dir": "analysis/binning_v2", "results": "results_full", "label": "full sample, binning v2"},
}
TARGETS = {"CxC": 12.011, "Cu": 63.546, "Sn": 118.71}
LABELS = {"CxC": "C", "Cu": "Cu", "Sn": "Sn"}
COLORS = {"LD2": "black", "CxC": "#1f77b4", "Cu": "#d62728", "Sn": "#2ca02c"}   # repo colours

if len(sys.argv) != 2 or sys.argv[1] not in SAMPLES:
    sys.exit(f"usage: python3 analysis/scripts/dpt2.py <{'|'.join(SAMPLES)}>")
S = sys.argv[1]
OUT = SAMPLES[S]["dir"]
RES = f"{OUT}/{SAMPLES[S]['results']}"


def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


d = {t: load(f"{RES}/broadening_{t}.csv") for t in TARGETS}

zedges = np.array(next(a["edges"] for a in json.load(open("config/binning/grid_B_z_pt2.json"))["axes"]
                       if a["name"] == "z"))
if max(int(v) for t in TARGETS for v in d[t]["i_z"]) >= len(zedges) - 1:
    sys.exit("the z bins in the csv go beyond config/binning/grid_B_z_pt2.json: wrong grid for this sample")

good = {}
for t in TARGETS:
    ok = np.isfinite(d[t]["delta"]) & (d[t]["sdelta"] > 0) & (d[t]["spt2_D"] > 0)
    good[t] = set(d[t]["bin3d"][ok].astype(int))
common = set.intersection(*good.values())

zbins = sorted(set(d["CxC"]["i_z"].astype(int)))
x = {t: A ** (1 / 3) for t, A in TARGETS.items()}
xs = np.array([x[t] for t in TARGETS])

print(f"=== Delta<pT2> [GeV2], sample {S} ({RES}), {len(common)} common 3D bins")
print("i_z  z range        n_bins   " + "   ".join(f"{LABELS[t]:>16s}" for t in TARGETS) + "   slope [GeV2]")
ncol = 3
nrow = int(np.ceil(len(zbins) / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow), sharex=True, sharey=True, squeeze=False)
for k, iz in enumerate(zbins):
    ax = axs.flat[k]
    vals, errs = [], []
    for t in TARGETS:
        b = d[t]
        sel = np.array([(int(kk) in common) and (int(z) == iz) for kk, z in zip(b["bin3d"], b["i_z"])])
        nb = sel.sum()
        w = 1.0 / b["spt2_D"][sel] ** 2
        vals.append(np.sum(w * b["delta"][sel]) / np.sum(w))
        errs.append(np.sqrt(np.sum(w ** 2 * b["sdelta"][sel] ** 2)) / np.sum(w))
    vals, errs = np.array(vals), np.array(errs)
    slope = np.sum(xs * vals / errs ** 2) / np.sum(xs ** 2 / errs ** 2)
    zr = f"[{zedges[iz]:.2f}, {zedges[iz + 1]:.2f})"
    print(f"{iz:3d}  {zr:13s}  {nb:6d}   " + "   ".join(f"{v:.5f} +- {e:.5f}" for v, e in zip(vals, errs))
          + f"   {slope:.5f}")
    for t, v, e in zip(TARGETS, vals, errs):
        ax.errorbar(x[t], v, yerr=e, fmt="o", color=COLORS[t], capsize=2, label=LABELS[t])
    xx = np.linspace(0, xs.max() * 1.05, 50)
    ax.plot(xx, slope * xx, color="gray", lw=0.8, ls="--")
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_title(rf"$z \in {zr}$", fontsize=10)
    ax.set_xlim(0, xs.max() * 1.1)
for ax in axs.flat[len(zbins):]:
    ax.axis("off")
for ax in axs[-1]:
    ax.set_xlabel(r"$A^{1/3}$")
for ax in axs[:, 0]:
    ax.set_ylabel(r"$\Delta\langle p_T^2\rangle$ [GeV$^2$]")
axs.flat[0].legend(frameon=False, fontsize=9, loc="upper left")
fig.suptitle(f"DIAGNOSTIC ({SAMPLES[S]['label']}, RG-A photon ID); common bins, D-error weights; "
             "dashed: line through the origin", fontsize=9, color="gray")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/dpt2_vs_A13_{S}.{ext}", dpi=150)
print(f"\nwrote {OUT}/dpt2_vs_A13_{S}.pdf/.png")
