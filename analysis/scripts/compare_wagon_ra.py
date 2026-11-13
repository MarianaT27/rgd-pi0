"""
compare_wagon_ra.py : R_A with the LD2 from the SIDIS train vs the LD2 from the RGDSidisWagon test (run 18431).

Reads analysis/cmp_wagon/results/ratio_<T>_{train,wagon}.csv (run_wagon_compare.sh). For one nuclear
target the two R_A share Y_A and N_DIS_A, so R_A(wagon) / R_A(train) = M_D(train) / M_D(wagon) with
M_D = Y_D / N_DIS_D, the same number for CxC, Cu and Sn. So:
  1. N_DIS_D(wagon) / N_DIS_D(train) per Q2 bin: shows the loss and the prescale
     (expected 0.708 / 5 = 0.142 below Q2 = 2.0 and 0.708 above).
  2. M_D(wagon) / M_D(train) per Q2 bin (weighted mean over the bins of that Q2 row) and chi2/ndf vs 1.
     Expected 1 in every Q2 bin.
  3. R_A(wagon) / R_A(train) per target: same test seen through the published quantity.
The wagon events are a subset of the train events, so the two samples are correlated and the errors
used here (both added in quadrature) are conservative: chi2/ndf below 1 is expected if they agree.

Run from the repo root:  python3 analysis/scripts/compare_wagon_ra.py
Writes analysis/cmp_wagon/ra_wagon_vs_train.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = Path("analysis/cmp_wagon")
Q2_EDGES = next(a for a in json.loads((D / "grid_A_q2_xb_q2cut2.json").read_text())["axes"]
                if a["name"] == "q2")["edges"]
Q2_CUT = 2.0
TARGETS = [T for T in ["CxC", "Cu", "Sn"] if (Path("analysis/cmp_wagon/results") / f"ratio_{T}_wagon.csv").exists()]
COLORS = {"LD2": "black", "CxC": "#1f77b4", "Cu": "#d62728", "Sn": "#2ca02c"}  # repo target colours


def load_csv(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


def join(a, b):
    """Rows of a and b with the same bin4d, in the same order."""
    common, ia, ib = np.intersect1d(a["bin4d"], b["bin4d"], return_indices=True)
    return {k: v[ia] for k, v in a.items()}, {k: v[ib] for k, v in b.items()}, len(a["bin4d"]), len(b["bin4d"])


def q2_label(i):
    return f"[{Q2_EDGES[i]:.2f}, {Q2_EDGES[i + 1]:.2f})"


# ---- 1 and 2: the LD2 columns, taken from the CxC files (identical in all three)
t, w, nt, nw = join(load_csv(D / "results/ratio_CxC_train.csv"), load_csv(D / "results/ratio_CxC_wagon.csv"))
print(f"bins in the CxC table: train {nt}, wagon {nw}, in both {len(t['bin4d'])}")
iq = t["i_q2"].astype(int)
md_t, md_w = t["Y_D"] / t["N_DIS_D"], w["Y_D"] / w["N_DIS_D"]
rho = md_w / md_t
s_rho = rho * np.sqrt((t["sY_D"] / t["Y_D"]) ** 2 + (w["sY_D"] / w["Y_D"]) ** 2)

print(f"\n{'Q2 bin':18s} {'prescaled':>9s} {'N_DIS wagon/train':>18s} {'bins':>5s} "
      f"{'M_D wagon/train (weighted mean)':>32s} {'chi2/ndf vs 1':>14s}")
rows = []
for i in range(len(Q2_EDGES) - 1):
    m = iq == i
    if not m.any():
        continue
    # N_DIS per Grid A cell is the same for every z, pT2 bin of the cell: take one per cell
    _, first = np.unique(t["cell_a"][m], return_index=True)
    ndis = w["N_DIS_D"][m][first].sum() / t["N_DIS_D"][m][first].sum()
    wt = 1 / s_rho[m] ** 2
    mean = np.sum(wt * rho[m]) / wt.sum()
    err = 1 / np.sqrt(wt.sum())
    chi2 = np.sum(((rho[m] - 1) / s_rho[m]) ** 2)
    pre = "yes" if Q2_EDGES[i + 1] <= Q2_CUT + 1e-9 else "no"
    rows.append((i, mean, err))
    print(f"{q2_label(i):18s} {pre:>9s} {ndis:18.4f} {m.sum():5d} {mean:22.4f} +- {err:.4f} "
          f"{chi2:8.1f}/{m.sum():<4d}")

# ---- 3: R_A per target
print(f"\n{'target':8s} {'bins':>5s} {'R_A wagon/train (weighted mean)':>32s} {'chi2/ndf vs 1':>14s}")
ra = {}
for T in TARGETS:
    a, b, _, _ = join(load_csv(D / f"results/ratio_{T}_train.csv"), load_csv(D / f"results/ratio_{T}_wagon.csv"))
    r = b["R"] / a["R"]
    sr = r * np.sqrt((a["sR"] / a["R"]) ** 2 + (b["sR"] / b["R"]) ** 2)
    ok = np.isfinite(r) & np.isfinite(sr) & (sr > 0)
    wt = 1 / sr[ok] ** 2
    print(f"{T:8s} {ok.sum():5d} {np.sum(wt * r[ok]) / wt.sum():22.4f} +- {1 / np.sqrt(wt.sum()):.4f} "
          f"{np.sum(((r[ok] - 1) / sr[ok]) ** 2):8.1f}/{ok.sum():<4d}")
    ra[T] = (a["R"][ok], a["sR"][ok], b["R"][ok], b["sR"][ok])

# ---- plot
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [1.4, 1]})
order = np.lexsort((t["bin4d"], iq))
x = np.arange(len(order))
a1.errorbar(x, rho[order], yerr=s_rho[order], fmt="o", ms=2, color=COLORS["LD2"], elinewidth=0.6)
a1.axhline(1, color="gray", lw=0.8)
for i in range(len(Q2_EDGES) - 1):
    sel = np.where(iq[order] == i)[0]
    if not len(sel):
        continue
    lo, hi = sel[0] - 0.5, sel[-1] + 0.5
    if Q2_EDGES[i + 1] <= Q2_CUT + 1e-9:
        a1.axvspan(lo, hi, color="#ff7f0e", alpha=0.08, lw=0)
    a1.axvline(hi, color="gray", lw=0.5, ls=":")
    a1.text((lo + hi) / 2, 1.27, rf"$Q^2\in${q2_label(i)}", ha="center", va="top", fontsize=7, rotation=90)
a1.set_ylim(0.7, 1.3)
a1.set_xlabel("4D bin, grouped by $Q^2$ bin (shaded: prescaled in the wagon)")
a1.set_ylabel(r"$M_D$(wagon) / $M_D$(train),  $M_D = Y_D/N_{DIS,D}$")
a1.set_title("LD2 multiplicity per 4D bin, run 18431")
lim = [np.inf, -np.inf]
for T, (rt, st, rw, sw) in ra.items():
    a2.errorbar(rt, rw, xerr=st, yerr=sw, fmt="o", ms=2.5, color=COLORS[T], elinewidth=0.5, label=T)
    lim = [min(lim[0], rt.min(), rw.min()), max(lim[1], rt.max(), rw.max())]
a2.plot(lim, lim, color="gray", lw=0.8)
a2.set_xlabel(r"$R_A$ with LD2 from the SIDIS train")
a2.set_ylabel(r"$R_A$ with LD2 from the wagon")
a2.set_title(r"$R_A$ per 4D bin")
a2.legend()
fig.tight_layout()
fig.savefig(D / "ra_wagon_vs_train.png", dpi=150)
print(f"\nplot: {D / 'ra_wagon_vs_train.png'}")
