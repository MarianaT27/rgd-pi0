"""
bsa_grid.py : A_LU^{sin phi_h} per 4D bin vs z, as in note/figures/ALU_grid_perQ2.pdf.

One page per Q2 bin; rows are x_B bins, columns are pT2 bins; the four targets are
overlaid in each panel. Every point is one bin exactly as pi0.bsa wrote it (column A_LU):
no averaging over bins, no ratio, no extra cut.

Usage (from the repo root):
    python3 analysis/scripts/bsa_grid.py full

Reads:   the pi0.bsa csv files of the sample (bsa_{LD2,CxC,Cu,Sn}.csv)
         <results>/ratio_CxC.csv for each bin's <z> (count-weighted mean, matched by bin4d);
         a bin with no <z> there is placed at the centre of its z bin
         config/binning/grid_A_q2_xb.json and grid_B_z_pt2.json for the bin ranges
Writes:  analysis/<sample>/bsa_grid_<sample>.pdf          (all Q2 pages)
         analysis/<sample>/bsa_grid_<sample>_q2_<i>.png   (one per Q2 bin)
Prints:  number of bins plotted per target
"""
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

SAMPLES = {
    "diag3": {"dir": "analysis/diag3", "results": "results3", "bsa": "analysis/diag3/bsa3_{t}.csv",
              "label": "3% sample"},
    "full": {"dir": "analysis/full", "results": "results_full", "bsa": "analysis/full/results_full/bsa_{t}.csv",
             "label": "full sample"},
    "full_z10": {"dir": "analysis/full_z10", "results": "results_full",
                 "bsa": "analysis/full_z10/results_full/bsa_{t}.csv", "label": "full sample, 10 z bins"},
    "full_z10_v2": {"dir": "analysis/full_z10_v2", "results": "results_full",
                    "bsa": "analysis/full_z10_v2/results_full/bsa_{t}.csv",
                    "label": "full sample, 10 z bins, pool from xB 0.06"},
        "binning_v2": {"dir": "analysis/binning_v2", "results": "results_full", "bsa": "analysis/binning_v2/results_full/bsa_{t}.csv", "label": "full sample, binning v2"},
}
TARGETS = ["LD2", "CxC", "Cu", "Sn"]
LABELS = {"LD2": r"LD$_2$", "CxC": "C", "Cu": "Cu", "Sn": "Sn"}
COLORS = {"LD2": "black", "CxC": "#1f77b4", "Cu": "#d62728", "Sn": "#2ca02c"}
MARKERS = {"LD2": "D", "CxC": "o", "Cu": "s", "Sn": "^"}
YLIM = 0.3

if len(sys.argv) != 2 or sys.argv[1] not in SAMPLES:
    sys.exit(f"usage: python3 analysis/scripts/bsa_grid.py <{'|'.join(SAMPLES)}>")
S = sys.argv[1]
OUT = SAMPLES[S]["dir"]
RES = f"{OUT}/{SAMPLES[S]['results']}"


def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


def grid_edges(path):
    g = json.load(open(path))
    return {a["name"]: np.array(a["edges"]) for a in g["axes"]}


E = {**grid_edges("config/binning/grid_A_q2_xb.json"), **grid_edges("config/binning/grid_B_z_pt2.json")}
nq, nx, nz, npt = (len(E[k]) - 1 for k in ("q2", "xb", "z", "pt2"))

d = {t: load(SAMPLES[S]["bsa"].format(t=t)) for t in TARGETS}
for t in TARGETS:  # the grid on disk must be the one this sample was binned with
    for col, k in (("i_q2", "q2"), ("i_xb", "xb"), ("i_z", "z"), ("i_pt2", "pt2")):
        if d[t][col].size and d[t][col].max() >= len(E[k]) - 1:
            sys.exit(f"{t}: {col} goes beyond the {k} edges in config/binning; wrong grid for sample {S}")

# <z> per bin4d from ratio_CxC.csv
r = load(f"{RES}/ratio_CxC.csv")
zmean = {int(b): z for b, z in zip(r["bin4d"], r["z_mean"]) if np.isfinite(z)}

ok = {t: np.isfinite(d[t]["A_LU"]) & np.isfinite(d[t]["sA_LU"]) & (d[t]["sA_LU"] > 0) for t in TARGETS}
print(f"=== A_LU per bin, sample {S}")
for t in TARGETS:
    print(f"  {t:4s}: {ok[t].sum()} bins plotted of {len(ok[t])} in the csv")

zc = 0.5 * (E["z"][:-1] + E["z"][1:])
dz = 0.004  # horizontal offset between targets, display only

with PdfPages(f"{OUT}/bsa_grid_{S}.pdf") as pdf:
    for iq in range(nq):
        fig, axes = plt.subplots(nx, npt, figsize=(2.6 * npt, 1.9 * nx), sharex=True, sharey=True,
                                 squeeze=False)
        for ix in range(nx):
            for ip in range(npt):
                ax = axes[ix, ip]
                has = False
                for k, t in enumerate(TARGETS):
                    b = d[t]
                    m = ok[t] & (b["i_q2"] == iq) & (b["i_xb"] == ix) & (b["i_pt2"] == ip)
                    if not m.any():
                        continue
                    has = True
                    z = np.array([zmean.get(int(bb), zc[int(iz)]) for bb, iz in zip(b["bin4d"][m], b["i_z"][m])])
                    ax.errorbar(z + (k - 1.5) * dz, b["A_LU"][m], yerr=b["sA_LU"][m], fmt=MARKERS[t],
                                color=COLORS[t], ms=3, capsize=1.5, lw=0.8, label=LABELS[t])
                ax.axhline(0, color="gray", lw=0.6, ls="--")
                ax.set_xlim(0, 1)
                ax.set_ylim(-YLIM, YLIM)
                if not has:
                    ax.set_facecolor("0.95")
                if ix == 0:
                    ax.set_title(rf"$p_T^2\in[{E['pt2'][ip]:.3f},{E['pt2'][ip + 1]:.3f}]$", fontsize=8)
                if ip == 0:
                    ax.set_ylabel(rf"$x_B\in[{E['xb'][ix]:.3f},{E['xb'][ix + 1]:.3f}]$", fontsize=8)
                if ix == nx - 1:
                    ax.set_xlabel(r"$z$", fontsize=8)
                ax.tick_params(labelsize=7)
        handles = [plt.Line2D([], [], marker=MARKERS[t], color=COLORS[t], ls="", label=LABELS[t]) for t in TARGETS]
        fig.legend(handles=handles, loc="upper right", ncol=4, frameon=False, fontsize=9)
        fig.suptitle(rf"$A_{{LU}}^{{\sin\phi_h}}$ vs $z$,  $Q^2\in[{E['q2'][iq]:.2f},{E['q2'][iq + 1]:.2f}]$ GeV$^2$"
                     f"  (rows $x_B$, columns $p_T^2$)   DIAGNOSTIC, {SAMPLES[S]['label']}",
                     x=0.02, ha="left", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        pdf.savefig(fig)
        fig.savefig(f"{OUT}/bsa_grid_{S}_q2_{iq}.png", dpi=130)
        plt.close(fig)
print(f"wrote {OUT}/bsa_grid_{S}.pdf and {nq} png pages")