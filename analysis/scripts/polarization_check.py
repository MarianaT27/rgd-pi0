"""
polarization_check.py : diagnostics for the per-run beam polarization in the BSA (v0.2.1).

Usage (from the repo root, with PYTHONPATH=python):
    python3 analysis/scripts/polarization_check.py <binned_dir> <results_dir> <outdir>

BEFORE running this script, these must already exist:
  1. Stage B v0.2.1 files <binned_dir>/{LD2,CxC,Cu,Sn}.root (they carry the "runs" and "bsa_pol" trees)
  2. pi0.bsa with the per-run table, for each target:
       <results_dir>/bsa_<T>.csv
  3. pi0.bsa with one fixed P = 0.790 for every bin, for each target:
       python3 -m pi0.bsa ... --polarization 0.790 --polarization-err 0.015 --out <results_dir>/bsa_<T>_P0790.csv

Reads:   cuts.json /bsa/polarization/periods (config/cuts.json)
Writes:  <outdir>/pol_runs.pdf/.png     (a) events per run vs run number, coloured by Moller period,
                                            with |P| of the period on the right axis; one panel per target
         <outdir>/pol_bins.pdf/.png     (b) histogram of the per-bin <P> (count-weighted over the periods),
                                            all targets, the period values marked as dotted lines
         <outdir>/pol_scale.pdf/.png    (c) A_LU(per-run P) / A_LU(P = 0.790) vs 0.790 / <P>_bin per bin;
                                            the points must lie on the diagonal
Prints:  per target: runs, events, BSA entries per period; range of <P>; largest deviation in (c).
"""
import sys
from pathlib import Path

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pi0.config import load_config

if len(sys.argv) != 4:
    sys.exit("usage: python3 analysis/scripts/polarization_check.py <binned_dir> <results_dir> <outdir>")
BINNED, RES, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)

TARGETS = ["LD2", "CxC", "Cu", "Sn"]
LABELS = {"LD2": "D", "CxC": "C", "Cu": "Cu", "Sn": "Sn"}
COLORS = {"LD2": "black", "CxC": "#1f77b4", "Cu": "#d62728", "Sn": "#2ca02c"}   # repo colours
P_FIXED = 0.790

cfg = load_config(Path("config"))
periods = cfg.bsa.polarization_periods
PK = np.array([pp.p for pp in periods])
PCOL = plt.get_cmap("tab10")


def load_csv(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


# ---------------------------------------------------------------- (a) events per run
fig, axs = plt.subplots(len(TARGETS), 1, figsize=(12, 2.6 * len(TARGETS)), sharex=True, squeeze=False)
bin_p = {}
for k, t in enumerate(TARGETS):
    with uproot.open(BINNED / f"{t}.root") as f:
        if "runs" not in f or "bsa_pol" not in f:
            sys.exit(f"{BINNED / t}.root has no 'runs' or 'bsa_pol' tree: rerun Stage B v0.2.1")
        runs = f["runs"].arrays(library="np")
        bp = f["bsa_pol"].arrays(library="np")
    ax = axs[k, 0]
    for j, pp in enumerate(periods):
        sel = runs["period"] == j
        if sel.any():
            ax.bar(runs["run"][sel], runs["events"][sel], width=1.0, color=PCOL(j % 10))
    ax.set_ylabel(f"{LABELS[t]}: events / run")
    ax2 = ax.twinx()
    for j, pp in enumerate(periods):
        ax2.plot([pp.run_min - 0.5, pp.run_max + 0.5], [pp.p, pp.p], color=PCOL(j % 10), lw=2)
    ax2.set_ylim(0.70, 0.88)
    ax2.set_ylabel(r"$|P|$")

    # per-bin <P> from bsa_pol, and the printout
    w = np.zeros((cfg.binning.n_4d, len(periods)))
    np.add.at(w, (bp["bin4d"].astype(int), bp["period"].astype(int)), bp["counts"])
    n = w.sum(axis=1)
    bin_p[t] = (w @ PK)[n > 0] / n[n > 0]
    print(f"=== {LABELS[t]}: {len(runs['run'])} runs, {int(runs['events'].sum())} events, "
          f"<P> per bin from {bin_p[t].min():.4f} to {bin_p[t].max():.4f}")
    print("    period  runs           |P|      runs used   events        BSA entries")
    for j, pp in enumerate(periods):
        sel = runs["period"] == j
        print(f"    {j:4d}    {pp.run_min}-{pp.run_max}  {pp.p:.5f}  {int(sel.sum()):6d}   "
              f"{int(runs['events'][sel].sum()):12d}  {int(w[:, j].sum()):12d}")
axs[-1, 0].set_xlabel("run number")
fig.suptitle("Events per run, coloured by Moller period; lines: |P| of the period (right axis)", fontsize=10)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"pol_runs.{ext}", dpi=130)

# ---------------------------------------------------------------- (b) per-bin <P>
fig, ax = plt.subplots(figsize=(7, 4.5))
edges = np.linspace(0.73, 0.86, 66)
for t in TARGETS:
    ax.hist(bin_p[t], bins=edges, histtype="step", color=COLORS[t], lw=1.5, label=LABELS[t])
for pp in periods:
    ax.axvline(pp.p, color="gray", ls=":", lw=0.8)
ax.set_xlabel(r"$\langle P\rangle$ of the 4D bin")
ax.set_ylabel("4D bins")
ax.set_title(r"Per-bin $\langle P\rangle$ (dotted: $|P|$ of each period)", fontsize=10)
ax.legend(frameon=False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"pol_bins.{ext}", dpi=130)

# ---------------------------------------------------------------- (c) scale check
fig, axs = plt.subplots(1, len(TARGETS), figsize=(4 * len(TARGETS), 4), sharex=True, sharey=True, squeeze=False)
print(f"\n=== (c) A_LU(per-run P) / A_LU(P = {P_FIXED}) vs {P_FIXED}/<P>_bin")
for k, t in enumerate(TARGETS):
    ax = axs[0, k]
    a = load_csv(RES / f"bsa_{t}.csv")
    b = load_csv(RES / f"bsa_{t}_P0790.csv")
    common, ia, ib = np.intersect1d(a["bin4d"], b["bin4d"], return_indices=True)
    ok = b["A_LU"][ib] != 0
    ratio = a["A_LU"][ia][ok] / b["A_LU"][ib][ok]
    expect = P_FIXED / a["P_mean"][ia][ok]
    ax.plot(expect, ratio, "o", ms=3, color=COLORS[t], label=f"{LABELS[t]}, {ok.sum()} bins")
    lo, hi = 0.90, 1.08
    ax.plot([lo, hi], [lo, hi], color="gray", lw=0.8, ls="--")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(rf"${P_FIXED}/\langle P\rangle_{{bin}}$")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    print(f"    {LABELS[t]:3s} {ok.sum():5d} bins, max |ratio - expected| = {np.max(np.abs(ratio - expect)):.2e}"
          f"   (bins only in one file: {len(a['bin4d']) + len(b['bin4d']) - 2 * len(common)})")
axs[0, 0].set_ylabel(rf"$A_{{LU}}$(per-run $P$) / $A_{{LU}}$($P$ = {P_FIXED})")
fig.suptitle("Each bin divided by its own <P>: points on the dashed diagonal", fontsize=10)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"pol_scale.{ext}", dpi=130)
print(f"\nwrote {OUT}/pol_runs, pol_bins, pol_scale (.pdf/.png)")
