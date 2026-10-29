"""
xf_z.py : background-subtracted pi0 distribution in (z, x_F), per target.

Purpose: see where each z bin sits in Feynman x, to judge the target-fragmentation
contribution (x_F < 0) as in the Lambda note (Chetry et al., Fig. 4.6 / App. S).

Usage (from the repo root):
    python3 analysis/scripts/xf_z.py <binned_dir> <out_dir>
e.g.
    python3 analysis/scripts/xf_z.py analysis/full/binned_full analysis/full/xf

BEFORE running: <binned_dir>/<T>.root for LD2, CxC, Cu, Sn must come from a
stageB_bin that writes the xf_z and xf_zb trees (x_F added Oct 2026).

Method, per target and per (z, x_F) cell, the same steps as the yield:
  1. alpha = same / mixed in the sideband (cuts.json /extraction), per cell;
     if the cell has no sideband counts, the alpha of its whole z row is used
  2. S = same - alpha * mixed, var = same + alpha^2 * mixed, per m_gg bin
  3. window = mu +- n_sigma_window * sigma from a Gaussian fit (pi0.extract.fit_core)
     of the subtracted spectrum summed over x_F in that z row; if that fit fails,
     the fit of the whole target is used
  4. cell yield = sum of S in the window, error = sqrt(sum of var in the window)
Only pairs in a 4D analysis bin enter (the same population as the yields).

Reads:   <binned_dir>/<T>.root  trees xf_z (50 z bins over [0,1]) and xf_zb (Grid B z bins)
         config/cuts.json, config/binning/grid_B_z_pt2.json
Writes:  <out_dir>/xf_z_map.pdf/.png      subtracted yield in (z, x_F), one panel per target
         <out_dir>/xf_per_zbin.pdf/.png   x_F distribution per Grid B z bin, targets overlaid (unit area)
Prints:  fraction of pi0 with x_F < 0 per Grid B z bin and target
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot

sys.path.insert(0, "python")
from pi0.config import load_config  # noqa: E402
from pi0.extract import fit_core, subtract  # noqa: E402

TARGETS = ["LD2", "CxC", "Cu", "Sn"]
LABELS = {"LD2": r"LD$_2$", "CxC": "C", "Cu": "Cu", "Sn": "Sn"}
COLORS = {"LD2": "black", "CxC": "#1f77b4", "Cu": "#d62728", "Sn": "#2ca02c"}
N_XF, XF_LO, XF_HI = 90, -3.0, 1.5   # must match stageB_bin (provenance xf_axis)
N_ZF = 50                            # fine z axis of tree xf_z, over [0, 1]

if len(sys.argv) != 3:
    sys.exit("usage: python3 analysis/scripts/xf_z.py <binned_dir> <out_dir>")
BIN_DIR, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

cfg = load_config(Path("config"))
ax, cuts = cfg.mgg, cfg.extraction
sb = ax.mask(cuts.sideband_min, cuts.sideband_max)
zb_edges = np.asarray(cfg.binning.grid_b.axes[0].edges)
xf_edges = np.linspace(XF_LO, XF_HI, N_XF + 1)
xf_c = 0.5 * (xf_edges[1:] + xf_edges[:-1])


def read(path, tree, nz):
    """Dense (nz, N_XF, n_mgg) same and mixed counts from one xf tree."""
    t = uproot.open(path)[tree].arrays(["iz", "ixf", "imgg", "n_same", "n_mixed"], library="np")
    same = np.zeros((nz, N_XF, ax.n))
    mixed = np.zeros((nz, N_XF, ax.n))
    same[t["iz"], t["ixf"], t["imgg"]] = t["n_same"]
    mixed[t["iz"], t["ixf"], t["imgg"]] = t["n_mixed"]
    return same, mixed


def window_of(s, var):
    """Boolean m_gg mask mu +- n_sigma * sigma from the core Gaussian fit, or None."""
    f = fit_core(ax, s, var, cuts)
    if not f.converged:
        return None
    c = ax.centres
    return (c >= f.mu - cuts.n_sigma_window * f.sigma) & (c < f.mu + cuts.n_sigma_window * f.sigma)


def subtracted(same, mixed):
    """Yield and variance per (z, x_F) cell, following the module docstring."""
    nz = same.shape[0]
    y = np.full((nz, N_XF), np.nan)
    v = np.full((nz, N_XF), np.nan)
    # target-wide fallback window
    a_all = same[..., sb].sum() / mixed[..., sb].sum()
    s_all, var_all = subtract(same.sum(axis=(0, 1)), mixed.sum(axis=(0, 1)), a_all)
    w_all = window_of(s_all, var_all)
    if w_all is None:
        sys.exit("the fit of the whole-target subtracted spectrum failed")
    for iz in range(nz):
        sb_same, sb_mixed = same[iz][:, sb].sum(), mixed[iz][:, sb].sum()
        if sb_mixed <= 0 or sb_same <= 0:
            continue
        a_row = sb_same / sb_mixed
        s_row, var_row = subtract(same[iz].sum(axis=0), mixed[iz].sum(axis=0), a_row)
        w = window_of(s_row, var_row)
        if w is None:
            w = w_all
        for ix in range(N_XF):
            n_sb_m, n_sb_s = mixed[iz, ix, sb].sum(), same[iz, ix, sb].sum()
            a = n_sb_s / n_sb_m if (n_sb_m > 0 and n_sb_s > 0) else a_row
            s, var = subtract(same[iz, ix], mixed[iz, ix], a)
            y[iz, ix] = s[w].sum()
            v[iz, ix] = var[w].sum()
    return y, v


res, resb = {}, {}
for t in TARGETS:
    path = BIN_DIR / f"{t}.root"
    same, mixed = read(path, "xf_z", N_ZF)
    res[t] = subtracted(same, mixed)
    nzb = len(zb_edges) - 1
    same_b, mixed_b = read(path, "xf_zb", nzb)
    resb[t] = subtracted(same_b, mixed_b)

# ---------------------------------------------------------------- 1. (z, x_F) maps
z_edges = np.linspace(0.0, 1.0, N_ZF + 1)
fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
for a, t in zip(axs.flat, TARGETS):
    y = res[t][0].copy()
    y[~(y > 0)] = np.nan  # empty, failed or negative cells left blank
    m = a.pcolormesh(xf_edges, z_edges, y, norm=LogNorm(), cmap="viridis")
    fig.colorbar(m, ax=a, label="subtracted yield")
    a.axvline(0, color="red", lw=1, ls="--")
    for e in zb_edges:
        a.axhline(e, color="white", lw=0.6, ls=":")
    a.set_title(LABELS[t])
for a in axs[1]:
    a.set_xlabel(r"$x_F$")
for a in axs[:, 0]:
    a.set_ylabel(r"$z$")
fig.suptitle(r"$\pi^0$ yield (mixed-event subtracted) in ($z$, $x_F$); dotted: Grid B $z$ edges", fontsize=10)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"xf_z_map.{ext}", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- 2. x_F per Grid B z bin
nzb = len(zb_edges) - 1
fig, axs = plt.subplots(1, nzb, figsize=(3.2 * nzb, 3.4), sharey=True)
axs = np.atleast_1d(axs)
print("=== fraction of pi0 with x_F < 0, per Grid B z bin (mixed-event subtracted)")
print("i_z   z range        " + "   ".join(f"{t:>14s}" for t in TARGETS))
for iz in range(nzb):
    row = []
    for t in TARGETS:
        y, v = resb[t][0][iz], resb[t][1][iz]
        ok = np.isfinite(y)
        tot = y[ok].sum()
        if tot > 0:
            axs[iz].errorbar(xf_c[ok], y[ok] / tot, yerr=np.sqrt(v[ok]) / tot, fmt="o", ms=2.5,
                             color=COLORS[t], label=LABELS[t], capsize=1, lw=0.8)
        neg, pos = ok & (xf_c < 0), ok & (xf_c > 0)
        a_, b_ = y[neg].sum(), y[pos].sum()
        va, vb = v[neg].sum(), v[pos].sum()
        if a_ + b_ > 0:
            f = a_ / (a_ + b_)
            ef = np.sqrt(b_**2 * va + a_**2 * vb) / (a_ + b_) ** 2
            row.append(f"{f:6.3f} +- {ef:5.3f}")
        else:
            row.append(f"{'n/a':>14s}")
    print(f"{iz:3d}   {zb_edges[iz]:.3f}-{zb_edges[iz + 1]:.3f}   " + "   ".join(row))
    axs[iz].axvline(0, color="gray", lw=0.8, ls="--")
    axs[iz].set_title(rf"${zb_edges[iz]:.3f} < z < {zb_edges[iz + 1]:.3f}$", fontsize=9)
    axs[iz].set_xlabel(r"$x_F$")
axs[0].set_ylabel("fraction per bin")
axs[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"xf_per_zbin.{ext}", dpi=150)
plt.close(fig)
print(f"\nwrote {OUT}/xf_z_map.pdf/.png and {OUT}/xf_per_zbin.pdf/.png")
