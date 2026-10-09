"""
ra.py : multiplicity ratio R_A for one sample (diag3, f10 or full).

Usage (from the repo root):
    python3 analysis/scripts/ra.py full

BEFORE running this script, these must already exist for the chosen sample:
  1. Stage A slims, make_grid (config/binning/grid_*.json) and Stage B
     (binned ROOT files, e.g. analysis/full/binned_full/<T>.root)
  2. pi0.ratio for CxC, Cu and Sn, which writes the input csv files:
       PYTHONPATH=python python3 -m pi0.ratio --target <binned>/<T>.root --ld2 <binned>/LD2.root \
           --config config --out <results>/ratio_<T>.csv --allow-unpublishable
     (run_full_chain.sh does steps 1 and 2 for the full sample)

Reads:   analysis/<sample>/<results>/ratio_{CxC,Cu,Sn}.csv
Writes:  analysis/<sample>/ra_vs_z_<sample>.pdf/.png      R_A vs <z>, integrated over Q2, xB, pT2
         analysis/<sample>/ra_vs_z_nu_<sample>.pdf/.png   same, in 3 nu panels (equal D DIS counts)
         analysis/<sample>/ra_compare_pip_<sample>.png    pi0 vs RG-D pi+ thesis, 4 (xB, Q2) panels
Prints:  integrated R_A per target and per z bin, and the numbers behind each plot.

R_A in a group of bins = sum(Y_A) / sum(Y_D * N_DIS_A / N_DIS_D)   (yields summed, not averaged)
<z> of a group = Y_A-weighted mean of z_mean.
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
TARGETS = ["CxC", "Cu", "Sn"]
LABELS = {"CxC": "C", "Cu": "Cu", "Sn": "Sn"}
COLORS = {"CxC": "tab:blue", "Cu": "tab:orange", "Sn": "tab:red"}
M_P = 0.938272

if len(sys.argv) != 2 or sys.argv[1] not in SAMPLES:
    sys.exit("usage: python3 analysis/scripts/ra.py <diag3|f10|full>")
S = sys.argv[1]
OUT = SAMPLES[S]["dir"]
RES = f"{OUT}/{SAMPLES[S]['results']}"
TAG = f"DIAGNOSTIC ({SAMPLES[S]['label']}, RG-A photon ID), stat. errors only"


def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    d = {n: data[:, i] for i, n in enumerate(names)}
    ok = np.isfinite(d["R"]) & (d["Y_D"] > 0)
    return {n: v[ok] for n, v in d.items()}


def ratio(d, sel):
    """summed-yield R_A, its stat. error and the Y_A-weighted <z> for the rows in sel"""
    ya, sya = d["Y_A"][sel], d["sY_A"][sel]
    yd, syd = d["Y_D"][sel], d["sY_D"][sel]
    f = d["N_DIS_A"][sel] / d["N_DIS_D"][sel]
    num, den = ya.sum(), (yd * f).sum()
    r = num / den
    err = np.sqrt(np.sum(sya**2) / den**2 + num**2 * np.sum((syd * f)**2) / den**4)
    zm = np.sum(ya * d["z_mean"][sel]) / num
    return zm, r, err


data = {t: load(f"{RES}/ratio_{t}.csv") for t in TARGETS}
zbins = sorted(set(data["CxC"]["i_z"].astype(int)))

# ---------------------------------------------------------------- summary
print(f"=== R_A summary, sample {S} ({RES})")
for t in TARGETS:
    d = data[t]
    _, r, e = ratio(d, np.ones(len(d["R"]), bool))
    print(f"\n{t}: {len(d['R'])} bins, integrated R = {r:.3f} +- {e:.3f}, median bin R = {np.median(d['R']):.3f}")
    for iz in zbins:
        m = d["i_z"] == iz
        zm, r, e = ratio(d, m)
        print(f"  z bin {iz} (<z> = {zm:.3f}): R = {r:.3f} +- {e:.3f}  median {np.median(d['R'][m]):.3f}  ({m.sum()} bins)")

# ---------------------------------------------------------------- R_A vs z
fig, ax = plt.subplots(figsize=(6, 4.5))
for i, t in enumerate(TARGETS):
    d = data[t]
    pts = [ratio(d, d["i_z"] == iz) for iz in zbins]
    zs, rs, es = map(np.array, zip(*pts))
    ax.errorbar(zs + (i - 1) * 0.006, rs, yerr=es, fmt="o", color=COLORS[t], label=LABELS[t], capsize=3, ms=6)
ax.axhline(1.0, color="gray", ls="--", lw=1)
ax.set_xlabel(r"$\langle z \rangle$")
ax.set_ylabel(r"$R_A^{\pi^0}$")
ax.set_title(r"RG-D $\pi^0$ multiplicity ratio vs $z$")
ax.text(0.5, 0.04, TAG, transform=ax.transAxes, ha="center", fontsize=8, color="gray")
ax.legend(frameon=False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/ra_vs_z_{S}.{ext}", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- R_A vs z in nu panels
NPANELS = 3
ref = data["CxC"]
cells = np.unique(ref["cell_a"]).astype(int)
nu_cell, ndis_cell = {}, {}
for c in cells:                       # nu of each Grid A cell from the D-yield-weighted Q2 and xB
    k = ref["cell_a"] == c
    w = ref["Y_D"][k]
    q2 = np.sum(w * ref["q2_mean"][k]) / np.sum(w)
    xb = np.sum(w * ref["xb_mean"][k]) / np.sum(w)
    nu_cell[c] = q2 / (2 * M_P * xb)
    ndis_cell[c] = ref["N_DIS_D"][k][0]
order = sorted(cells, key=lambda c: nu_cell[c])   # split cells into panels with ~equal D DIS counts
cum = np.cumsum([ndis_cell[c] for c in order]) / sum(ndis_cell.values())
panel_of = {c: min(int(f * NPANELS - 1e-9), NPANELS - 1) for c, f in zip(order, cum)}
edges = [(min(nu_cell[c] for c in cells if panel_of[c] == p), max(nu_cell[c] for c in cells if panel_of[c] == p))
         for p in range(NPANELS)]

print(f"\n=== R_A vs z in nu panels\npanel  nu range (GeV)   target   <z>    R_A     error")
fig, axes = plt.subplots(1, NPANELS, figsize=(4 * NPANELS, 4), sharey=True)
for p, ax in enumerate(axes):
    for i, t in enumerate(TARGETS):
        d = data[t]
        pan = np.array([panel_of.get(int(c), -1) for c in d["cell_a"]])
        zs, rs, es = [], [], []
        for iz in zbins:
            sel = (pan == p) & (d["i_z"] == iz)
            if sel.sum() == 0:
                continue
            zm, r, e = ratio(d, sel)
            zs.append(zm); rs.append(r); es.append(e)
            print(f"{p:5d}  {edges[p][0]:4.1f}-{edges[p][1]:4.1f}        {t:6s}  {zm:.3f}  {r:.3f}  {e:.3f}")
        ax.errorbar(np.array(zs) + (i - 1) * 0.006, rs, yerr=es, fmt="o", color=COLORS[t], label=LABELS[t], capsize=2)
    ax.axhline(1, color="gray", lw=0.8, ls="--")
    ax.set_title(rf"{edges[p][0]:.1f} < $\nu$ < {edges[p][1]:.1f} GeV")
    ax.set_xlabel(r"$\langle z \rangle$")
axes[0].set_ylabel(r"$R_A$")
axes[0].legend(frameon=False)
axes[0].text(0.02, 0.03, TAG, transform=axes[0].transAxes, fontsize=7, color="gray")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/ra_vs_z_nu_{S}.{ext}", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- pi0 vs pi+ thesis
# pi+ RG-D thesis values read off Fig. 4.1 at z ~ 0.34, 0.42, 0.50, 0.58, 0.66
PT2 = (0.24, 0.48)
ZR = (0.30, 0.70)
ZT = [0.34, 0.42, 0.50, 0.58, 0.66]
PANELS = [
    ("(a)", (0.11, 0.13), (1.1, 1.5), {"CxC": [.915, .895, .875, .85, .825], "Cu": [.755, .71, .695, .67, .63], "Sn": [.665, .615, .57, .55, .535]}),
    ("(b)", (0.20, 0.25), (1.2, 1.7), {"CxC": [.91, .89, .87, .85, .835], "Cu": [.755, .715, .685, .67, .65], "Sn": [.655, .62, .59, .545, .53]}),
    ("(c)", (0.25, 0.36), (1.6, 2.4), {"CxC": [.955, .89, .845, .85, .81], "Cu": [.78, .695, .655, .645, .62], "Sn": [.70, .615, .56, .54, .515]}),
    ("(d)", (0.36, 1.00), (1.8, 3.0), {"CxC": [1.03, .90, .83, .81, .81], "Cu": [.855, .715, .635, .60, .59], "Sn": [.76, .62, .565, .53, .51]}),
]
PIP_COLORS = {"CxC": "black", "Cu": "tab:green", "Sn": "tab:red"}


def inside(v, r):
    return (v >= r[0]) & (v < r[1])


print(f"\n=== pi0 vs pi+ thesis, pT2 {PT2} GeV2, z {ZR}")
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharey=True)
for ax, (name, xbr, q2r, thesis) in zip(axes.flat, PANELS):
    print(f"\nPanel {name}: xB {xbr}, Q2 {q2r}")
    print("  target  i_z  n_bins   <z>    R_pi0   err     R_pi+ (thesis, interpolated)")
    for t in TARGETS:
        d = data[t]
        base = (inside(d["xb_mean"], xbr) & inside(d["q2_mean"], q2r)
                & inside(d["pt2_mean"], PT2) & inside(d["z_mean"], ZR))
        zs, rs, es = [], [], []
        for iz in sorted(set(d["i_z"][base].astype(int))):
            sel = base & (d["i_z"] == iz)
            zm, r, e = ratio(d, sel)
            print(f"  {t:6s}  {iz:3d}  {sel.sum():6d}   {zm:.3f}  {r:.3f}  {e:.3f}   {np.interp(zm, ZT, thesis[t]):.3f}")
            zs.append(zm); rs.append(r); es.append(e)
        ax.plot(ZT, thesis[t], "^", mfc="none", color=PIP_COLORS[t], label=f"{LABELS[t]} pi+ (thesis)")
        if zs:
            ax.errorbar(zs, rs, yerr=es, fmt="o", color=PIP_COLORS[t], capsize=2, label=f"{LABELS[t]} pi0")
    ax.axhline(1, color="gray", lw=0.6, ls="--")
    ax.set_title(f"{name}  xB {xbr}, Q2 {q2r} GeV2", fontsize=9)
    ax.set_xlabel("z")
    ax.set_xlim(0.28, 0.72)
axes[0, 0].set_ylabel(r"$R_A$"); axes[1, 0].set_ylabel(r"$R_A$")
axes[0, 0].legend(fontsize=7, frameon=False, ncol=2)
fig.suptitle(f"RG-D: pi0 ({SAMPLES[S]['label']}, diagnostic) vs pi+ (thesis), pT2 0.24-0.48 GeV2", fontsize=10)
fig.tight_layout()
fig.savefig(f"{OUT}/ra_compare_pip_{S}.png", dpi=150)
plt.close(fig)

print(f"\nwrote {OUT}/ra_vs_z_{S}.pdf/.png, {OUT}/ra_vs_z_nu_{S}.pdf/.png, {OUT}/ra_compare_pip_{S}.png")
