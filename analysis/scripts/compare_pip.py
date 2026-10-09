import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TARGETS = ["CxC", "Cu", "Sn"]
COLORS = {"CxC": "black", "Cu": "tab:green", "Sn": "tab:red"}
PT2 = (0.24, 0.48)
ZR = (0.30, 0.70)

# pi+ thesis (RG-D), values read off Fig. 4.1, z ~ 0.34, 0.42, 0.50, 0.58, 0.66
ZT = [0.34, 0.42, 0.50, 0.58, 0.66]
PANELS = [
    ("(a)", (0.11, 0.13), (1.1, 1.5), {"CxC": [.915, .895, .875, .85, .825], "Cu": [.755, .71, .695, .67, .63], "Sn": [.665, .615, .57, .55, .535]}),
    ("(b)", (0.20, 0.25), (1.2, 1.7), {"CxC": [.91, .89, .87, .85, .835], "Cu": [.755, .715, .685, .67, .65], "Sn": [.655, .62, .59, .545, .53]}),
    ("(c)", (0.25, 0.36), (1.6, 2.4), {"CxC": [.955, .89, .845, .85, .81], "Cu": [.78, .695, .655, .645, .62], "Sn": [.70, .615, .56, .54, .515]}),
    ("(d)", (0.36, 1.00), (1.8, 3.0), {"CxC": [1.03, .90, .83, .81, .81], "Cu": [.855, .715, .635, .60, .59], "Sn": [.76, .62, .565, .53, .51]}),
]

def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}

def inside(v, r):
    return (v >= r[0]) & (v < r[1])

data = {t: load(f"raf10_{t}.csv") for t in TARGETS}

fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharey=True)
for ax, (name, xbr, q2r, thesis) in zip(axes.flat, PANELS):
    print(f"\nPanel {name}: xB {xbr}, Q2 {q2r}, pT2 {PT2}")
    print("  target  i_z  n_bins   <z>    R_pi0   err     R_pi+ (thesis, interpolated)")
    for t in TARGETS:
        d = data[t]
        base = inside(d["xb_mean"], xbr) & inside(d["q2_mean"], q2r) & inside(d["pt2_mean"], PT2) & inside(d["z_mean"], ZR) & (d["Y_D"] > 0)
        zs, rs, es = [], [], []
        for iz in sorted(set(d["i_z"][base].astype(int))):
            sel = base & (d["i_z"] == iz)
            ya, sya, yd, syd = d["Y_A"][sel], d["sY_A"][sel], d["Y_D"][sel], d["sY_D"][sel]
            f = d["N_DIS_A"][sel] / d["N_DIS_D"][sel]
            num, den = ya.sum(), (yd * f).sum()
            r = num / den
            e = np.sqrt(np.sum(sya**2) / den**2 + num**2 * np.sum((syd * f)**2) / den**4)
            zm = np.sum(ya * d["z_mean"][sel]) / num
            pip = np.interp(zm, ZT, thesis[t])
            print(f"  {t:6s}  {iz:3d}  {sel.sum():6d}   {zm:.3f}  {r:.3f}  {e:.3f}   {pip:.3f}")
            zs.append(zm); rs.append(r); es.append(e)
        ax.plot(ZT, thesis[t], "^", mfc="none", color=COLORS[t], label=f"{t} pi+ (thesis)")
        if zs:
            ax.errorbar(zs, rs, yerr=es, fmt="o", color=COLORS[t], capsize=2, label=f"{t} pi0 (10%)")
    ax.axhline(1, color="gray", lw=0.6, ls="--")
    ax.set_title(f"{name}  xB {xbr}, Q2 {q2r} GeV2", fontsize=9)
    ax.set_xlabel("z")
    ax.set_xlim(0.28, 0.72)
axes[0, 0].set_ylabel(r"$R_A$"); axes[1, 0].set_ylabel(r"$R_A$")
axes[0, 0].legend(fontsize=7, frameon=False, ncol=2)
fig.suptitle("RG-D: pi0 (this analysis, 10%, diagnostic) vs pi+ (thesis), pT2 0.24-0.48 GeV2", fontsize=10)
fig.tight_layout()
fig.savefig("compare_pip.png", dpi=150)
print("\nwrote compare_pip.png")
