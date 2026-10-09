import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M_P = 0.938272
NPANELS = 3
TARGETS = ["CxC", "Cu", "Sn"]
COLORS = {"CxC": "tab:blue", "Cu": "tab:orange", "Sn": "tab:red"}

def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}

# nu of each Grid A cell, from the deuterium-yield-weighted mean Q2 and xB of its rows
ref = load("raf10_CxC.csv")
cells = np.unique(ref["cell_a"]).astype(int)
nu_cell, ndis_cell = {}, {}
for c in cells:
    k = ref["cell_a"] == c
    w = ref["Y_D"][k]
    q2 = np.sum(w * ref["q2_mean"][k]) / np.sum(w)
    xb = np.sum(w * ref["xb_mean"][k]) / np.sum(w)
    nu_cell[c] = q2 / (2 * M_P * xb)
    ndis_cell[c] = ref["N_DIS_D"][k][0]

# split cells into NPANELS nu ranges with ~equal deuterium DIS counts
order = sorted(cells, key=lambda c: nu_cell[c])
cum = np.cumsum([ndis_cell[c] for c in order]) / sum(ndis_cell.values())
panel_of = {c: min(int(f * NPANELS - 1e-9), NPANELS - 1) for c, f in zip(order, cum)}
edges = []
for p in range(NPANELS):
    nus = [nu_cell[c] for c in cells if panel_of[c] == p]
    edges.append((min(nus), max(nus)))

def ratio(d, sel):
    ya, sya = d["Y_A"][sel], d["sY_A"][sel]
    yd, syd = d["Y_D"][sel], d["sY_D"][sel]
    f = d["N_DIS_A"][sel] / d["N_DIS_D"][sel]
    num, den = ya.sum(), (yd * f).sum()
    r = num / den
    err = np.sqrt(np.sum(sya**2) / den**2 + num**2 * np.sum((syd * f)**2) / den**4)
    zm = np.sum(ya * d["z_mean"][sel]) / num
    return zm, r, err

fig, axes = plt.subplots(1, NPANELS, figsize=(4 * NPANELS, 4), sharey=True)
print("panel  nu range (GeV)   target   <z>    R_A     error")
for p, ax in enumerate(axes):
    for i, t in enumerate(TARGETS):
        d = load(f"raf10_{t}.csv")
        pan = np.array([panel_of.get(int(c), -1) for c in d["cell_a"]])
        zs, rs, es = [], [], []
        for iz in range(5):
            sel = (pan == p) & (d["i_z"] == iz) & (d["Y_D"] > 0)
            if sel.sum() == 0:
                continue
            zm, r, e = ratio(d, sel)
            zs.append(zm); rs.append(r); es.append(e)
            print(f"{p:5d}  {edges[p][0]:4.1f}-{edges[p][1]:4.1f}        {t:6s}  {zm:.3f}  {r:.3f}  {e:.3f}")
        ax.errorbar(np.array(zs) + (i - 1) * 0.006, rs, yerr=es, fmt="o",
                    color=COLORS[t], label=t, capsize=2)
    ax.axhline(1, color="gray", lw=0.8, ls="--")
    ax.set_title(rf"{edges[p][0]:.1f} < $\nu$ < {edges[p][1]:.1f} GeV")
    ax.set_xlabel(r"$\langle z \rangle$")
axes[0].set_ylabel(r"$R_A$")
axes[0].legend(frameon=False)
axes[0].text(0.02, 0.03, "DIAGNOSTIC (10% farm sample, RG-A photon ID)", transform=axes[0].transAxes,
             fontsize=8, color="gray")
fig.tight_layout()
fig.savefig("raf10_vs_z_nu.pdf")
fig.savefig("raf10_vs_z_nu.png", dpi=150)
print("wrote raf10_vs_z_nu.pdf / .png")
