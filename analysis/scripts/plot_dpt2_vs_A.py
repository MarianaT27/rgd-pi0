import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TARGETS = {"CxC": 12.011, "Cu": 63.546, "Sn": 118.71}
RES = "results_f10"

def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}

d = {t: load(f"{RES}/broadening_{t}.csv") for t in TARGETS}

# bins present (finite, positive errors) in all three targets
good = {}
for t in TARGETS:
    ok = np.isfinite(d[t]["delta"]) & (d[t]["sdelta"] > 0) & (d[t]["spt2_D"] > 0)
    good[t] = set(d[t]["bin3d"][ok].astype(int))
common = set.intersection(*good.values())

zbins = sorted(set(d["CxC"]["i_z"].astype(int)))
x = {t: A ** (1 / 3) for t, A in TARGETS.items()}

fig, ax = plt.subplots(figsize=(6.5, 4.8))
colors = plt.cm.viridis(np.linspace(0, 0.9, len(zbins)))
print("i_z  n_bins   " + "   ".join(f"{t:>16s}" for t in TARGETS) + "   slope [GeV^2]")
for iz, col in zip(zbins, colors):
    vals, errs = [], []
    nb = 0
    for t in TARGETS:
        b = d[t]
        sel = np.array([(int(k) in common) and (int(z) == iz) for k, z in zip(b["bin3d"], b["i_z"])])
        nb = sel.sum()
        w = 1.0 / b["spt2_D"][sel] ** 2
        vals.append(np.sum(w * b["delta"][sel]) / np.sum(w))
        errs.append(np.sqrt(np.sum(w ** 2 * b["sdelta"][sel] ** 2)) / np.sum(w))
    xs = np.array([x[t] for t in TARGETS])
    vals, errs = np.array(vals), np.array(errs)
    slope = np.sum(xs * vals / errs ** 2) / np.sum(xs ** 2 / errs ** 2)
    print(f"{iz:3d}  {nb:6d}   " + "   ".join(f"{v:.5f} +- {e:.5f}" for v, e in zip(vals, errs)) + f"   {slope:.5f}")
    ax.errorbar(xs, vals, yerr=errs, fmt="o", color=col, capsize=2, label=f"z bin {iz}")
    xx = np.linspace(0, xs.max() * 1.05, 50)
    ax.plot(xx, slope * xx, color=col, lw=0.8, ls="--")

for t in TARGETS:
    ax.annotate(t, (x[t], 1), xycoords=("data", "axes fraction"), xytext=(0, 4), textcoords="offset points", ha="center", fontsize=9)
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel(r"$A^{1/3}$")
ax.set_ylabel(r"$\Delta\langle p_T^2\rangle$ [GeV$^2$]")
ax.set_xlim(0, xs.max() * 1.1)
ax.legend(frameon=False, fontsize=9)
ax.text(0.02, 0.97, "DIAGNOSTIC (10% sample, RG-A photon ID); common bins, D-error weights",
        transform=ax.transAxes, va="top", fontsize=8, color="gray")
fig.tight_layout()
fig.savefig("dpt2_vs_A13.pdf")
fig.savefig("dpt2_vs_A13.png", dpi=150)
print("wrote dpt2_vs_A13.pdf / .png")
