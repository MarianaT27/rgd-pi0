import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

style = {"CxC": ("C", "o", "#1f77b4"), "Cu": ("Cu", "s", "#d62728"), "Sn": ("Sn", "^", "#2ca02c")}

fig, ax = plt.subplots(figsize=(6, 4.5))
for i, (t, (lab, mk, col)) in enumerate(style.items()):
    lines = [l for l in open(f"ra3_{t}.csv") if not l.startswith("#")]
    d = np.genfromtxt(lines, delimiter=",", names=True)
    d = d[np.isfinite(d["R"])]
    zs, Rs, eRs = [], [], []
    for iz in range(5):
        s = d[d["i_z"] == iz]
        if len(s) == 0:
            continue
        k = s["N_DIS_A"] / s["N_DIS_D"]
        num = np.sum(s["Y_A"])
        den = np.sum(s["Y_D"] * k)
        R = num / den
        eR = R * np.sqrt(np.sum(s["sY_A"]**2) / num**2 + np.sum((s["sY_D"] * k)**2) / den**2)
        zs.append(np.mean(s["z_mean"]) + (i - 1) * 0.006)
        Rs.append(R)
        eRs.append(eR)
    ax.errorbar(zs, Rs, yerr=eRs, fmt=mk, color=col, label=lab, capsize=3, ms=6)

ax.axhline(1.0, color="gray", ls="--", lw=1)
ax.set_xlabel(r"$\langle z \rangle$")
ax.set_ylabel(r"$R_A^{\pi^0}$")
ax.set_title(r"RG-D $\pi^0$ multiplicity ratio vs $z$")
ax.text(0.5, 0.04, "DIAGNOSTIC: stat. errors only, no corrections", transform=ax.transAxes, ha="center", fontsize=8, color="red")
ax.legend()
fig.tight_layout()
fig.savefig("ra_vs_z.pdf")
fig.savefig("ra_vs_z.png", dpi=150)
print("wrote ra_vs_z.pdf and ra_vs_z.png")