import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TARGETS = ["LD2", "CxC", "Cu", "Sn"]
COLORS = {"LD2": "black", "CxC": "tab:blue", "Cu": "tab:orange", "Sn": "tab:red"}

def load(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}

def z_centers():
    # use <z> from the R_A output if it has a mean-z column, otherwise fall back to the bin index
    try:
        d = load("ra3_CxC.csv")
        for key in d:
            if "z" in key.lower() and ("mean" in key.lower() or "avg" in key.lower()):
                return [np.mean(d[key][d["i_z"] == iz]) for iz in range(5)], r"$\langle z \rangle$"
    except Exception:
        pass
    return list(range(5)), "z bin index"

def per_z(t):
    d = load(f"bsa3_{t}.csv")
    ok = np.isfinite(d["A_LU"]) & (d["sA_LU"] > 0) & (np.abs(d["A_LU"]) < 0.5)
    a, s, iz = d["A_LU"][ok], d["sA_LU"][ok], d["i_z"][ok].astype(int)
    m, e = [], []
    for k in range(5):
        w = 1.0 / s[iz == k]**2
        m.append(np.sum(w * a[iz == k]) / np.sum(w))
        e.append(1.0 / np.sqrt(np.sum(w)))
    return np.array(m), np.array(e)

z, xlabel = z_centers()
z = np.array(z, dtype=float)
dx = 0.006 if xlabel != "z bin index" else 0.06
res = {t: per_z(t) for t in TARGETS}

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 7), sharex=True,
                               gridspec_kw={"height_ratios": [3, 2], "hspace": 0.05})
for i, t in enumerate(TARGETS):
    m, e = res[t]
    ax1.errorbar(z + (i - 1.5) * dx, m, yerr=e, fmt="o", color=COLORS[t], label=t, capsize=2)
ax1.axhline(0, color="gray", lw=0.8, ls="--")
ax1.set_ylabel(r"$A_{LU}^{\sin\phi_h}$")
ax1.legend(frameon=False, ncol=4)
ax1.text(0.02, 0.95, "DIAGNOSTIC (10% sample, RG-A photon ID)", transform=ax1.transAxes,
         va="top", fontsize=8, color="gray")

mD, eD = res["LD2"]
for i, t in enumerate(["CxC", "Cu", "Sn"]):
    m, e = res[t]
    r = m / mD
    er = np.abs(r) * np.sqrt((e / m)**2 + (eD / mD)**2)
    ax2.errorbar(z + (i - 1) * dx, r, yerr=er, fmt="o", color=COLORS[t], capsize=2)
ax2.axhline(1, color="gray", lw=0.8, ls="--")
ax2.set_ylabel(r"$A_{LU}^{A} / A_{LU}^{D}$")
ax2.set_xlabel(xlabel)
ax2.set_ylim(-1.5, 3.5)

fig.savefig("bsa_vs_z.pdf", bbox_inches="tight")
fig.savefig("bsa_vs_z.png", dpi=150, bbox_inches="tight")
print("wrote bsa_vs_z.pdf / .png  (x axis:", xlabel + ")")
