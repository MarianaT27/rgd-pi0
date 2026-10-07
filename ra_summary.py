import numpy as np
for t in ["CxC", "Cu", "Sn"]:
    lines = [l for l in open(f"ra3_{t}.csv") if not l.startswith("#")]
    d = np.genfromtxt(lines, delimiter=",", names=True)
    d = d[np.isfinite(d["R"])]
    def integ(s):
        num = np.sum(s["Y_A"])
        den = np.sum(s["Y_D"] * s["N_DIS_A"] / s["N_DIS_D"])
        return num / den
    print(f"\n{t}: {len(d)} bins, integrated R = {integ(d):.3f}, median R = {np.median(d['R']):.3f}")
    for iz in range(5):
        m = d["i_z"] == iz
        if m.sum() == 0:
            continue
        print(f"  z bin {iz} (<z> = {np.mean(d['z_mean'][m]):.2f}): R = {integ(d[m]):.3f}  median {np.median(d['R'][m]):.3f}  ({m.sum()} bins)")