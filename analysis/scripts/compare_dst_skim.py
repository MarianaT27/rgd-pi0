"""
compare_dst_skim.py : DST vs SIDIS-train skim for the same runs (output of run_dst_compare.sh).

Usage (from the repo root; run_dst_compare.sh calls it at the end):
    python3 analysis/scripts/compare_dst_skim.py

BEFORE running: run_dst_compare.sh must have produced
    analysis/cmp_dst/binned_full/<T>.root  and  analysis/cmp_skim/binned_full/<T>.root
and, for part 3, analysis/cmp_<dst|skim>/results_full/ratio_<T>.csv.

Prints:
  1. per target: total N_DIS and total same-event gamma-gamma pairs in the pi0 window
     (0.10 < m_gg < 0.17 GeV, before background subtraction), dst vs skim
  2. per Grid A cell: N_DIS(skim) / N_DIS(dst) for each target
  3. per z bin: R_A summed over all bins, dst vs skim

If the skim keeps every DIS event, N_DIS(skim)/N_DIS(dst) = 1 in every cell.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "python")
from pi0 import io  # noqa: E402
from pi0.config import load_config  # noqa: E402

TARGETS = ["LD2", "CxC", "Cu", "Sn"]
cfg = load_config(Path("config"))
mh = cfg.raw["mgg_histogram"]
edges = np.linspace(mh["min_gev"], mh["max_gev"], mh["bins"] + 1)
centres = 0.5 * (edges[1:] + edges[:-1])
win = (centres > 0.10) & (centres < 0.17)

d = {(s, t): io.load(Path(f"analysis/cmp_{s}/binned_full/{t}.root"), cfg, allow_unpublishable=True)
     for s in ("dst", "skim") for t in TARGETS}

print("\n=== 1. totals per target")
print("target   N_DIS dst      N_DIS skim     skim/dst   pairs dst     pairs skim    skim/dst   pairs/DIS dst  pairs/DIS skim")
for t in TARGETS:
    nd, ns = d["dst", t].n_dis.sum(), d["skim", t].n_dis.sum()
    pd, ps = d["dst", t].n_same[:, win].sum(), d["skim", t].n_same[:, win].sum()
    print(f"{t:6s}  {nd:12d}  {ns:12d}   {ns / nd:8.4f}  {pd:12d}  {ps:12d}   {ps / pd:8.4f}   {pd / nd:12.5f}  {ps / ns:12.5f}")

print("\n=== 2. N_DIS(skim) / N_DIS(dst) per Grid A cell (cell = i_q2 * n_xb + i_xb)")
b = d["dst", "LD2"].binning
print("cell  i_q2 i_xb   " + "   ".join(f"{t:>7s}" for t in TARGETS))
for c in range(b.n_a):
    vals = []
    for t in TARGETS:
        nd, ns = d["dst", t].n_dis[c], d["skim", t].n_dis[c]
        vals.append(f"{ns / nd:7.4f}" if nd > 0 else "    n/a")
    print(f"{c:4d}  {c // b.n_xb:4d} {c % b.n_xb:4d}   " + "   ".join(vals))


def load_csv(fn):
    rows = [l for l in open(fn) if not l.startswith("#")]
    names = rows[0].strip().split(",")
    data = np.array([[float(x) for x in l.strip().split(",")] for l in rows[1:] if l.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


print("\n=== 3. R_A summed over all bins of each z bin: dst vs skim")
print("target  i_z   R dst            R skim           skim/dst")
for t in ["CxC", "Cu", "Sn"]:
    r = {s: load_csv(f"analysis/cmp_{s}/results_full/ratio_{t}.csv") for s in ("dst", "skim")}
    for iz in sorted(set(r["dst"]["i_z"].astype(int))):
        out = []
        for s in ("dst", "skim"):
            x = r[s]
            k = (x["i_z"] == iz) & np.isfinite(x["R"]) & (x["Y_D"] > 0)
            f = x["N_DIS_A"][k] / x["N_DIS_D"][k]
            num, den = x["Y_A"][k].sum(), (x["Y_D"][k] * f).sum()
            err = np.sqrt(np.sum(x["sY_A"][k] ** 2) / den ** 2 + num ** 2 * np.sum((x["sY_D"][k] * f) ** 2) / den ** 4)
            out.append((num / den, err))
        print(f"{t:6s}  {iz:3d}   {out[0][0]:.3f} +- {out[0][1]:.3f}   {out[1][0]:.3f} +- {out[1][1]:.3f}   {out[1][0] / out[0][0]:.3f}")