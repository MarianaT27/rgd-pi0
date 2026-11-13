"""
compare_wagon_train.py : RGDSidisWagon skim vs SIDIS train, run 18431 (LD2), at the slim level.

Both slims were made by stageA_skim with the same config/cuts.json (sha256 d2e18f5f...):
    train : /work/clas12/mtenorio/pi0/stageA_v021/LD2/018431/*/slim_*.root
    wagon : /work/clas12/mtenorio/pi0/wagon_test/slim_wagon_018431.root

The wagon keeps every event with Q2 > 2.0 GeV^2 and, for 0.8 < Q2 < 2.0, only events with
RUN::config.event % 5 == 0. Stage A's cuts (Q2 > 1, W > 2, y < 0.85, ...) are tighter than the
wagon's (Q2 > 0.8, W > 1.7), so every train slim event should be:
    Q2 > 2.0           -> in the wagon slim
    Q2 < 2.0, ev%5==0  -> in the wagon slim
    Q2 < 2.0, ev%5!=0  -> not in the wagon slim
Events are matched by event number. Events with event number 0 (RUN::config not filled) are
counted separately: the wagon rejects events without a RUN::config row.

Prints the counts for each case and writes analysis/cmp_wagon/q2_wagon_vs_train.png.
Run from the repo root:  python3 analysis/scripts/compare_wagon_train.py
"""
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

Q2_CUT = 2.0
PRESCALE = 5
TRAIN = sorted(glob.glob("/work/clas12/mtenorio/pi0/stageA_v021/LD2/018431/*/slim_*.root"))
WAGON = "/work/clas12/mtenorio/pi0/wagon_test/slim_wagon_018431.root"
OUT = Path("analysis/cmp_wagon")
OUT.mkdir(parents=True, exist_ok=True)


def load(files):
    ev, q2 = [], []
    for f in files:
        t = uproot.open(f)["events"].arrays(["event", "q2"], library="np")
        ev.append(t["event"])
        q2.append(t["q2"])
    return np.concatenate(ev), np.concatenate(q2)


print("train slims:", len(TRAIN))
for f in TRAIN:
    print("   ", f)
ev_t, q2_t = load(TRAIN)
ev_w, q2_w = load([WAGON])
print(f"train events: {len(ev_t)}")
print(f"wagon events: {len(ev_w)}")

# event number 0 = RUN::config not filled; cannot be matched
z_t, z_w = ev_t == 0, ev_w == 0
print(f"event number 0: train {z_t.sum()}, wagon {z_w.sum()} (not matched below)")
ev_t, q2_t = ev_t[~z_t], q2_t[~z_t]
ev_w, q2_w = ev_w[~z_w], q2_w[~z_w]

print(f"duplicate event numbers: train {len(ev_t) - len(np.unique(ev_t))}, wagon {len(ev_w) - len(np.unique(ev_w))}")

in_w = np.isin(ev_t, ev_w)
hi = q2_t > Q2_CUT
keep = (ev_t % PRESCALE) == 0

print("\n=== train events, by where they should be")
rows = [
    (f"Q2 > {Q2_CUT}                    (expect all in wagon)", hi),
    (f"Q2 < {Q2_CUT}, event % {PRESCALE} == 0    (expect all in wagon)", ~hi & keep),
    (f"Q2 < {Q2_CUT}, event % {PRESCALE} != 0    (expect none in wagon)", ~hi & ~keep),
]
print(f"{'case':52s} {'train':>10s} {'in wagon':>10s} {'fraction':>9s}")
for label, m in rows:
    n, k = m.sum(), (m & in_w).sum()
    print(f"{label:52s} {n:10d} {k:10d} {k / n if n else float('nan'):9.5f}")

extra = ~np.isin(ev_w, ev_t)
print(f"\nwagon events not in the train slim: {extra.sum()}")

# events that broke the expectation, with their Q2 (are they at the 2.0 boundary?)
bad = (hi & ~in_w) | (~hi & keep & ~in_w) | (~hi & ~keep & in_w)
if bad.any():
    print(f"\nfirst {min(20, bad.sum())} train events that break the expectation (event, Q2, in wagon):")
    for e, q, w in list(zip(ev_t[bad], q2_t[bad], in_w[bad]))[:20]:
        print(f"   {e:12d}  {q:8.4f}  {bool(w)}")

# ---- Q2 distribution: train vs wagon, raw and with weight 5 below Q2_CUT
w_w = np.where(q2_w < Q2_CUT, PRESCALE, 1)
print(f"\nwagon weighted total (x{PRESCALE} below Q2 = {Q2_CUT}): {w_w.sum()}   train total: {len(ev_t)}   "
      f"ratio {w_w.sum() / len(ev_t):.5f}")

edges = np.linspace(1.0, 9.0, 41)
h_t, _ = np.histogram(q2_t, edges)
h_w, _ = np.histogram(q2_w, edges)
h_ww, _ = np.histogram(q2_w, edges, weights=w_w)
c = 0.5 * (edges[1:] + edges[:-1])

fig, (a, r) = plt.subplots(2, 1, figsize=(7, 7), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
a.stairs(h_t, edges, color="black", lw=1.5, label="SIDIS train")
a.stairs(h_w, edges, color="gray", lw=1.2, ls="--", label="wagon (raw)")
a.stairs(h_ww, edges, color="#ff7f0e", lw=1.2, label=f"wagon (x{PRESCALE} for $Q^2<{Q2_CUT}$)")
a.axvline(Q2_CUT, color="red", ls=":", lw=1)
a.set_yscale("log")
a.set_ylabel("events / bin")
a.set_title("Run 18431 (LD2), Stage A slims")
a.legend()
with np.errstate(divide="ignore", invalid="ignore"):
    r.errorbar(c, h_ww / h_t, yerr=np.sqrt(h_ww * np.where(c < Q2_CUT, PRESCALE, 1)) / h_t, fmt="o",
               ms=3, color="#ff7f0e", label="wagon weighted / train")
    r.plot(c, h_w / h_t, "s", ms=3, color="gray", label="wagon raw / train")
r.axhline(1.0, color="black", lw=0.8)
r.axhline(1.0 / PRESCALE, color="black", lw=0.8, ls="--")
r.axvline(Q2_CUT, color="red", ls=":", lw=1)
r.set_xlabel(r"$Q^2$ (GeV$^2$)")
r.set_ylabel("ratio")
r.set_ylim(0, 1.3)
r.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "q2_wagon_vs_train.png", dpi=150)
print(f"\nplot: {OUT / 'q2_wagon_vs_train.png'}")
