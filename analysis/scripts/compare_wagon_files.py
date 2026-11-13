"""
compare_wagon_files.py : wagon vs SIDIS train, file by file, for wagon outputs known to be complete.

Each argument is a Stage A slim of ONE wagon output file (one DST file). A DST file covers a
contiguous range of event numbers, so the train events in [min, max] of that slim are the ones
the wagon saw. For those, every train event must be:
    Q2 > 2.0           -> in the wagon slim
    Q2 < 2.0, ev%5==0  -> in the wagon slim
    Q2 < 2.0, ev%5!=0  -> not in the wagon slim
and every wagon event must be in the train slim.

Run from the repo root:
    python3 analysis/scripts/compare_wagon_files.py <wagon_slim_1.root> [<wagon_slim_2.root> ...]
"""
import glob
import sys

import numpy as np
import uproot

Q2_CUT = 2.0
PRESCALE = 5
TRAIN = sorted(glob.glob("/work/clas12/mtenorio/pi0/stageA_v021/LD2/018431/*/slim_*.root"))


def load(f):
    t = uproot.open(f)["events"].arrays(["event", "q2"], library="np")
    return t["event"], t["q2"]


ev_t, q2_t = map(np.concatenate, zip(*[load(f) for f in TRAIN]))
ok_all = True
for f in sys.argv[1:]:
    ev_w, q2_w = load(f)
    ev_w, q2_w = ev_w[ev_w > 0], q2_w[ev_w > 0]
    lo, hi_ev = ev_w.min(), ev_w.max()
    sel = (ev_t >= lo) & (ev_t <= hi_ev)
    e, q = ev_t[sel], q2_t[sel]
    in_w = np.isin(e, ev_w)
    hi = q > Q2_CUT
    keep = (e % PRESCALE) == 0
    rows = [
        (f"Q2 > {Q2_CUT}               (expect 1)", hi),
        (f"Q2 < {Q2_CUT}, ev % {PRESCALE} == 0  (expect 1)", ~hi & keep),
        (f"Q2 < {Q2_CUT}, ev % {PRESCALE} != 0  (expect 0)", ~hi & ~keep),
    ]
    extra = (~np.isin(ev_w, e)).sum()
    w = np.where(q2_w < Q2_CUT, PRESCALE, 1).sum()
    print(f"\n=== {f}")
    print(f"event range {lo} to {hi_ev}: train {len(e)}, wagon {len(ev_w)}")
    print(f"{'case':42s} {'train':>9s} {'in wagon':>9s} {'fraction':>9s}")
    for label, m in rows:
        n, k = m.sum(), (m & in_w).sum()
        print(f"{label:42s} {n:9d} {k:9d} {k / n if n else float('nan'):9.5f}")
    print(f"wagon events not in the train slim: {extra}")
    print(f"wagon weighted (x{PRESCALE} below {Q2_CUT}) / train: {w} / {len(e)} = {w / len(e):.5f}")
    bad = (hi & ~in_w) | (~hi & keep & ~in_w) | (~hi & ~keep & in_w)
    if bad.any() or extra:
        ok_all = False
        print("events that break the expectation (event, Q2, in wagon):")
        for x, y, z in list(zip(e[bad], q[bad], in_w[bad]))[:20]:
            print(f"   {x:12d}  {y:8.4f}  {bool(z)}")
print("\nALL FILES PASS" if ok_all else "\nSOME FILES FAIL (see above)")
