"""
make_grid_q2cut2.py : test copy of Grid A with a Q2 edge at the wagon's prescale threshold, 2.0 GeV^2.

The committed grid (config/binning/grid_A_q2_xb.json) has the Q2 bin [1.836, 2.225), which straddles
the RGDSidisWagon prescale threshold (Q2 < 2.0 prescaled by 5). Here the edge 1.836 is moved to 2.0,
so every Q2 bin is either entirely prescaled or entirely not:
    ... [1.574, 2.0), [2.0, 2.225) ...
All other edges, and the whole x_B axis, are unchanged. The committed grid is NOT modified.

Run from the repo root:  python3 analysis/scripts/make_grid_q2cut2.py
Writes analysis/cmp_wagon/grid_A_q2_xb_q2cut2.json
"""
import json
from pathlib import Path

SRC = Path("config/binning/grid_A_q2_xb.json")
OUT = Path("analysis/cmp_wagon/grid_A_q2_xb_q2cut2.json")
OLD, NEW = 1.8357110485458725, 2.0

g = json.loads(SRC.read_text())
q2 = next(a for a in g["axes"] if a["name"] == "q2")
i = min(range(len(q2["edges"])), key=lambda k: abs(q2["edges"][k] - OLD))
if abs(q2["edges"][i] - OLD) > 1e-9:
    raise SystemExit(f"edge {OLD} not found in {SRC}: {q2['edges']}")
q2["edges"][i] = NEW
g["_wagon_test"] = (f"TEST COPY of {SRC} for the RGDSidisWagon check: Q2 edge {OLD} moved to {NEW} (the wagon "
                    "prescale threshold) by analysis/scripts/make_grid_q2cut2.py. Everything else unchanged. "
                    "Not the production binning.")
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(g, indent=2) + "\n")
print("q2 edges:", [round(e, 4) for e in q2["edges"]])
print("wrote", OUT)
