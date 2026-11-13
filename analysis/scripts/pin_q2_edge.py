"""
pin_q2_edge.py : move the Q2 edge of a Grid A file that is nearest to a given value onto that value.

Used after make_grid (equal-statistics edges) so that a Q2 edge sits exactly at the RGDSidisWagon
prescale threshold, 2.0 GeV^2: every Q2 bin is then either entirely prescaled or entirely not.
Only that one edge changes; the outer edges are never moved, and the x_B axis is untouched. The two
Q2 bins next to the moved edge are no longer exactly equal-statistics; the file records what was moved.

Usage (from the repo root):
    python3 analysis/scripts/pin_q2_edge.py <grid_A.json> [value, default 2.0]
The file is rewritten in place.
"""
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
value = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0

g = json.loads(path.read_text())
q2 = next(a for a in g["axes"] if a["name"] == "q2")
e = q2["edges"]
inner = range(1, len(e) - 1)
i = min(inner, key=lambda k: abs(e[k] - value))
if not (e[i - 1] < value < e[i + 1]):
    raise SystemExit(f"{value} does not lie between the neighbours of edge {i} ({e[i - 1]}, {e[i + 1]}): {e}")
old = e[i]
e[i] = value
g["_q2_edge_pinned"] = (f"Q2 edge {i} moved from {old} (make_grid, equal statistics) to {value} by "
                        f"analysis/scripts/pin_q2_edge.py, so that no Q2 bin straddles the RGDSidisWagon "
                        f"prescale threshold. The bins [{e[i - 1]:.4f}, {value}) and [{value}, {e[i + 1]:.4f}) "
                        f"are therefore not exactly equal-statistics.")
path.write_text(json.dumps(g, indent=2) + "\n")
print(f"q2 edge {i}: {old:.4f} -> {value}")
print("q2 edges:", [round(x, 4) for x in e])
