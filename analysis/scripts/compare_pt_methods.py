"""
compare_pt_methods.py : compare two ways of splitting the pi0 momentum into longitudinal and
transverse parts with respect to the virtual photon, on real slim events.

  Method A (what compute_sidis does now): projection on q-hat
      p_L = p . q_hat,   p_T^2 = |p - p_L q_hat|^2
  Method B (cross product, as suggested): transverse part from the cross product
      p_T^2 = |p x q_hat|^2,   p_L = sqrt(|p|^2 - p_T^2)   (as written in the Lambda note, no sign)

x_F is computed from each p_L with the same formula (Eq. 13 of the x_F document).

Usage (from the repo root):
    python3 analysis/scripts/compare_pt_methods.py <out_dir> slim1.root [slim2.root ...]
e.g.
    python3 analysis/scripts/compare_pt_methods.py analysis/ptcheck farm_ndis/LD2/018431/*/slim_0.root

Pairs: every pair of selected photons in an event with theta_e_gamma > 8 deg (the pairing cut
on single photons) and 0.10 < m_gg < 0.17 GeV. This is a comparison of two formulas on the same
vectors, so the exact pair selection does not matter; it only needs realistic kinematics.

Prints:  number of pairs, largest |A - B| for p_T^2, p_L and x_F, and the fraction of pairs
         with p_L < 0 in the lab (the only place the two p_L can differ)
Writes:  <out_dir>/pt_methods.pdf/.png   p_T^2, p_L and x_F distributions for A and B overlaid,
         and the per-pair differences B - A
"""
import json
import sys
from pathlib import Path

import awkward as ak
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

if len(sys.argv) < 3:
    sys.exit("usage: python3 analysis/scripts/compare_pt_methods.py <out_dir> slim1.root [slim2.root ...]")
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
FILES = sys.argv[2:]

cuts = json.load(open("config/cuts.json"))
EB = cuts["beam"]["energy_gev"]
MP, MPI = 0.938272, 0.1349768
EGAMMA_MIN = cuts["pairing"]["e_gamma_min_angle_deg"]

cols = ["nu", "w", "ex", "ey", "ez", "gpx", "gpy", "gpz", "g_e_gamma_deg"]
parts = [uproot.open(f)["events"].arrays(cols) for f in FILES]
ev = ak.concatenate(parts)

# photons passing the single-photon e-gamma cut, then all pairs per event
keep = ev.g_e_gamma_deg > EGAMMA_MIN
g = ak.zip({"x": ev.gpx[keep], "y": ev.gpy[keep], "z": ev.gpz[keep]})
pairs = ak.combinations(g, 2, fields=["a", "b"])
px = pairs.a.x + pairs.b.x
py = pairs.a.y + pairs.b.y
pz = pairs.a.z + pairs.b.z
e = np.sqrt(pairs.a.x**2 + pairs.a.y**2 + pairs.a.z**2) + np.sqrt(pairs.b.x**2 + pairs.b.y**2 + pairs.b.z**2)
m = np.sqrt(np.maximum(e**2 - px**2 - py**2 - pz**2, 0.0))
win = (m > 0.10) & (m < 0.17)

# virtual photon per event, broadcast to the pairs
qx, qy, qz = -ev.ex, -ev.ey, EB - ev.ez
qm = np.sqrt(qx**2 + qy**2 + qz**2)
ux, uy, uz = qx / qm, qy / qm, qz / qm
nu, w = ev.nu, ev.w
ux, uy, uz, nu, w, qm = (ak.broadcast_arrays(a, px)[0] for a in (ux, uy, uz, nu, w, qm))

sel = lambda a: ak.to_numpy(ak.flatten(a[win]))
px, py, pz, e, ux, uy, uz, nu, w, qm = (sel(a) for a in (px, py, pz, e, ux, uy, uz, nu, w, qm))
p2 = px**2 + py**2 + pz**2

# Method A: projection
pl_a = px * ux + py * uy + pz * uz
tx, ty, tz = px - pl_a * ux, py - pl_a * uy, pz - pl_a * uz
pt2_a = tx**2 + ty**2 + tz**2

# Method B: cross product
cx, cy, cz = py * uz - pz * uy, pz * ux - px * uz, px * uy - py * ux
pt2_b = cx**2 + cy**2 + cz**2
pl_b = np.sqrt(np.maximum(p2 - pt2_b, 0.0))


def xf(pl):
    a = w**2 + MPI**2 - MP**2
    return 2.0 * ((nu + MP) * pl - qm * e) / np.sqrt(a**2 - 4.0 * w**2 * MPI**2)


xf_a, xf_b = xf(pl_a), xf(pl_b)

n = len(pt2_a)
print(f"=== {len(FILES)} slim file(s), {len(ev)} events, {n} photon pairs with 0.10 < m_gg < 0.17 GeV")
print(f"max |pT2_B - pT2_A| = {np.max(np.abs(pt2_b - pt2_a)):.3e} GeV^2   (largest pT2 = {np.max(pt2_a):.3f})")
print(f"max |pL_B  - pL_A|  = {np.max(np.abs(pl_b - pl_a)):.3e} GeV")
print(f"max |xF_B  - xF_A|  = {np.nanmax(np.abs(xf_b - xf_a)):.3e}")
neg = pl_a < 0
print(f"pairs with p_L < 0 in the lab: {neg.sum()} of {n} ({100.0 * neg.mean():.4f}%)")
if neg.any():
    print(f"  excluding them: max |pL_B - pL_A| = {np.max(np.abs(pl_b[~neg] - pl_a[~neg])):.3e} GeV, "
          f"max |xF_B - xF_A| = {np.nanmax(np.abs(xf_b[~neg] - xf_a[~neg])):.3e}")

fig, ax = plt.subplots(2, 3, figsize=(13, 7))
for k, (va, vb, lab, rng) in enumerate([(pt2_a, pt2_b, r"$p_T^2$ [GeV$^2$]", (0, 1.5)),
                                        (pl_a, pl_b, r"$p_L$ (lab) [GeV]", (-1, 8)),
                                        (xf_a, xf_b, r"$x_F$", (-3, 1.5))]):
    bins = np.linspace(*rng, 151)
    d = vb - va
    d = d[np.isfinite(d)]
    va, vb = va[np.isfinite(va)], vb[np.isfinite(vb)]
    ax[0, k].hist(va, bins=bins, histtype="step", lw=1.5, label="A: projection")
    ax[0, k].hist(vb, bins=bins, histtype="step", lw=1.0, ls="--", label="B: cross product")
    ax[0, k].set_xlabel(lab)
    ax[0, k].set_yscale("log")
    ax[0, k].legend(frameon=False, fontsize=8)
    ax[1, k].hist(d, bins=200, histtype="step")
    ax[1, k].set_xlabel("B - A, " + lab)
    ax[1, k].set_yscale("log")
fig.suptitle(f"p_T / p_L with respect to q: projection vs cross product ({n} pairs)", fontsize=10)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"pt_methods.{ext}", dpi=130)
print(f"\nwrote {OUT}/pt_methods.pdf/.png")
