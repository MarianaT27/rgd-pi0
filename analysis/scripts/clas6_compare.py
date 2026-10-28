"""
clas6_compare.py : RG-D pi0 multiplicity ratios in the bins of the CLAS6 (EG2, 5 GeV) pi0
paper, overlaid with its published values.

Reference: T. Mineeva et al. (CLAS), "Suppression of neutral pion production in deep-inelastic
scattering off nuclei with the CLAS detector", arXiv:2406.04513, Appendix B, Tables I and II
(analysis note v16, 2020-07-13, Table 3.9 for the bin edges). Their values are in
analysis/scripts/mineeva2024_tables.csv, transcribed from the arXiv v3 PDF (162 rows:
Table I 18 blocks x 6 pT2 bins, Table II 9 blocks x 6 z bins). Blocks were assigned in table
order (nu fastest). One header in the PDF does not follow that order: the 9th block of Table I
is labelled "2.2 < nu < 3.2, 0.5 < z < 0.6", a label already used by the 7th; by its position
it is 3.73 < nu < 4.25, 0.5 < z < 0.6, which is how it is stored.

Usage (from the repo root):
    python3 analysis/scripts/clas6_compare.py <binned_dir> <out_dir>
e.g.
    python3 analysis/scripts/clas6_compare.py analysis/clas6_q2nu/binned_full analysis/clas6_q2nu/compare

<binned_dir>/<T>.root (T = LD2, CxC, Cu, Sn) must come from stageB_bin with
    --grid-a config/binning/grid_A_q2_nu.json --grid-b config/binning/grid_B_z_pt2.json

Their two sets of 3D bins, and how each RG-D point is built:
  1. (nu, z, pT2), integrated over Q2 in 1.0-4.1 GeV2      (their Table I / Fig. 3)
  2. (Q2, nu, z), integrated over pT2 in 0-1.5 GeV2       (their Table II / Fig. 2)
For each of their bins, the same- and mixed-event m_gg spectra of every RG-D 4D bin inside it
are SUMMED, and the summed spectrum goes through the production extraction
(pi0.extract.extract_bin: sideband scale, subtraction, Gaussian fit, +-n_sigma window sum).
N_DIS is summed over the same (Q2, nu) cells. Then R = (Y_A/N_DIS_A) / (Y_D/N_DIS_D), with
the same error formula as pi0.ratio. Ratios are NOT averaged over bins.

Differences that remain (not corrected here):
  - their pT2 upper limit in set 2 is set per target from the pT2 tail (note sec. 3.6); here 1.5
  - their photon cuts: > 12 deg from the electron and E_pi0 > 0.5 GeV; here the RG-D cuts
  - their nu edge is written 4.25 in Table I and 4.2 in Table II; both are taken as 4.25
  - beam energy 5.014 GeV (theirs) vs 10.54 GeV, so the same (Q2, nu) bin is at a different y
  - nuclei: C vs C is direct; their Fe (A=56) and Pb (A=208) are drawn next to Cu (A=63)
    and Sn (A=120), which are different nuclei
  - error bars: statistical only, for both

Reads:   <binned_dir>/{LD2,CxC,Cu,Sn}.root, config/, config/binning/grid_A_q2_nu.json,
         analysis/scripts/mineeva2024_tables.csv
Writes:  <out_dir>/clas6_compare.csv                 every bin: their values and RG-D values
         <out_dir>/clas6_q2_nu_z.pdf/.png             set 2, R vs z, one panel per (Q2, nu)
         <out_dir>/clas6_nu_z_pt2.pdf/.png            set 1, R vs pT2, one panel per (z, nu)
Prints:  number of bins computed and failed per set and target
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "python")
from pi0 import io  # noqa: E402
from pi0.config import load_config  # noqa: E402
from pi0.extract import extract_bin  # noqa: E402

TABLES = Path("analysis/scripts/mineeva2024_tables.csv")
GRID_A_NU = Path("config/binning/grid_A_q2_nu.json")
OURS = ["CxC", "Cu", "Sn"]
THEIRS = {"CxC": "C", "Cu": "Fe", "Sn": "Pb"}          # drawn side by side
LABEL = {"CxC": "C", "Cu": "Cu", "Sn": "Sn"}
COLOR = {"CxC": "tab:blue", "Cu": "tab:orange", "Sn": "tab:green"}
Q2_RANGE = (1.0, 4.1)       # set 1 integrates Q2 over their full range
PT2_RANGE = (0.0, 1.5)      # set 2 integrates pT2 over their full range
EPS = 1e-6


def norm_nu(x):
    """Their Table II writes the top nu edge as 4.2; Table I and the note write 4.25."""
    return 4.25 if abs(x - 4.2) < EPS else x


def indices_in(edges, lo, hi):
    """Indices of the bins of `edges` that exactly tile [lo, hi]; error if they do not."""
    e = np.asarray(edges)
    i0 = np.where(np.abs(e - lo) < EPS)[0]
    i1 = np.where(np.abs(e - hi) < EPS)[0]
    if len(i0) != 1 or len(i1) != 1:
        sys.exit(f"[{lo}, {hi}] is not on the RG-D edges {list(e)}; the grids must contain the CLAS6 edges")
    return list(range(int(i0[0]), int(i1[0])))


def load_tables():
    rows = list(csv.DictReader(open(TABLES)))
    out = []
    for r in rows:
        d = {"set": r["set"]}
        for k in ("a_lo", "a_hi", "b_lo", "b_hi", "c_lo", "c_hi"):
            d[k] = float(r[k])
        if d["set"] == "nu_z_pt2":            # a = nu, b = z, c = pT2
            d["a_lo"], d["a_hi"] = norm_nu(d["a_lo"]), norm_nu(d["a_hi"])
        else:                                 # a = Q2, b = nu, c = z
            d["b_lo"], d["b_hi"] = norm_nu(d["b_lo"]), norm_nu(d["b_hi"])
        for t in ("C", "Fe", "Pb"):
            v = r[f"R_{t}"]
            d[t] = (float(v), float(r[f"stat_{t}"]), float(r[f"sys_{t}"])) if v else None
        out.append(d)
    return out


def bins_of(row, b):
    """RG-D 4D bins and Grid A cells inside one CLAS6 bin."""
    q2e, nue = b.grid_a.axes[0].edges, b.grid_a.axes[1].edges
    ze, pte = b.grid_b.axes[0].edges, b.grid_b.axes[1].edges
    if row["set"] == "nu_z_pt2":
        iq = indices_in(q2e, *Q2_RANGE)
        inu = indices_in(nue, row["a_lo"], row["a_hi"])
        iz = indices_in(ze, row["b_lo"], row["b_hi"])
        ip = indices_in(pte, row["c_lo"], row["c_hi"])
    else:
        iq = indices_in(q2e, row["a_lo"], row["a_hi"])
        inu = indices_in(nue, row["b_lo"], row["b_hi"])
        iz = indices_in(ze, row["c_lo"], row["c_hi"])
        ip = indices_in(pte, *PT2_RANGE)
    b4 = [b.encode_4d(a, n, z, p) for a in iq for n in inu for z in iz for p in ip]
    cells = sorted({a * b.n_xb + n for a in iq for n in inu})
    return b4, cells


def extract_sum(data, b4, cfg):
    """Production extraction on the summed spectrum of a set of 4D bins."""
    sel = np.asarray(b4)
    n_same = data.n_same[sel].sum(axis=0)
    n_mixed = data.n_mixed[sel].sum(axis=0)
    sums = {k: v[sel].sum(axis=0) for k, v in data.sums.items()}
    return extract_bin(-1, n_same, n_mixed, sums, cfg, box=None)


def ratio(ya, yd, nda, ndd):
    if not (ya.ok and yd.ok) or nda <= 0 or ndd <= 0 or ya.value <= 0 or yd.value <= 0:
        return np.nan, np.nan
    r = (ya.value / nda) / (yd.value / ndd)
    rel = np.sqrt((yd.error / yd.value) ** 2 + (ya.error / ya.value) ** 2 + 1.0 / ndd + 1.0 / nda)
    return r, r * rel


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: python3 analysis/scripts/clas6_compare.py <binned_dir> <out_dir>")
    bdir, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    cfg = load_config(Path("config"), grid_a=GRID_A_NU)
    b = cfg.binning
    if b.grid_a.axes[1].name != "nu":
        sys.exit(f"{GRID_A_NU} does not have a nu axis")
    data = {t: io.load(bdir / f"{t}.root", cfg, allow_unpublishable=True) for t in ["LD2"] + OURS}

    rows = load_tables()
    nfail = {}
    for row in rows:
        b4, cells = bins_of(row, b)
        yd = extract_sum(data["LD2"], b4, cfg)
        ndd = int(data["LD2"].n_dis[cells].sum())
        row["zbar"] = yd.abscissae.get("z", np.nan) if yd.ok else np.nan
        row["pt2bar"] = yd.abscissae.get("pt2", np.nan) if yd.ok else np.nan
        for t in OURS:
            ya = extract_sum(data[t], b4, cfg)
            nda = int(data[t].n_dis[cells].sum())
            row[t] = ratio(ya, yd, nda, ndd)
            row[t + "_why"] = ya.reason if not ya.ok else (yd.reason if not yd.ok else "ok")
            if not np.isfinite(row[t][0]):
                nfail[(row["set"], t)] = nfail.get((row["set"], t), 0) + 1

    print("=== RG-D ratios in the CLAS6 pi0 bins")
    for s in ("q2_nu_z", "nu_z_pt2"):
        n = sum(r["set"] == s for r in rows)
        print(f"  {s}: {n} bins; failed per target: "
              + ", ".join(f"{LABEL[t]} {nfail.get((s, t), 0)}" for t in OURS))

    with open(out / "clas6_compare.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["# their Table I: a=nu, b=z, c=pT2 (Q2 1.0-4.1); Table II: a=Q2, b=nu, c=z (pT2 0-1.5). "
                    "R stat-only for RG-D. reason = extraction status of the RG-D target bin"])
        w.writerow(["set", "a_lo", "a_hi", "b_lo", "b_hi", "c_lo", "c_hi",
                    "C_R", "C_stat", "C_sys", "Fe_R", "Fe_stat", "Fe_sys", "Pb_R", "Pb_stat", "Pb_sys",
                    "rgd_C_R", "rgd_C_err", "rgd_C_reason", "rgd_Cu_R", "rgd_Cu_err", "rgd_Cu_reason",
                    "rgd_Sn_R", "rgd_Sn_err", "rgd_Sn_reason", "rgd_zmean_LD2", "rgd_pt2mean_LD2"])
        for r in rows:
            theirs = []
            for t in ("C", "Fe", "Pb"):
                theirs += list(r[t]) if r[t] else ["", "", ""]
            ours = []
            for t in OURS:
                ours += [r[t][0], r[t][1], r[t + "_why"]]
            w.writerow([r["set"], r["a_lo"], r["a_hi"], r["b_lo"], r["b_hi"], r["c_lo"], r["c_hi"]]
                       + theirs + ours + [r["zbar"], r["pt2bar"]])

    def draw(ax, sub, xkey_lo, xkey_hi):
        for k, t in enumerate(OURS):
            th = THEIRS[t]
            xc = np.array([0.5 * (r[xkey_lo] + r[xkey_hi]) for r in sub])
            dx = 0.01 * (k - 1) * (max(r[xkey_hi] for r in sub) - min(r[xkey_lo] for r in sub))
            y = np.array([r[t][0] for r in sub])
            e = np.array([r[t][1] for r in sub])
            ok = np.isfinite(y)
            ax.errorbar(xc[ok] + dx, y[ok], yerr=e[ok], fmt="o", ms=3.5, color=COLOR[t], capsize=1.5,
                        lw=0.8, label=f"RG-D {LABEL[t]}")
            hv = [(xc[i], r[th][0], r[th][1]) for i, r in enumerate(sub) if r[th]]
            if hv:
                hx, hy, he = map(np.array, zip(*hv))
                ax.errorbar(hx + dx + 0.004 * (max(r[xkey_hi] for r in sub)), hy, yerr=he, fmt="s",
                            mfc="none", ms=4, color=COLOR[t], capsize=1.5, lw=0.8, label=f"CLAS6 {th}")
        ax.axhline(1.0, color="gray", lw=0.6, ls="--")
        ax.set_ylim(0.0, 1.6)

    # set 2: (Q2, nu, z) -> panels rows Q2, cols nu, x = z
    s2 = [r for r in rows if r["set"] == "q2_nu_z"]
    q2b = sorted({(r["a_lo"], r["a_hi"]) for r in s2})
    nub = sorted({(r["b_lo"], r["b_hi"]) for r in s2})
    fig, axs = plt.subplots(len(q2b), len(nub), figsize=(4 * len(nub), 3 * len(q2b)), sharex=True, sharey=True,
                            squeeze=False)
    for i, q in enumerate(q2b):
        for j, n in enumerate(nub):
            ax = axs[len(q2b) - 1 - i, j]
            draw(ax, [r for r in s2 if (r["a_lo"], r["a_hi"]) == q and (r["b_lo"], r["b_hi"]) == n], "c_lo", "c_hi")
            ax.set_title(rf"${q[0]}<Q^2<{q[1]}$, ${n[0]}<\nu<{n[1]}$", fontsize=9)
            if i == 0:
                ax.set_xlabel(r"$z$")
            if j == 0:
                ax.set_ylabel(r"$R_A$")
    axs[0, -1].legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle(r"$\pi^0$ $R_A$ in the CLAS6 (Q$^2$, $\nu$, $z$) bins, $p_T^2<1.5$ GeV$^2$; "
                 "RG-D (filled) vs Mineeva et al. arXiv:2406.04513 (open); stat. errors", fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"clas6_q2_nu_z.{ext}", dpi=140)
    plt.close(fig)

    # set 1: (nu, z, pT2) -> panels rows z, cols nu, x = pT2
    s1 = [r for r in rows if r["set"] == "nu_z_pt2"]
    nub = sorted({(r["a_lo"], r["a_hi"]) for r in s1})
    zb = sorted({(r["b_lo"], r["b_hi"]) for r in s1})
    fig, axs = plt.subplots(len(zb), len(nub), figsize=(4 * len(nub), 2.6 * len(zb)), sharex=True, sharey=True,
                            squeeze=False)
    for i, z in enumerate(zb):
        for j, n in enumerate(nub):
            ax = axs[i, j]
            draw(ax, [r for r in s1 if (r["a_lo"], r["a_hi"]) == n and (r["b_lo"], r["b_hi"]) == z], "c_lo", "c_hi")
            ax.set_title(rf"${n[0]}<\nu<{n[1]}$, ${z[0]}<z<{z[1]}$", fontsize=9)
            if i == len(zb) - 1:
                ax.set_xlabel(r"$p_T^2$ [GeV$^2$]")
            if j == 0:
                ax.set_ylabel(r"$R_A$")
    axs[0, -1].legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle(r"$\pi^0$ $R_A$ in the CLAS6 ($\nu$, $z$, $p_T^2$) bins, $1<Q^2<4.1$ GeV$^2$; "
                 "RG-D (filled) vs Mineeva et al. arXiv:2406.04513 (open); stat. errors", fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"clas6_nu_z_pt2.{ext}", dpi=140)
    plt.close(fig)
    print(f"\nwrote {out}/clas6_compare.csv, clas6_q2_nu_z.pdf/.png, clas6_nu_z_pt2.pdf/.png")


if __name__ == "__main__":
    main()
