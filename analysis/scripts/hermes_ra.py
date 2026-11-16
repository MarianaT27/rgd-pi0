"""
hermes_ra.py : RG-D pi0 multiplicity ratio R_A in the 1D bins of HERMES, with the HERMES cuts.
RG-D only: the HERMES values are not drawn (yet).

Reference: A. Airapetian et al. (HERMES), "Hadronization in semi-inclusive deep-inelastic
scattering on nuclei", Nucl. Phys. B 780 (2007) 1, arXiv:0704.3270. Bin edges from their
Table A.1; cuts from Sec. 5: "Unless specified otherwise the data are shown with the following
constraints: nu > 6.0 GeV, z > 0.2, and x_F > 0". Each panel applies the cuts on the OTHER
variables, not on its own (as their Table A.1 shows: the nu bins start at 4 GeV, the z bins at 0.1):

    panel    x axis (HERMES edges)                       cuts
    R vs nu  nu  = 4, 6, 8, 10 GeV                       z > 0.2,             x_F > 0
    R vs z   z   = 0.1, 0.2, ..., 0.9, 1.3               nu > 6 GeV,          x_F > 0
    R vs Q2  Q2  = 1, 1.5, 2, 3, 4, 5, 6, 8, 25 GeV2     nu > 6 GeV, z > 0.2, x_F > 0
    R vs pT2 pT2 = 0, 0.05, 0.1, 0.3, ..., 1.95, 5 GeV2  nu > 6 GeV, z > 0.2, x_F > 0

The nu and z cuts are bin edges of the HERMES grids, so they are applied by selecting bins.
x_F > 0 is applied per pi0 pair in Stage B (stageB_bin --xf-min 0, v0.2.4); the script reads
the cut from the files' provenance and writes it in the title. The cuts apply to THIS plot only:
the grids and the x_F cut are used nowhere else.

Usage (from the repo root):
    python3 analysis/scripts/hermes_ra.py <binned_dir> <out_dir>
e.g.
    python3 analysis/scripts/hermes_ra.py analysis/hermes/binned analysis/hermes
(run_hermes_chain.sh runs Stage B and then this script.)

<binned_dir>/<T>.root (T = LD2, CxC, Cu, Sn) must come from stageB_bin with
    --grid-a config/binning/grid_A_q2_nu_hermes.json --grid-b config/binning/grid_B_z_pt2_hermes.json
    --xf-min 0

Method, as clas6_compare.py: for each HERMES bin, the same- and mixed-event m_gg spectra of every
RG-D 4D bin inside it are SUMMED and the summed spectrum goes through the production extraction
(pi0.extract.extract_bin). N_DIS is summed over the same (Q2, nu) cells. Then
R = (Y_A/N_DIS_A) / (Y_D/N_DIS_D), statistical errors as in clas6_compare.py. Each point is drawn
at the target's own mean of the variable (sideband-subtracted). For nu that mean is
<Q2> / (2 M <x_B>), because Stage B stores sums of Q2 and x_B but not of nu: an approximation.

Differences from HERMES that remain:
  - beam 10.54 GeV vs 27.6 GeV: RG-D covers nu < 8.74 GeV (nu = E - E' with the electron cut
    p > 1.8 GeV; y < 0.85 alone would allow 8.96), so only their first nu bins; with nu > 6 GeV the
    other panels use 6 < nu < 8.74 GeV, the top of the RG-D range (y > 0.57)
  - nuclei: C, Cu, Sn vs their He, Ne, Kr, Xe
  - pi0 identification, acceptance, radiative corrections: RG-D's (none applied here), not theirs

HERMES values (v0.2.6): analysis/scripts/hermes2007_pi0.csv, the 16 pi0 tables of HEPData
ins749249 v1 (He, Ne, Kr, Xe vs nu, z, Q2, pT2), drawn as open squares at their own mean,
stat. errors only (as ours).

Layout (v0.2.6): as their pi0 figure, 4 panels in one row with a shared y axis and no gaps;
axis ranges read off their figure by eye: y 0.4-1.2, nu 1-25 GeV (linear), z 0-1.1 (linear),
Q2 0.8-16 GeV2 (log), pT2 0.01-4 GeV2 (log). RG-D points above 1.2 fall outside.

Reads:   <binned_dir>/{LD2,CxC,Cu,Sn}.root, config/, the two HERMES grids,
         analysis/scripts/hermes2007_pi0.csv
Writes:  <out_dir>/hermes_ra.csv        every RG-D point: bin, cuts, mean, R, error, extraction status
         <out_dir>/hermes_ra.pdf/.png   1x4 panels: R_A vs nu, z, Q2, pT2, RG-D and HERMES
Prints:  every point, and the bins that failed per target
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker  # noqa: E402,F401

sys.path.insert(0, "python")
from pi0 import io  # noqa: E402
from pi0.config import load_config  # noqa: E402
from pi0.extract import extract_bin  # noqa: E402

GRID_A = Path("config/binning/grid_A_q2_nu_hermes.json")
GRID_B = Path("config/binning/grid_B_z_pt2_hermes.json")
OURS = ["CxC", "Cu", "Sn"]
LABEL = {"CxC": "C", "Cu": "Cu", "Sn": "Sn"}
COLOR = {"CxC": "#1f77b4", "Cu": "#d62728", "Sn": "#2ca02c"}   # standard RG-D target colours
HERMES_CSV = Path("analysis/scripts/hermes2007_pi0.csv")
HERMES_NUCLEI = ["He", "Ne", "Kr", "Xe"]
# none shared with RG-D (blue, red, green) or with CLAS6 Fe/Pb (orange, gold)
HERMES_COLOR = {"He": "#9467bd", "Ne": "#8c564b", "Kr": "#e377c2", "Xe": "#7f7f7f"}
M_P = 0.938272
NU_CUT = 6.0
Z_CUT = 0.2
EPS = 1e-6

# (variable, cuts on the other variables {variable: lower limit}, x label, x range, x scale, ticks)
# Ranges and ticks as the HERMES pi0 figure, read off it by eye (v0.2.6).
PANELS = [
    ("nu", {"z": Z_CUT}, r"$\nu$ [GeV]", (1.0, 25.0), "linear", [10, 20]),
    ("z", {"nu": NU_CUT}, r"$z$", (0.0, 1.1), "linear", [0.2, 0.6, 1.0]),
    ("q2", {"nu": NU_CUT, "z": Z_CUT}, r"$Q^2$ [GeV$^2$]", (0.8, 16.0), "log", [1, 10]),
    ("pt2", {"nu": NU_CUT, "z": Z_CUT}, r"$p_T^2$ [GeV$^2$]", (0.01, 4.0), "log", [0.1, 1]),
]
Y_RANGE = (0.4, 1.2)
CUT_TEXT = {"nu": rf"$\nu>{NU_CUT:g}$ GeV", "z": rf"$z>{Z_CUT:g}$"}


def load_hermes(path=HERMES_CSV):
    """HERMES pi0 R_A from the HEPData copy: {(nucleus, var): [(mean, R, stat), ...]}."""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            out.setdefault((r["nucleus"], r["var"]), []).append((float(r["mean"]), float(r["R"]), float(r["stat"])))
    return out


def axes_of(b):
    """Edges per variable, in the order of the 4D index (q2, nu, z, pt2)."""
    return {"q2": b.grid_a.axes[0].edges, "nu": b.grid_a.axes[1].edges,
            "z": b.grid_b.axes[0].edges, "pt2": b.grid_b.axes[1].edges}


def above(edges, cut):
    """Indices of the bins at or above `cut`, which must be one of the edges."""
    e = np.asarray(edges)
    i0 = np.where(np.abs(e - cut) < EPS)[0]
    if len(i0) != 1:
        sys.exit(f"cut {cut} is not an edge of {list(e)}; the HERMES grids must contain it")
    return list(range(int(i0[0]), len(e) - 1))


def selection(b, var, k, cuts):
    """RG-D 4D bins and Grid A cells for bin k of `var`, with `cuts` on the other variables."""
    edges = axes_of(b)
    idx = {v: list(range(len(e) - 1)) for v, e in edges.items()}
    for v, c in cuts.items():
        idx[v] = above(edges[v], c)
    idx[var] = [k]
    b4 = [b.encode_4d(a, n, z, p) for a in idx["q2"] for n in idx["nu"] for z in idx["z"] for p in idx["pt2"]]
    cells = sorted({a * b.n_xb + n for a in idx["q2"] for n in idx["nu"]})
    return b4, cells


def extract_sum(data, b4, cfg):
    """Production extraction on the summed spectrum of a set of 4D bins (as clas6_compare.py)."""
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


def mean_of(y, var):
    """The target's own sideband-subtracted mean of `var` in the bin (nu from <Q2>, <x_B>)."""
    if not y.ok:
        return np.nan
    a = y.abscissae
    if var == "nu":
        q2, xb = a.get("q2", np.nan), a.get("xb", np.nan)
        return q2 / (2 * M_P * xb) if xb and np.isfinite(xb) else np.nan
    return a.get(var, np.nan)


def compute(data, cfg):
    """Every point of every panel: list of dicts."""
    b = cfg.binning
    edges = axes_of(b)
    rows = []
    for var, cuts, *_ in PANELS:
        for k in range(len(edges[var]) - 1):
            b4, cells = selection(b, var, k, cuts)
            yd = extract_sum(data["LD2"], b4, cfg)
            ndd = int(data["LD2"].n_dis[cells].sum())
            for t in OURS:
                ya = extract_sum(data[t], b4, cfg)
                nda = int(data[t].n_dis[cells].sum())
                r, e = ratio(ya, yd, nda, ndd)
                why = ya.reason if not ya.ok else (yd.reason if not yd.ok else "ok")
                rows.append({"var": var, "lo": edges[var][k], "hi": edges[var][k + 1],
                             "cuts": " ".join(f"{v}>{c:g}" for v, c in cuts.items()),
                             "target": t, "mean": mean_of(ya, var), "R": r, "err": e,
                             "Y_A": ya.value if ya.ok else np.nan, "Y_D": yd.value if yd.ok else np.nan,
                             "N_DIS_A": nda, "N_DIS_D": ndd, "n_4d_bins": len(b4), "reason": why})
    return rows


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: python3 analysis/scripts/hermes_ra.py <binned_dir> <out_dir>")
    bdir, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    cfg = load_config(Path("config"), grid_a=GRID_A, grid_b=GRID_B)
    if cfg.binning.grid_a.axes[1].name != "nu":
        sys.exit(f"{GRID_A} does not have a nu axis")
    data = {t: io.load(bdir / f"{t}.root", cfg, allow_unpublishable=True) for t in ["LD2"] + OURS}

    xf = {t: (d.prov.get("pairs.xf_min") or "none").strip() for t, d in data.items()}
    if len({v.split(" -- ")[0] for v in xf.values()}) != 1:
        sys.exit(f"the four files were binned with different x_F cuts: {xf}")
    xf_cut = xf["LD2"].split(" -- ")[0]
    if xf_cut.lower().startswith("none"):
        print("WARNING: these files were binned WITHOUT an x_F cut (stageB_bin --xf-min 0 not given). "
              "The HERMES cut x_F > 0 is NOT applied; the title says so.")
        xf_text = r"no $x_F$ cut"
    else:
        xf_text = xf_cut.replace("x_F", r"$x_F$")
    label = io.preliminary_label(*data.values())

    rows = compute(data, cfg)

    print(f"=== RG-D pi0 R_A in the HERMES bins ({xf_cut}); {label}")
    print("var   bin            target   mean     R       error    N_DIS_D    status")
    for r in rows:
        print(f"{r['var']:4s}  {r['lo']:5.2f}-{r['hi']:<6.2f}  {r['target']:6s}  {r['mean']:7.3f}  "
              f"{r['R']:6.3f}  {r['err']:6.3f}  {r['N_DIS_D']:9d}  {r['reason']}")

    with open(out / "hermes_ra.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([f"# RG-D pi0 R_A in the HERMES 1D bins (arXiv:0704.3270 Table A.1), {xf_cut}; {label}; "
                    "stat. errors only; mean = target's own sideband-subtracted mean (nu from <Q2>/(2M<xB>))"])
        keys = ["var", "lo", "hi", "cuts", "target", "mean", "R", "err", "Y_A", "Y_D",
                "N_DIS_A", "N_DIS_D", "n_4d_bins", "reason"]
        w.writerow(keys)
        for r in rows:
            w.writerow([r[k] for k in keys])

    hermes = load_hermes()
    fig, axs = plt.subplots(1, 4, figsize=(15.0, 4.6), sharey=True, gridspec_kw={"wspace": 0.0})
    for ax, (var, cuts, xlabel, xr, scale, ticks) in zip(axs, PANELS):
        sub = [r for r in rows if r["var"] == var]
        ax.set_xscale(scale)
        for i, t in enumerate(OURS):
            pts = [(r["mean"], r["R"], r["err"]) for r in sub
                   if r["target"] == t and np.isfinite(r["R"]) and np.isfinite(r["mean"])]
            if not pts:
                continue
            x, y, e = map(np.array, zip(*pts))
            # small sideways offset between targets: multiplicative on a log axis
            x = x * (1 + 0.03 * (i - 1)) if scale == "log" else x + (i - 1) * 0.008 * (xr[1] - xr[0])
            ax.errorbar(x, y, yerr=e, fmt="o", ms=6, color=COLOR[t], capsize=2, lw=1.0,
                        label=f"RG-D {LABEL[t]}")
        for n in HERMES_NUCLEI:
            pts = hermes.get((n, var), [])
            if pts:
                x, y, e = map(np.array, zip(*pts))
                ax.errorbar(x, y, yerr=e, fmt="s", ms=6, mfc="none", mew=1.2, color=HERMES_COLOR[n],
                            capsize=2, lw=1.0, label=f"HERMES {n}")
        ax.axhline(1.0, color="gray", lw=0.8, ls="--")
        ax.set_xlim(*xr)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{v:g}" for v in ticks])
        if scale == "linear":
            ax.minorticks_on()
        else:
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_title(", ".join([CUT_TEXT[v] for v in cuts] + [xf_text]), fontsize=9.5)
    axs[0].set_ylabel(r"$R_A^{\pi^0}$", fontsize=13)
    axs[0].set_ylim(*Y_RANGE)
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=11)
    fig.suptitle(r"$\pi^0$ multiplicity ratio in the HERMES bins and cuts: RG-D (filled circles) vs HERMES "
                 r"(open squares, Airapetian et al., NPB 780 (2007) 1, HEPData ins749249); stat. errors",
                 fontsize=11)
    fig.text(0.99, 0.0, label, ha="right", va="top", fontsize=9, color="0.4")
    fig.tight_layout(rect=(0, 0.0, 1, 0.95))
    for ext in ("pdf", "png"):
        fig.savefig(out / f"hermes_ra.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)

    for t in OURS:
        bad = [r for r in rows if r["target"] == t and not np.isfinite(r["R"])]
        print(f"{LABEL[t]}: {len(bad)} of {sum(r['target'] == t for r in rows)} points failed"
              + ("" if not bad else ": " + "; ".join(f"{r['var']} {r['lo']:g}-{r['hi']:g} ({r['reason']})"
                                                       for r in bad)))
    print(f"\nwrote {out}/hermes_ra.csv, hermes_ra.pdf/.png")


if __name__ == "__main__":
    main()
