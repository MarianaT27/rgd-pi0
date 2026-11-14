"""HERMES comparison (v0.2.4): grids, bin selection of the cuts, and the x_F provenance warning."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest

from pi0.config import load_config
from pi0.io import Provenance, check_provenance

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def hermes():
    spec = importlib.util.spec_from_file_location("hermes_ra", REPO / "analysis/scripts/hermes_ra.py")
    mod = importlib.util.module_from_spec(spec)
    cwd = os.getcwd()
    os.chdir(REPO)  # the script resolves config/ from the repo root
    try:
        spec.loader.exec_module(mod)
        cfg = load_config(REPO / "config", grid_a=REPO / mod.GRID_A, grid_b=REPO / mod.GRID_B)
    finally:
        os.chdir(cwd)
    return mod, cfg


def test_grids_have_the_hermes_edges(hermes):
    mod, cfg = hermes
    e = mod.axes_of(cfg.binning)
    assert list(e["q2"]) == [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 25.0]
    assert list(e["nu"]) == [4.0, 6.0, 8.0, 10.0]
    assert list(e["z"]) == [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.3]
    assert list(e["pt2"]) == [0.0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.4, 1.95, 5.0]
    assert cfg.binning.grid_a.axes[1].name == "nu"


def test_z_panel_cuts_nu_not_z(hermes):
    mod, cfg = hermes
    b = cfg.binning
    b4, cells = mod.selection(b, "z", 0, {"nu": 6.0})        # 0.1 < z < 0.2 is kept
    assert cells == sorted(a * b.n_xb + n for a in range(b.n_q2) for n in (1, 2))
    assert len(b4) == b.n_q2 * 2 * 1 * b.n_pt2


def test_nu_panel_cuts_z_not_nu(hermes):
    mod, cfg = hermes
    b = cfg.binning
    b4, cells = mod.selection(b, "nu", 0, {"z": 0.2})       # 4 < nu < 6 is kept
    assert cells == [a * b.n_xb + 0 for a in range(b.n_q2)]
    assert len(b4) == b.n_q2 * 1 * (b.n_z - 1) * b.n_pt2


def test_q2_panel_cuts_both(hermes):
    mod, cfg = hermes
    b = cfg.binning
    b4, cells = mod.selection(b, "q2", 2, {"nu": 6.0, "z": 0.2})
    assert cells == [2 * b.n_xb + 1, 2 * b.n_xb + 2]
    assert len(b4) == 2 * (b.n_z - 1) * b.n_pt2


def test_a_cut_off_the_edges_is_refused(hermes):
    mod, cfg = hermes
    with pytest.raises(SystemExit):
        mod.selection(cfg.binning, "z", 0, {"nu": 7.0})


def _blockers(xf):
    cfg = load_config(REPO / "config")
    prov = Provenance("provenance", {} if xf is None else {"pairs.xf_min": xf})
    return [b for b in check_provenance(prov, Provenance("a", {"gbt.fallback_used": "false"}), cfg)
            if b.key == "provenance/pairs.xf_min"]


def test_xf_cut_file_is_flagged_not_blocked():
    (w,) = _blockers("x_F > 0 -- --xf-min WAS GIVEN")
    assert not w.fatal


def test_no_xf_cut_is_not_flagged():
    assert _blockers("none (the default: no x_F cut)") == []
    assert _blockers(None) == []


# ---------------------------------------------------------------------------
# v0.2.6: the HERMES values (HEPData ins749249 v1, Tables 97 to 112)
# ---------------------------------------------------------------------------


def test_hermes_values_file_matches_the_grids(hermes):
    mod, cfg = hermes
    data = mod.load_hermes(REPO / mod.HERMES_CSV)
    e = mod.axes_of(cfg.binning)
    n_bins = {"nu": 9, "z": 9, "q2": 8, "pt2": 9}       # their tables; our pT2 grid adds 1.95-5
    assert sorted({n for n, _ in data}) == sorted(mod.HERMES_NUCLEI)
    for (n, var), pts in data.items():
        assert len(pts) == n_bins[var], (n, var)
        assert all(0.3 < r < 1.2 and 0 < s < 0.2 for _, r, s in pts), (n, var)
    # bin edges in the file are the HERMES grid edges (nu only up to 10 in our grid)
    import csv as _csv
    with open(REPO / mod.HERMES_CSV) as f:
        rows = list(_csv.DictReader(line for line in f if not line.startswith("#")))
    for r in rows:
        edges = list(e[r["var"]])
        lo, hi = float(r["lo"]), float(r["hi"])
        if r["var"] == "nu" and lo >= 10.0:
            continue
        assert any(abs(lo - x) < 1e-9 for x in edges) and any(abs(hi - x) < 1e-9 for x in edges), r


def test_hermes_spot_values():
    """Copied unchanged: Kr vs z first bin is 0.853 +- 0.013 in HEPData Table 106."""
    import csv as _csv
    with open(REPO / "analysis/scripts/hermes2007_pi0.csv") as f:
        rows = list(_csv.DictReader(line for line in f if not line.startswith("#")))
    kr = [r for r in rows if r["nucleus"] == "Kr" and r["var"] == "z"][0]
    assert (kr["lo"], kr["hi"], kr["R"], kr["stat"], kr["hepdata_table"]) == ("0.1", "0.2", "0.853", "0.013", "106")
