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
