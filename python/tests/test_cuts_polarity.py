"""v1.0.0: config/cuts.json holds one block per torus polarity.

* /polarities/outbending must hold exactly the v0.2.8 cut values (no outbending cut changed).
* /polarities/inbending must hold the RG-D common aspects analysis note (v3.1) values:
  DC edges Table 5 (IB), v_z windows Table 6 (CuSn IB, CxC IB), the rest Table 7.
* cuts_for_polarity lifts one block to the top level and refuses anything else.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pi0.config import ConfigError, cuts_for_polarity

CONFIG_DIR = Path(__file__).resolve().parents[2] / "config"
DATA = Path(__file__).resolve().parent / "data"


def _cuts() -> dict:
    return json.loads((CONFIG_DIR / "cuts.json").read_text())


def _values(o):
    """Strip every _comment-style key, recursively: only the values the code reads remain."""
    if isinstance(o, dict):
        return {k: _values(v) for k, v in o.items() if not k.startswith("_")}
    if isinstance(o, list):
        return [_values(x) for x in o]
    return o


def test_outbending_values_are_exactly_those_of_v028():
    want = json.loads((DATA / "cuts_v028_outbending_values.json").read_text())
    ob = _cuts()["polarities"]["outbending"]
    for block in ("electron", "photon", "pairing", "vertex"):
        assert _values(ob[block]) == want[block], block


def test_the_polarity_blocks_replace_the_top_level_ones():
    c = _cuts()
    for block in ("electron", "photon", "pairing", "vertex"):
        assert block not in c, f"/{block} should only exist under /polarities since v1.0.0"
    assert "polarity" not in c["beam"]
    assert set(k for k in c["polarities"] if not k.startswith("_")) == {"outbending", "inbending"}


def test_inbending_electron_values_are_the_notes():
    e = _cuts()["polarities"]["inbending"]["electron"]
    # Table 5, IB regional averages
    assert (e["dc_edge_r1_cm"], e["dc_edge_r2_cm"], e["dc_edge_r3_cm"]) == (2.58, 2.51, 1.37)
    # Table 7: the same list and values as outbending
    ob = _cuts()["polarities"]["outbending"]["electron"]
    for k in ("chi2pid_min", "chi2pid_max", "min_momentum_gev", "pcal_lv_min_cm", "pcal_lw_min_cm"):
        assert e[k] == ob[k], k
    assert e["min_momentum_gev"] == 1.8
    assert _values(e["trigger"]) == _values(ob["trigger"])
    assert e["sampling_fraction"]["n_sigma"] == ob["sampling_fraction"]["n_sigma"] == 3.5


def test_inbending_vertex_windows_are_the_notes_table_6():
    t = _cuts()["polarities"]["inbending"]["vertex"]["targets"]
    assert t["cu"]["peaks"] == [{"mu_cm": -7.765, "sigma_cm": 0.321}]
    assert t["sn"]["peaks"] == [{"mu_cm": -2.812, "sigma_cm": 0.315}]
    assert t["cxc"]["peaks"] == [{"mu_cm": -7.67, "sigma_cm": 0.313}, {"mu_cm": -2.713, "sigma_cm": 0.307}]
    for k in ("cu", "sn", "cxc"):
        assert t[k]["rule"] == "corrected_peaks" and t[k]["correction_enabled"] is True and t[k]["n_sigma"] == 3.0
    assert (t["ld2"]["rule"], t["ld2"]["vz_min_cm"], t["ld2"]["vz_max_cm"]) == ("raw_window", -15.0, 5.0)


def test_inbending_photon_and_pairing_equal_outbending_for_now():
    p = _cuts()["polarities"]
    assert _values(p["inbending"]["photon"]) == _values(p["outbending"]["photon"])
    assert _values(p["inbending"]["pairing"]) == _values(p["outbending"]["pairing"])


def test_cuts_for_polarity_lifts_one_block():
    c = _cuts()
    ib = cuts_for_polarity(c, "inbending")
    assert ib["electron"]["dc_edge_r1_cm"] == 2.58
    assert "polarities" not in ib
    assert ib["beam"] == c["beam"] and ib["bsa"] == c["bsa"]
    ob = cuts_for_polarity(c, "outbending")
    assert ob["electron"]["dc_edge_r1_cm"] == 1.88
    assert ob["vertex"]["targets"]["cu"]["peaks"] == [{"mu_cm": -7.861, "sigma_cm": 0.415}]


@pytest.mark.parametrize("bad", ["", "sideways", "Outbending", None])
def test_cuts_for_polarity_refuses_an_unknown_polarity(bad):
    with pytest.raises(ConfigError):
        cuts_for_polarity(_cuts(), bad)


def test_cuts_for_polarity_refuses_a_missing_block():
    c = _cuts()
    del c["polarities"]["inbending"]["vertex"]
    with pytest.raises(ConfigError, match="vertex"):
        cuts_for_polarity(c, "inbending")
    cuts_for_polarity(c, "outbending")  # the other block is untouched
