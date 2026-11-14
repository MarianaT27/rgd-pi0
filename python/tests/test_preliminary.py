"""The Preliminary label (v0.2.3) that replaced the DIAGNOSTIC wording."""

from __future__ import annotations

from types import SimpleNamespace

from pi0.io import Provenance, preliminary_label


def _data(**entries):
    return SimpleNamespace(prov_stage_a=Provenance(directory="provenance_stageA_[1 input(s)]", entries=entries))


def test_version_and_rga_fallback():
    d = _data(**{"stageA_skim.code_versions": "0.2.1", "gbt.fallback_used": "TRUE -- RG-A model"})
    assert preliminary_label(d) == "Preliminary (v0.2.1, RG-A photon ID)"


def test_no_fallback_no_rga():
    d = _data(**{"stageA_skim.code_versions": "0.2.1", "gbt.fallback_used": "false"})
    assert preliminary_label(d) == "Preliminary (v0.2.1)"


def test_two_files_two_versions_name_the_oldest():
    """v0.2.6: the label names the production (oldest version), not every version."""
    a = _data(**{"stageA_skim.code_versions": "0.2.1", "gbt.fallback_used": "TRUE"})
    b = _data(**{"stageA_skim.code_version": "0.2.2"})
    assert preliminary_label(a, b) == "Preliminary (v0.2.1, RG-A photon ID)"


def test_nothing_known():
    assert preliminary_label(_data()) == "Preliminary"


# v0.2.6: the stamp Stage A really writes is the full meson string (meson.build, pi0_code_version)
MESON = "{} (meson project version; source tree is not under VCS, so no commit is recorded)"


def test_real_meson_string_gives_only_the_number():
    d = _data(**{"stageA_skim.code_version": MESON.format("0.2.1"), "gbt.fallback_used": "TRUE"})
    assert preliminary_label(d) == "Preliminary (v0.2.1, RG-A photon ID)"


def test_mixed_production_is_named_by_its_oldest_version():
    a = _data(**{"stageA_skim.code_version": MESON.format("0.2.1"), "gbt.fallback_used": "TRUE"})
    b = _data(**{"stageA_skim.code_version": MESON.format("0.2.2")})
    assert preliminary_label(a, b) == "Preliminary (v0.2.1, RG-A photon ID)"
    assert preliminary_label(_data(**{"stageA_skim.code_versions": "0.2.1, 0.2.2"})) == "Preliminary (v0.2.1)"


def test_garbled_v023_csv_header_is_rebuilt():
    """The exact first header line pi0.ratio wrote from v0.2.3 to v0.2.5 for the v0.2.1 slims."""
    from pi0.io import preliminary_label_from_text
    old = ("# *** Preliminary (v0.2.1 (meson project version; source tree is not under VCS, "
           "v0.2.2 (meson project version; source tree is not under VCS, vso no commit is recorded), "
           "RG-A photon ID) -- provenance blockers present, see the prov lines ***")
    assert preliminary_label_from_text(old) == "Preliminary (v0.2.1, RG-A photon ID)"
    assert preliminary_label_from_text("# provenance: target CxC") == ""
