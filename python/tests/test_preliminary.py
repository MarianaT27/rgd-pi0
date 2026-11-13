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


def test_two_files_two_versions_are_both_named():
    a = _data(**{"stageA_skim.code_versions": "0.2.1", "gbt.fallback_used": "TRUE"})
    b = _data(**{"stageA_skim.code_version": "0.2.2"})
    assert preliminary_label(a, b) == "Preliminary (v0.2.1, v0.2.2, RG-A photon ID)"


def test_nothing_known():
    assert preliminary_label(_data()) == "Preliminary"
