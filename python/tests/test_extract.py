"""Tests for the extraction stage, on synthetic spectra with known truth.

These tests exist to pin the four things the superseded analysis got wrong or
left undefined:

1. the yield is the **window sum**, not the Gaussian integral (its own docstring
   claimed the latter);
2. the abscissa is the **signal's** count-weighted mean, not the full-window
   mean and certainly not a box centre;
3. the quality gates actually reject;
4. ``bsa`` uses the measured per-run beam polarization and never a default.

Everything is built from an analytic model with known truth, so a failure here
means the code is wrong, not that the data is hard.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from pi0 import bsa, extract
from pi0.config import ConfigError, MggAxis, load_config, require

CONFIG_DIR = Path(__file__).resolve().parents[2] / "config"


@pytest.fixture(scope="module")
def cfg():
    """The real repository configuration -- no test-local cut values."""
    return load_config(CONFIG_DIR)


# --------------------------------------------------------------------------
# synthetic spectrum construction
# --------------------------------------------------------------------------
def gauss_counts(axis: MggAxis, n_total: float, mu: float, sigma: float) -> np.ndarray:
    """Counts per mgg bin of a Gaussian peak holding ``n_total`` entries."""
    m = axis.centres
    return n_total * axis.width / (sigma * math.sqrt(2 * math.pi)) * np.exp(
        -0.5 * ((m - mu) / sigma) ** 2
    )


def make_spectra(
    axis: MggAxis,
    *,
    n_signal: float,
    mu: float,
    sigma: float,
    mixed_per_bin: float,
    alpha_true: float,
    z_signal: float = 0.6,
    z_bkg: float = 0.3,
    tail_frac: float = 0.0,
    tail_offset_mult: float = 2.5,
    tail_sigma: float = 0.002,
    rng: np.random.Generator | None = None,
):
    """Build a same/mixed pair with a known signal on a known flat background.

    The background in the same-event spectrum is exactly ``alpha_true`` times
    the mixed spectrum, which is what the sideband scale must recover.

    Parameters
    ----------
    n_signal : float
        Entries in the injected core peak.
    mixed_per_bin : float
        Flat mixed-event level.
    alpha_true : float
        The true same/mixed background ratio.
    z_signal, z_bkg : float
        The ``<z>`` carried by signal pairs and by background pairs. Making
        them differ is the whole point of the abscissa test.
    tail_frac : float
        Entries in a **non-Gaussian** extra component, as a fraction of
        ``n_signal``, placed as a narrow shoulder at ``mu + tail_offset_mult *
        sigma``. It sits OUTSIDE the core fit range ``[mu - 1.5s, mu + 1.5s]``
        -- so the Gaussian fit cannot absorb it and its analytic integral cannot
        contain it -- but INSIDE the +-3 sigma window, so the window sum must.
        That separation is what makes the two estimators distinguishable.
    tail_offset_mult, tail_sigma : float
        Where the shoulder sits (in units of the core sigma) and how wide it is.
    rng : numpy.random.Generator, optional
        If given, Poisson noise is applied. If ``None`` the spectra are exact,
        so estimators can be tested to numerical precision.

    Returns
    -------
    n_same, n_mixed, sums, truth : tuple
    """
    core = gauss_counts(axis, n_signal, mu, sigma)
    tail = (
        gauss_counts(axis, n_signal * tail_frac, mu + tail_offset_mult * sigma, tail_sigma)
        if tail_frac
        else np.zeros_like(core)
    )
    signal = core + tail
    n_mixed = np.full(axis.n, float(mixed_per_bin))
    bkg = alpha_true * n_mixed
    n_same = signal + bkg
    if rng is not None:
        n_mixed = rng.poisson(n_mixed).astype(float)
        n_same = rng.poisson(n_same).astype(float)
        # Recompute the truth-carrying decomposition on the noisy spectra by
        # scaling the noiseless shares -- the sums are accumulators, not counts.
        share = np.divide(signal, signal + bkg, out=np.zeros_like(signal),
                          where=(signal + bkg) > 0)
        sig_part, bkg_part = n_same * share, n_same * (1 - share)
    else:
        sig_part, bkg_part = signal, bkg

    sums = {
        "q2": sig_part * 2.0 + bkg_part * 2.0,   # same in signal and background
        "xb": sig_part * 0.2 + bkg_part * 0.2,
        "z": sig_part * z_signal + bkg_part * z_bkg,   # DIFFERENT -- the test
        "pt2": sig_part * 0.15 + bkg_part * 0.4,       # DIFFERENT too
    }
    truth = {
        "signal": signal, "bkg": bkg, "z_signal": z_signal, "z_bkg": z_bkg,
        "mu": mu, "sigma": sigma, "alpha": alpha_true, "n_signal": n_signal,
    }
    return n_same, n_mixed, sums, truth


# --------------------------------------------------------------------------
# 1. the sideband scale recovers a known alpha
# --------------------------------------------------------------------------
def test_sideband_scale_recovers_known_alpha(cfg):
    """alpha = sum(same in SB)/sum(mixed in SB) must return the injected ratio."""
    axis = cfg.mgg
    n_same, n_mixed, _, truth = make_spectra(
        axis, n_signal=20000, mu=0.135, sigma=0.012, mixed_per_bin=2000.0, alpha_true=0.05
    )
    sb = axis.mask(cfg.extraction.sideband_min, cfg.extraction.sideband_max)
    alpha, why = extract.sideband_scale(n_same, n_mixed, sb)
    assert why == extract.OK
    # The recovered alpha is biased slightly HIGH by the signal's tail leaking
    # into the sideband -- a real, known property of a single high-mass
    # sideband, not a bug. It is well under a percent here.
    assert alpha == pytest.approx(truth["alpha"], rel=0.01)
    assert alpha > truth["alpha"], "signal leakage into the SB can only bias alpha up"


def test_sideband_scale_is_exact_with_no_signal(cfg):
    """With no signal at all, alpha must be recovered exactly."""
    axis = cfg.mgg
    n_mixed = np.full(axis.n, 1234.0)
    n_same = 0.037 * n_mixed
    sb = axis.mask(cfg.extraction.sideband_min, cfg.extraction.sideband_max)
    alpha, why = extract.sideband_scale(n_same, n_mixed, sb)
    assert why == extract.OK
    assert alpha == pytest.approx(0.037, rel=1e-12)


def test_sideband_scale_undefined_when_no_mixed(cfg):
    """An empty mixed sideband means no background estimate at all -- not alpha = 0."""
    axis = cfg.mgg
    sb = axis.mask(cfg.extraction.sideband_min, cfg.extraction.sideband_max)
    alpha, why = extract.sideband_scale(np.ones(axis.n), np.zeros(axis.n), sb)
    assert not np.isfinite(alpha)
    assert why == extract.NO_SIDEBAND_MIXED


# --------------------------------------------------------------------------
# 2. the yield recovers a known injected signal within statistics
# --------------------------------------------------------------------------
def test_yield_recovers_injected_gaussian_on_flat_background(cfg):
    """A known Gaussian on a known flat background: Y must recover the injection."""
    axis = cfg.mgg
    rng = np.random.default_rng(20260716)
    n_signal, mu, sigma = 20000.0, 0.135, 0.012
    n_same, n_mixed, sums, _ = make_spectra(
        axis, n_signal=n_signal, mu=mu, sigma=sigma, mixed_per_bin=2000.0,
        alpha_true=0.05, rng=rng,
    )
    by = extract.extract_bin(0, n_same, n_mixed, sums, cfg)
    assert by.ok, by.reason
    assert by.fit.converged
    assert by.fit.mu == pytest.approx(mu, abs=0.001)
    assert by.fit.sigma == pytest.approx(sigma, abs=0.001)

    # The +-3 sigma window holds erf(3/sqrt2) = 99.73% of a Gaussian.
    expected = n_signal * math.erf(cfg.extraction.n_sigma_window / math.sqrt(2))
    pull = (by.value - expected) / by.error
    assert abs(pull) < 4.0, f"yield {by.value:.1f} vs expected {expected:.1f}, pull {pull:.2f}"
    # And the quoted error must be of a sane size: sqrt(S + alpha^2 B)-ish.
    assert 0.003 < by.error / by.value < 0.03


# --------------------------------------------------------------------------
# 3. THE WINDOW SUM IS USED, NOT THE GAUSSIAN INTEGRAL
# --------------------------------------------------------------------------
def test_yield_is_window_sum_not_gaussian_integral(cfg):
    """Inject a non-Gaussian tail: the window sum must see it, the integral must not.

    The superseded script's docstring claimed the yield was
    ``amp * sigma * sqrt(2 pi) / binW``. It was not -- that path existed only in
    a demonstration plot. This test injects a component the Gaussian cannot
    describe and asserts the two estimators DISAGREE, pinning which one the code
    returns: the window sum, which contains the tail, not the analytic integral,
    which by construction cannot.

    The shoulder is placed at ``mu + 2.5 sigma``: outside the core fit range
    (``mu +- 1.5 sigma``), so the iterative core restriction of step 4 refits
    without it and the fitted Gaussian stays a clean description of the core --
    but inside the ``+-3 sigma`` window, so the window sum contains it. The two
    estimators are then separated by construction, and only one of them can be
    right about the yield.
    """
    axis = cfg.mgg
    mu, sigma = 0.135, 0.012
    n_core, tail_frac = 20000.0, 0.25
    n_same, n_mixed, sums, _ = make_spectra(
        axis, n_signal=n_core, mu=mu, sigma=sigma, mixed_per_bin=2000.0,
        alpha_true=0.05, tail_frac=tail_frac, tail_offset_mult=2.5, tail_sigma=0.002,
    )
    by = extract.extract_bin(0, n_same, n_mixed, sums, cfg)
    assert by.ok, by.reason
    # The core restriction must have kept the fit on the core, not on the shoulder.
    assert by.fit.mu == pytest.approx(mu, abs=0.001)
    assert by.fit.sigma == pytest.approx(sigma, abs=0.001)

    integral = extract.gaussian_integral_yield(by.fit, axis)

    # The truth: everything actually inside the +-3 sigma window.
    core = gauss_counts(axis, n_core, mu, sigma)
    tail = gauss_counts(axis, n_core * tail_frac, mu + 2.5 * sigma, 0.002)
    true_in_window = float(np.sum((core + tail)[by.window]))

    assert by.value == pytest.approx(true_in_window, rel=0.02), (
        "the window sum must recover everything in the window, shoulder included"
    )
    assert by.value > 1.10 * integral, (
        f"window sum {by.value:.0f} must EXCEED the Gaussian integral {integral:.0f}: the "
        f"tail is in the window but cannot be in the analytic integral of the core Gaussian. "
        f"If these agree, the code is using the integral -- the very bug this pins."
    )
    assert abs(by.value - integral) / by.value > 0.10


# --------------------------------------------------------------------------
# 4. THE ABSCISSA -- the central test of the rewrite
# --------------------------------------------------------------------------
def test_abscissa_recovers_signal_mean_not_full_window_mean(cfg):
    """Signal <z> and background <z> differ: the abscissa must return the SIGNAL's.

    This is the defect the rewrite exists to fix. The old analysis reported a
    geometric box centre; a naive fix reports the full-window mean, which the
    combinatorial background drags toward its own <z>. Only the sideband-
    subtracted, count-weighted mean returns where the pi0 actually are.
    """
    axis = cfg.mgg
    z_signal, z_bkg = 0.6, 0.3
    n_same, n_mixed, sums, _ = make_spectra(
        axis, n_signal=20000.0, mu=0.135, sigma=0.012, mixed_per_bin=2000.0,
        alpha_true=0.05, z_signal=z_signal, z_bkg=z_bkg,
    )
    by = extract.extract_bin(0, n_same, n_mixed, sums, cfg)
    assert by.ok, by.reason

    got = by.abscissae["z"]
    full_window_mean = float(np.sum(sums["z"][by.window])) / float(np.sum(n_same[by.window]))

    assert got == pytest.approx(z_signal, abs=0.005), (
        f"the extracted abscissa {got:.4f} must recover the SIGNAL's <z> = {z_signal}"
    )
    assert abs(full_window_mean - z_signal) > 0.05, (
        "the test is only meaningful if the naive full-window mean is genuinely wrong"
    )
    assert abs(got - z_signal) < abs(full_window_mean - z_signal) / 10.0

    # pT2 too, where the background sits ABOVE the signal rather than below --
    # the correction must work in both directions, not just pull one way.
    assert by.abscissae["pt2"] == pytest.approx(0.15, abs=0.005)
    # And an axis where signal and background agree must be untouched.
    assert by.abscissae["q2"] == pytest.approx(2.0, abs=1e-6)


def test_abscissa_never_returns_a_box_centre(cfg):
    """A bin whose abscissa is not computable yields nan, so the caller must drop it."""
    axis = cfg.mgg
    n_same = np.zeros(axis.n)
    sb = axis.mask(cfg.extraction.sideband_min, cfg.extraction.sideband_max)
    win = axis.mask(0.11, 0.16)
    val = extract.abscissa(np.zeros(axis.n), n_same, np.ones(axis.n), 0.1, win, sb)
    assert not np.isfinite(val), "an uncomputable abscissa must be nan, never a fallback centre"


def test_abscissa_out_of_box_bin_is_dropped(cfg):
    """An abscissa outside the bin's own box is a broken estimate; drop the bin.

    A sample confined to a box has its mean inside that box. A value outside is
    therefore impossible, and a bin reported at an impossible place is worse
    than a missing bin.
    """
    axis = cfg.mgg
    n_same, n_mixed, sums, _ = make_spectra(
        axis, n_signal=20000.0, mu=0.135, sigma=0.012, mixed_per_bin=2000.0,
        alpha_true=0.05, z_signal=0.6,
    )
    box = {"q2": (1.0, 11.0), "xb": (0.1, 0.7), "z": (0.0, 0.2), "pt2": (0.0, 1.5)}
    by = extract.extract_bin(0, n_same, n_mixed, sums, cfg, box=box)
    assert not by.ok
    assert by.reason == extract.ABSCISSA_OUT_OF_BOX


# --------------------------------------------------------------------------
# 5. the quality gates reject what they should
# --------------------------------------------------------------------------
def test_gate_min_stats_rejects_thin_bin(cfg):
    """Fewer than min_stats raw same-event entries in the fit range -> rejected."""
    axis = cfg.mgg
    n_same = np.zeros(axis.n)
    fit_mask = axis.mask(cfg.extraction.fit_range_min, cfg.extraction.fit_range_max)
    n_same[np.flatnonzero(fit_mask)[:1]] = cfg.extraction.min_stats - 1
    n_mixed = np.full(axis.n, 100.0)
    sums = {k: np.zeros(axis.n) for k in ("q2", "xb", "z", "pt2")}
    by = extract.extract_bin(0, n_same, n_mixed, sums, cfg)
    assert not by.ok
    assert by.reason == extract.LOW_STATS


def test_gate_min_peak_amp_rejects_pure_background(cfg):
    """Pure background with no peak: the subtracted amplitude gate must reject it."""
    axis = cfg.mgg
    n_mixed = np.full(axis.n, 1000.0)
    n_same = 0.05 * n_mixed  # exactly the background, no signal at all
    sums = {k: np.zeros(axis.n) for k in ("q2", "xb", "z", "pt2")}
    by = extract.extract_bin(0, n_same, n_mixed, sums, cfg)
    assert not by.ok
    assert by.reason == extract.LOW_PEAK


def test_gate_min_stats_boundary_is_inclusive(cfg):
    """A bin exactly at min_stats passes the stats gate (it fails later, if at all)."""
    axis = cfg.mgg
    n_same = np.zeros(axis.n)
    n_same[np.flatnonzero(axis.mask(cfg.extraction.fit_range_min,
                                    cfg.extraction.fit_range_max))[:1]] = cfg.extraction.min_stats
    n_mixed = np.full(axis.n, 100.0)
    sums = {k: np.zeros(axis.n) for k in ("q2", "xb", "z", "pt2")}
    by = extract.extract_bin(0, n_same, n_mixed, sums, cfg)
    assert by.reason != extract.LOW_STATS


def test_no_sideband_mixed_is_rejected_not_defaulted(cfg):
    """No mixed statistics in the sideband -> no alpha -> the bin is dropped."""
    axis = cfg.mgg
    n_same = gauss_counts(axis, 5000.0, 0.135, 0.012)
    n_mixed = np.zeros(axis.n)
    sums = {k: np.zeros(axis.n) for k in ("q2", "xb", "z", "pt2")}
    by = extract.extract_bin(0, n_same, n_mixed, sums, cfg)
    assert not by.ok
    assert by.reason == extract.NO_SIDEBAND_MIXED


# --------------------------------------------------------------------------
# 6. bsa: per-run beam polarization (v0.2.1), never a default
# --------------------------------------------------------------------------
def test_polarization_table_is_read(cfg):
    """cuts.json carries the 9 Moller periods of RG-D, contiguous from 18305 to 19131."""
    pp = cfg.bsa.polarization_periods
    assert len(pp) == 9
    assert (pp[0].run_min, pp[0].run_max, pp[0].p, pp[0].p_err, pp[0].hwp) == (18305, 18401, 0.79228, 0.01498, "IN")
    assert (pp[-1].run_min, pp[-1].run_max, pp[-1].p, pp[-1].p_err) == (19021, 19131, 0.80262, 0.01463)
    assert cfg.bsa.polarization_period(18401) == 0
    assert cfg.bsa.polarization_period(18402) == 1
    assert cfg.bsa.polarization_period(18957) == 7
    assert cfg.bsa.polarization_period(18304) == -1
    assert cfg.bsa.polarization_period(19132) == -1


def test_no_override_means_use_the_table():
    """Neither --polarization nor --polarization-err: the per-run table is used."""
    assert bsa._resolve_polarization(None, None) is None


def test_override_needs_both_values():
    """P without sigma_P, or sigma_P without P, is refused."""
    with pytest.raises(bsa.PolarizationError, match="without its uncertainty"):
        bsa._resolve_polarization(0.86, None)
    with pytest.raises(bsa.PolarizationError, match="without --polarization"):
        bsa._resolve_polarization(None, 0.03)


def test_override_is_accepted():
    """An explicit (P, sigma_P) is accepted and returned unchanged."""
    assert bsa._resolve_polarization(0.86, 0.03) == (0.86, 0.03)


@pytest.mark.parametrize("bad", [0.0, -0.5, 1.5, 79.0])
def test_override_rejects_unphysical_polarization(bad):
    """P outside (0, 1] is not a polarization (79.0: a percentage)."""
    with pytest.raises(bsa.PolarizationError):
        bsa._resolve_polarization(bad, 0.03)


def _fake_data(cfg, bsa_pol):
    from types import SimpleNamespace
    return SimpleNamespace(path="fake.root", bsa_pol=bsa_pol)


def test_bin_polarization_by_hand(cfg):
    """<P>_bin = sum_k N_k |P_k| / sum_k N_k, the same weighting for sigma_P.

    Bin 0: 3 entries in period 0 (P = 0.79228 +- 0.01498) and 1 in period 7
    (P = 0.73788 +- 0.01463):
        <P>     = (3 * 0.79228 + 0.73788) / 4 = 0.77868
        sigma_P = (3 * 0.01498 + 0.01463) / 4 = 0.0148925
    Bin 1: only period 8, so exactly that period's P.
    Bin 2: no entries, so nan (never a default).
    """
    n4d, nper = cfg.binning.n_4d, len(cfg.bsa.polarization_periods)
    w = np.zeros((n4d, nper), dtype=np.int64)
    w[0, 0], w[0, 7] = 3, 1
    w[1, 8] = 5
    p, e = bsa.bin_polarization(_fake_data(cfg, w), cfg)
    assert p[0] == pytest.approx(0.77868, abs=1e-12)
    assert e[0] == pytest.approx(0.0148925, abs=1e-12)
    assert p[1] == pytest.approx(0.80262, abs=1e-12)
    assert e[1] == pytest.approx(0.01463, abs=1e-12)
    assert np.isnan(p[2]) and np.isnan(e[2])


def test_bin_polarization_override(cfg):
    """With an override every bin gets that value, even without a bsa_pol tree."""
    p, e = bsa.bin_polarization(_fake_data(cfg, None), cfg, (0.79, 0.015))
    assert np.all(p == 0.79) and np.all(e == 0.015)


def test_bin_polarization_refuses_an_old_file(cfg):
    """A Stage B file from before v0.2.1 has no bsa_pol tree: refuse, do not guess."""
    with pytest.raises(bsa.PolarizationError, match="rerun Stage B v0.2.1|Rerun Stage B v0.2.1"):
        bsa.bin_polarization(_fake_data(cfg, None), cfg)


def test_bin_polarization_refuses_a_different_table(cfg):
    """bsa_pol with a different number of periods than cuts.json is refused."""
    w = np.zeros((cfg.binning.n_4d, 3), dtype=np.int64)
    with pytest.raises(bsa.PolarizationError, match="not the one Stage B used"):
        bsa.bin_polarization(_fake_data(cfg, w), cfg)


# --------------------------------------------------------------------------
# 7. config and index formulas
# --------------------------------------------------------------------------
def test_missing_config_key_fails_loudly():
    """A missing key raises; it is never defaulted."""
    with pytest.raises(ConfigError, match="required config key"):
        require({"a": {"b": 1}}, "a/c")


def test_null_config_key_fails_loudly():
    """An explicit null is not a value either, for keys that require one."""
    with pytest.raises(ConfigError, match="present but null"):
        require({"a": {"b": None}}, "a/b")


def test_4d_index_round_trip(cfg):
    """encode_4d and decode_4d must invert each other over the whole grid."""
    b = cfg.binning
    for i in range(b.n_4d):
        assert b.encode_4d(*b.decode_4d(i)) == i


def test_3d_index_comes_from_the_decode_not_a_division(cfg):
    """bin3d = cell_a * n_z + i_z, always computed by decoding the 4D index.

    Stage B's provenance warns "NOT bin4d / n_pt2". Worth being precise about
    what that warning is and is not:

    * As INTEGER arithmetic the two coincide identically, because
      ``n_b = n_z * n_pt2`` by construction, so ``bin4d // n_pt2 ==
      cell_a * n_z + i_z == bin3d`` for every bin. The test below asserts that
      identity rather than pretending otherwise.
    * In *Python* the warning still bites hard: ``bin4d / n_pt2`` is TRUE
      division and returns a float (``7.4``), which is not an index at all and
      would raise or silently mis-index depending on where it lands.

    So the identity is not a licence to divide: the decode is the contract, the
    division is a coincidence of the shipped grid shape, and a future Grid B
    that is not a plain ``n_z x n_pt2`` product would break it while the decode
    would keep working.
    """
    b = cfg.binning
    for i in (0, 37, 501, b.n_4d - 1):
        i_q2, i_xb, i_z, _ = b.decode_4d(i)
        assert b.bin3d_of_4d(i) == b.encode_3d(i_q2, i_xb, i_z)
    # The integer identity, asserted rather than assumed.
    assert all(b.bin3d_of_4d(i) == i // b.n_pt2 for i in range(b.n_4d))
    # But true division is a float and is never an index.
    assert isinstance(37 / b.n_pt2, float)


def test_all_b_cells_of_an_a_cell_share_one_cell_a(cfg):
    """All Grid B cells of one Grid A cell map to the same N_DIS cell."""
    b = cfg.binning
    for cell in (0, 17, b.n_a - 1):
        idx = [cell * b.n_b + j for j in range(b.n_b)]
        assert {int(b.cell_a_of_4d(i)) for i in idx} == {cell}


def test_mgg_bin_centre_formula(cfg):
    """mgg centre = (imgg + 0.5) * width, matching Stage B's axis exactly."""
    axis = cfg.mgg
    assert axis.width == pytest.approx(0.3 / 200)
    assert axis.centres[0] == pytest.approx(0.5 * 0.3 / 200)
    assert axis.centres[-1] == pytest.approx(199.5 * 0.3 / 200)
