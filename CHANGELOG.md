# Changelog

One entry per patch, newest first. The version is the one in `meson.build` and
`python/pi0/__init__.py`, and the git tag of the commit that applies the patch.

## v0.2.1 (2026-10-29), patch `v0.2.1.patch`

Applies on top of v0.2.0. **Stage B must be rerun**: the BSA needs the new `bsa_pol` tree.

### Added: per-run beam polarization in the BSA
- **Rule** (LE note 11/14/24, following Stepan): each Moller measurement applies to the runs
  from its first run until the next measurement. No averaging over the run period.
- **`config/cuts.json`**: `bsa.polarization` is now a table `periods`, 9 rows
  (runs, |P|, σP, HWP, logbook entry), runs 18305 to 19131. The old `value`/`error` keys are removed.
  - |P| (the magnitude) is used, because `REC::Event.helicity` is already HWP-corrected.
  - The HWP column is information only; it is not used.
- **C++ loader** (`src/config/Cuts.hpp/.cpp`): reads the table and refuses an empty table,
  P not in (0, 1], negative σP, HWP not IN/OUT, overlapping periods or a gap between periods.
- **Stage B** (`src/stageB_bin/main.cxx`):
  - Refuses to start if any input run is not in the table, and stops if an event's run is not.
  - New tree `bsa_pol` (bin4d, period, counts): BSA entries of each 4D bin per period
    (both helicities, all m_γγ and φ_h bins of the bin).
  - New tree `runs` (run, period, P, events): events read per run.
  - The cutflow printout lists events and BSA entries per period; provenance records the table.
- **`pi0.bsa`**: each 4D bin is divided by its own beam polarization,
  ⟨P⟩_bin = Σ_k N_k |P_k| / Σ_k N_k (N_k = BSA entries of the bin in period k).
  - The bin's σP is the same weighted mean of the period σP. It is a fully-correlated scale
    uncertainty: it is included in `sA_LU` as before, and also written on its own.
  - The BSA CSV has two new columns, `P_mean` and `P_err`.
  - `--polarization P --polarization-err σP` is kept only as an override (one value for every bin,
    both must be given).
  - A Stage B file without `bsa_pol` (made before v0.2.1) is refused unless the override is given.
  - The weights N_k use all m_γγ bins of the 4D bin, not only the fitted mass window.
- **`pi0.config`**: `BsaCuts.polarization_periods` (same checks as C++), replaces
  `polarization` / `polarization_err`. `pi0.io` reads `bsa_pol` and `runs` into `StageBData`.
- **`analysis/scripts/polarization_check.py`**: diagnostic plots, in the standard target colours:
  (a) events per run vs run number, coloured by period, with |P| of the period;
  (b) histogram of the per-bin ⟨P⟩; (c) A_LU(per-run P) / A_LU(P = 0.790) vs 0.790/⟨P⟩_bin,
  which must lie on the diagonal. Needs a second `pi0.bsa` run with
  `--polarization 0.790 --polarization-err 0.015`.

### Changed
- **`python/pi0/plots_results.py`**: the number of z bins is read from the result tables instead
  of the fixed `N_Z = 5`, so R_A vs z and Δ⟨pT²⟩ vs z show every Grid B z bin (9 with binning v2;
  before, only z < 0.5 was drawn).
- **Standard target colours**: D black, C `#1f77b4` (blue), Cu `#d62728` (red), Sn `#2ca02c` (green).
  - `python/pi0/plots_kinematics.py` (before: LD2 blue, C teal)
  - `analysis/scripts/clas6_compare.py` (before: Cu orange). CLAS6 points: C in the same blue
    (same nucleus), Fe `maroon` (next to Cu), Pb `darkgreen` (next to Sn), open markers.
- **`.gitignore`**: `/farm_ndis/`, `/patches/`, `/analysis/binning_v2/`, `/analysis/clas6_q2nu/`,
  `/analysis/ndis/`, `/*.tar.gz`.
- **Version 0.2.0 → 0.2.1** in `meson.build` and `python/pi0/__init__.py`.

### Tests
- `test_selection`: the polarization table is read (9 periods, values, run → period mapping, runs
  outside the table give -1); bad tables are refused (overlap, gap, negative P, P as a percentage,
  empty table).
- `test_extract.py`: the 3 BSA tests that failed since cuts.json got a polarization value are
  replaced by tests of the new behaviour: table read, override needs both values, out-of-range P
  refused, ⟨P⟩ and σP of a bin computed by hand, override, refusal of an old file and of a table
  that does not match the file.
- All C++ unit tests pass. Python: 124 pass, 0 fail.

### Not included (later)
- Separate A_LU for HWP IN and HWP OUT runs.

---

## v0.2.0 (2026-10-28), patch `v0.2.0.patch`

Applies on top of `binning_clas6.patch`.

### Changed
- **Beam energy 10.53 → 10.54 GeV** (RG-D value).
  - `config/cuts.json` (`beam.energy_gev`)
  - `src/core/Constants.hpp` (`kBeamEnergyGeV`). Cuts::load refuses a `cuts.json` whose
    beam energy differs from this constant, so both change together.
  - `python/pi0/plots_kinematics.py` (`EBEAM`, fallback only, used when `cuts.json` has no value)
  - Comments that quote 10.53 as what the old analysis did are unchanged.
- **Plots without connecting lines**: `python/pi0/plots_results.py`, R_A vs z and
  Δ⟨pT²⟩ vs z (`ls="none"`). The A^n fit curve in Δ⟨pT²⟩ vs A is kept.
- **`python/pi0/extract.py`**: the check that a bin's mean value lies inside its own box
  now uses the Grid A axis name (`xb` or `nu`). Without it, every bin of a (Q², ν) grid
  was dropped, because ⟨x_B⟩ was compared with the ν edges.
- **Version 0.1.0 → 0.2.0** in `meson.build` and `python/pi0/__init__.py`.

### Added
- `analysis/scripts/run_binning_v2_chain.sh`: Stages B–D on the `farm_ndis` slims, two passes:
  (Q², x_B) × Grid B into `analysis/binning_v2/`, and (Q², ν) × Grid B into
  `analysis/clas6_q2nu/` followed by the comparison. Does not run `make_grid`; refuses to
  start if `config/binning` has uncommitted changes.
- `analysis/scripts/clas6_compare.py`: RG-D π0 R_A in the bins of Mineeva et al. (CLAS EG2 π0,
  arXiv:2406.04513), (ν, z, pT²) for 1 < Q² < 4.1 and (Q², ν, z) for pT² < 1.5. For each of
  their bins the m_γγ spectra of the RG-D 4D bins inside it are summed and extracted; N_DIS
  is summed over the same cells. Writes a CSV and two plots.
- `analysis/scripts/mineeva2024_tables.csv`: their Tables I and II (C, Fe, Pb; value, stat,
  sys). The 9th block of Table I is labelled "2.2 < ν < 3.2, 0.5 < z < 0.6" in the PDF; by
  its position it is 3.73 < ν < 4.25 and is stored that way.
- `CHANGELOG.md` (this file).

### Tests
- `test_selection`: expected values now match `cuts.json`: beam 10.54, electron p > 1.8 GeV,
  DC edges 1.88 / 2.08 / 8.62 cm, cutflow label "1.8".
- `test_donorpool`: pool grid 8 × 8 × 4 = 256 bins (x_B edge at 0.06 from commit 905cef6).
- `test_core`: the hand-computed DIS point recomputed at 10.54 GeV.
- All C++ unit tests pass. Python tests: 116 pass; the same 3 BSA tests in
  `test_extract.py` that failed before this patch still fail (not touched here).

### Notes
- Running Stages B–D on the current `farm_ndis` slims (Stage A at 10.53) mixes two beam
  energies: Q², ν and x_B come from the slims (10.53), while pT², φ_h and x_F are computed
  in Stage B with 10.54. The difference is 0.01 GeV. Stage B also warns that the slims were
  skimmed against a different `cuts.json`. Both go away after the Stage A rerun.

---

## Before v0.2.0 (applied on top of 0.1.0, commit 905cef6; not tagged)

### `binning_clas6.patch` (2026-10-28)
- Grid B (`config/binning/grid_B_z_pt2.json`) set by hand: z = 0, 0.1, 0.2, 0.3, 0.4, 0.5,
  0.6, 0.7, 0.8, 1.0 and pT² = 0, 0.1, 0.25, 0.4, 0.55, 0.75, 1.5, 10 GeV². Inside
  0.3 < z < 1.0 and pT² < 1.5 these are the edges of Mineeva et al.; the rest are extra bins.
- New `config/binning/grid_A_q2_nu.json` for the comparison: Q² = 1.0, 1.33, 1.76, 4.1, 11 GeV²
  and ν = 2.2, 3.2, 3.73, 4.25, 5.25, 6.25, 7.25, 9.0 GeV.
- Grid A's second axis may be `xb` or `nu` (C++ loader and Stage B; Python `load_config`).
- `--grid-a` / `--grid-b` options in `pi0.ratio`, `pi0.broadening`, `pi0.bsa`, `pi0.qa`.
- Grids with `provenance.source` starting with "manual" are not treated as placeholders.
- `cuts.json` grid B shape 9 × 7; `test_binning` updated to the shipped grids.

### `stageA_ndis.patch` (2026-10-21)
- Stage A keeps every DIS event, including those with no selected photon, so N_DIS is
  inclusive. Before, only events with at least one photon were written, which biased R_A
  upward by p_D/p_A (1.10 C, 1.29 Cu, 1.43 Sn).
- Provenance `photon.min_multiplicity` and the cutflow row label updated accordingly.

### `xf_range.patch` (2026-10-20)
- x_F axis of the Stage B x_F trees widened to −3 < x_F < 1.5 (90 bins).

### `xf_study.patch` (2026-10-20)
- x_F added to `SidisKin` (`src/core/Types.hpp`) and computed in `compute_sidis`
  (`src/core/Kinematics.cpp`), x_F = p_L*/p_L,max* with the proton mass; NaN below
  W = M_p + m_π.
- Stage B writes trees `xf_z` (50 fine z bins) and `xf_zb` (Grid B z bins) with same- and
  mixed-event counts per (z, x_F, m_γγ).
- Two x_F unit tests in `test_core`.
