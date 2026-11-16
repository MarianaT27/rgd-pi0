# Changelog

One entry per patch, newest first. The version is the one in `meson.build` and
`python/pi0/__init__.py`, and the git tag of the commit that applies the patch.

## v1.0.0 (2026-11-16), patch `v1.0.0.patch`

Applies on top of v0.2.8. Inbending is now supported: every polarity-dependent cut has an
outbending and an inbending value, and the inbending values are those of the RG-D common aspects
analysis note (v3.1). The outbending values do not change. `cuts.json` changes shape, so its
sha256 changes: slims and binned files made with v0.2.8 or earlier show a config-hash mismatch
when run through v1.0.0 Stage B or the Python extraction.

### Changed: `config/cuts.json` holds one block per torus polarity
- `/beam/polarity` is gone. `/electron`, `/photon`, `/pairing` and `/vertex` move to
  `/polarities/outbending/...` (values unchanged, comments updated) and `/polarities/inbending/...`.
  Everything else (`beam`, `dis`, `binning`, `mgg_histogram`, `bsa`, `broadening`, `extraction`,
  `mixing`) is shared.
- Inbending values, common aspects note v3.1:
  - DC fiducial edges, Table 5 (IB, regional averages): R1 2.58, R2 2.51, R3 1.37 cm (OB: 1.88, 2.08, 8.62).
  - v_z windows, Table 6, mu -+ 3 sigma after the v_z correction: Cu -7.765 / 0.321, Sn -2.812 / 0.315
    (CuSn IB); CxC -7.670 / 0.313 and -2.713 / 0.307 (CxC IB). LD2 raw [-15, 5] cm, as outbending.
  - Sampling fraction: Table 23, the Average coefficients (`SamplingFraction.cpp kInbending`; the old
    rows were the legacy RG-A table and did not match the note). Outbending keeps Table 22, Average.
  - Everything else as outbending (Table 7): p > 1.8 GeV, chi2pid within +-5, status in
    (-4000, -2000], SF 3.5 sigma, PCAL lv, lw > 9 cm. Photon and pair cuts: the same values as
    outbending, in their own block so they can change apart.
- Comments: `mixing._comment` describes the frozen donor pool (the mix-then-insert text was stale);
  the DIS comment says 10.54 GeV (was 10.53).

### Changed: the code takes the polarity
- `Cuts::load(path, polarity)` reads only `/polarities/<polarity>`; an unknown polarity or a missing
  key is refused. The old refusal of inbending is removed.
- **`stageA_skim --polarity <outbending|inbending>` is required.** It is written to the slim
  provenance as `polarity` and picks the sampling-fraction table and the v_z corrector variant.
  `pi0.batch` and `pi0.local_batch` pass the farm config's polarity. Stage A does not check the run
  against `runs.json`; the run-list filter of `pi0.batch` does, as before.
- `stageB_bin` and `make_grid` take the polarity from their inputs' Stage A provenance and refuse
  inputs that disagree or carry none (`chain_polarity` in `InputChain.hpp`). The grid provenance
  records `polarity`.
- Python: `pi0.config.cuts_for_polarity(cuts, polarity)` lifts one block to the top level; `batch`,
  `io` (RG-A fallback check), `plots_selection`, `plots_kinematics` and
  `analysis/scripts/compare_pt_methods.py` use it with the polarity of the slims or farm config.
- The stale comment "This is 2.0 GeV" in `ElectronSelection.cpp` now says p > 1.8 GeV from
  `cuts.json` (note Table 7); `Cuts.hpp` and `Kinematics.hpp` comments updated likewise.

### Added: inbending farm configs, `only_runs`
- `config/farm/v100_IB/{LD2,CxC,Cu,Sn}.farm.json` (polarity inbending, output
  `farm_v100_IB/<T>`) and `Cu_specialDC`, `Sn_specialDC` for the special-DC CuSn runs 18358 to
  18366 and 18368: the main Cu and Sn configs exclude them, the special-DC ones keep only them,
  so they land in their own folder. Accepted runs: LD2 24, CxC 13, Cu 27, Sn 27, special DC 10.
- New optional farm key `only_runs`: keep only these runs of the list (an overlap with
  `exclude_runs` is refused; it never adds a run of the other polarity).
- `analysis/scripts/run_v100_IB_chain.sh`: the inbending Stage B, extraction and plots, written to
  `analysis/v100_IB/`, apart from outbending. Default: the shipped Grid A and Grid B (the outbending
  bins); `MAKE_IB_GRID=1` builds an inbending Grid A instead. It never writes `config/binning/`.

### Changed: run 18753 skips its foreign events (Cu, Sn)
- In the v022 part2 production, Cu and Sn job 105 failed with Stage A exit 5: the SIDIS-train
  file of run 18753 holds 2 events of run 18752. `config/farm/v022/{Cu,Sn}.farm.json`:
  `skip_foreign_run_runs` is `[18753, 19102]`.

### Fixed: figures and text
- The Preliminary label sits under the figure, not on the x-axis title, the points or the legend:
  `pi0.plots_results` (all three figures), `analysis/scripts/ra.py` (R_A vs z, R_A vs z in nu panels).
- `pi0.qa`: the Preliminary label and the provenance line were on one baseline and ran into each
  other; now two lines under the figure.
- `analysis/scripts/clas6_compare.py`: the titles cite Mineeva et al., PRC 112, 035203 (was the arXiv number).
- `analysis/scripts/hermes_ra.py` docstring and the v0.2.6 entry below: RG-D reaches nu < 8.74 GeV
  (E - E' with p_e > 1.8 GeV), not 8.96 (that is y < 0.85 alone).

### Tests
- C++: `test_selection` reads both blocks, refuses an unknown polarity, reads only the requested
  block, and pins all six inbending SF rows; `test_inputchain` covers `chain_polarity`;
  `test_donorpool` loads with outbending.
- Python: new `test_cuts_polarity.py` (outbending equals v0.2.8 value by value,
  `python/tests/data/cuts_v028_outbending_values.json`; inbending equals the note); `test_farm.py`
  covers `--polarity` in the job script, `only_runs`, the v100_IB configs and the v022 skip lists.

### Version
- Version 1.0.0.

## v0.2.8 (2026-11-16), patch `v0.2.8.patch`

Applies on top of v0.2.7. Farm configs only: no change to the code, the cuts or the run lists.

### Changed: run 19102 is skimmed with `--skip-foreign-run --run 19102`
- In the v022 part1 production (2026-11-15) the Cu and Sn job 66 failed with Stage A exit 5:
  the SIDIS-train file `sidis_019102.hipo` holds 2 events of run 19101 next to 11721840 events
  of run 19102. Same issue as run 18540 (v0.2.2, v0.2.3).
- `config/farm/v022/Cu.farm.json` and `Sn.farm.json`: `skip_foreign_run_runs` goes from `[]` to
  `[19102]`, and its comment says why. Cu and Sn read the same CuSn files, so both need it.
- Stage A then drops the 2 events of run 19101 and keeps every event of run 19102. No other run
  changes. The 207 jobs that succeeded each read a single-run file, so they stay as they are.

### Tests (`test_farm.py`)
- New `test_the_v022_configs_list_the_known_mixed_runs`: the four v022 configs list exactly
  LD2 18540, Cu 19102, Sn 19102 (CxC none), and the Cu wrapper adds the option for 19102 only.
- `test_v022_farm_configs_load_and_read_the_new_sidis_train`: expects 19102 for Cu and Sn (was empty).

### Version
- Version 0.2.8.

## v0.2.7 (2026-11-15), patch `v0.2.7.patch`

Applies on top of v0.2.6. Run lists only: no change to the cuts, Stage A, Stage B or the extraction.
Existing slims and results stay as they are (every run used so far is in the new lists).

### Changed: `config/runs.json` now holds L. El Fassi's run lists
- The inbending and outbending lists for LD2, CxC and CuSn are replaced with L. El Fassi's pass1
  lists, read on 2026-11-15 from the IB and OB links of the SIDIS train
  (`/work/clas12/rg-d/production/pass1/recon/<T>/dst/train/sidis/<IB|OB>/`).
- Counts: outbending LD2 136, CxC 116, CuSn 290; inbending LD2 24, CxC 13, CuSn 37 (616 runs,
  542 outbending). Up to v0.2.6: 574 runs from clas-framework `clas12/Runs.hpp`.
- 42 runs added, none removed:
  - outbending LD2: 18432
  - outbending CxC: 18482, 18485, 18811, 18824
  - outbending CuSn: 18560, 18562, 18637, 18664, 18674, 18687, 18688, 18689, 18726, 18737, 18739,
    18874, 18887, 18897, 18899, 18901, 18917, 18931, 19061, 19075, 19124, 19128
  - inbending CxC: 18369, 18370, 18371, 18400, 18401
  - inbending CuSn: 18358 to 18366 and 18368 (the special-DC runs of `special_DC_runs_CuSn`).
    Whether to keep them is to be decided before the inbending production.
- Run 18432 (`LD2_trigger_rgd_v2_2_Q2_2_5`, a Q^2 > 2.5 trigger configuration) was excluded up to
  v0.2.6. Decided 2026-11-15: included as a normal outbending LD2 run, as in L. El Fassi's list.
- `special_DC_runs_CuSn` and `LD2_trigger_sets` are unchanged. `_comment` and `_generated` record the
  new source.
- The note's `tab:runs` no longer matches `runs.json`.

### Tests (`test_farm.py`)
- `test_runs_json_matches_the_note` becomes `test_runs_json_matches_lamiaas_run_lists` (new counts).
- New `test_no_run_is_in_two_lists`; outbending total 515 -> 542.
- `test_a_run_in_neither_list_is_rejected` now uses run 18418 (18432 is in a list now).
- New `test_18432_is_in_the_outbending_ld2_list`.
- 149 Python tests pass.

### Docs and version
- README: run-list section updated (source, counts, 18432).
- Version 0.2.7.

## v0.2.6 (2026-11-14), patch `v0.2.6.patch`

Applies on top of v0.2.5. Plots and labels only: no change to Stage A, Stage B, the extraction or any
number. Rerun only the plot scripts on the existing files (see README_v0.2.6.md).

### Fixed: the Preliminary label (bug of v0.2.3)
- Stage A stamps `stageA_skim.code_version` as the full meson string
  `0.2.1 (meson project version; source tree is not under VCS, so no commit is recorded)`.
  v0.2.3 split it at its commas, so labels came out garbled, and the long text widened the
  `results_*` figures (saved with `bbox_inches="tight"`) and pushed the plot to the right.
- `pi0.io`: only the x.y.z number is taken from the stamp (`stageA_skim.code_versions` now holds
  e.g. `0.2.1, 0.2.2`). The label names the OLDEST version among the inputs, the production's:
  `Preliminary (v0.2.1, RG-A photon ID)` for the v0.2.1 slims, of which at least one input carries
  a v0.2.2 stamp.
- New `pi0.io.preliminary_label_from_text`: rebuilds the label from a CSV header line.
  `pi0.plots_results` uses it, so the CSVs written by v0.2.3 to v0.2.5 work without being remade.
- Tests with the real meson string and the exact garbled header line (`test_preliminary.py`).

### Fixed: BSA label on top of the legend
- `bsa.py`: the label moves from the top left of the upper panel (where the 4-column legend is) to
  below the figure, bottom right.

### Changed: CLAS6 comparison (`clas6_compare.py`)
- CLAS6 C has the same blue as RG-D C (same nucleus, told apart by the marker: RG-D filled circles,
  CLAS6 open squares); CLAS6 Fe orange `#ff7f0e`, Pb gold `#d4b000`. RG-D colours unchanged.
- Legend font 7 -> 10, markers in the legend x1.4.
- New `clas6_fig2.pdf/.png`, as Mineeva et al. Fig. 2: R_A vs z, all 9 (Q^2, nu) bins in one panel,
  x 0.27-0.92, y 0-1.1 (read off their figure by eye), RG-D and CLAS6 together, the 9 bins spread by
  +-0.02 inside each z bin. x of a bin = RG-D LD2 mean z of that bin, for RG-D and CLAS6 (their table
  has no mean z). The two existing figures stay.

### Changed: HERMES comparison (`hermes_ra.py`)
- HERMES points added: new `analysis/scripts/hermes2007_pi0.csv`, the 16 pi0 tables of HEPData
  ins749249 v1 (doi:10.17182/hepdata.13387.v1, Tables 97 to 112): He, Ne, Kr, Xe vs nu, z, Q^2, p_T^2,
  copied unchanged. Open squares at their mean, stat. errors only; colours He `#9467bd`, Ne `#8c564b`,
  Kr `#e377c2`, Xe `#7f7f7f`.
- Layout as their pi0 figure: 4 panels in one row, shared y, no gaps. Ranges read off their figure by
  eye: y 0.4-1.2, nu 1-25 GeV (linear), z 0-1.1 (linear), Q^2 0.8-16 GeV^2 (log), p_T^2 0.01-4 GeV^2
  (log). RG-D points above 1.2 fall outside. Legend to the right, font 11.
- Tests: the HEPData copy against the HERMES grid edges and a spot value (`test_hermes.py`).

- Version 0.2.6 in `meson.build` and `python/pi0/__init__.py`.

## v0.2.5 (2026-11-13), patch `v0.2.5.patch`

Applies on top of v0.2.4. Test-only: no code, grid or config change, nothing to rerun.

### Fixed: `tests/test_binning.cpp` follows the v0.2.2 Grid A
- The binning test had failed since v0.2.2, which remade Grid A (Q^2 edge pinned at 2.0) without updating it.
- "the shipped grids have the geometry the note specifies": the Grid A Q^2 and x_B edges are now the ones in
  `config/binning/grid_A_q2_xb.json`.
- "an edge literal too large for a double": the Q^2 edge it overwrites with `1e999` is now `1.0883105529903487`
  (was `1.1036395052848946`, an edge of the v0.2.1 Grid A that no longer exists in the file).
- Comments quoting the bin boxes of Q^2 = 2.5 and x_B = 0.25 updated; the bin indices they check are unchanged.
- Version 0.2.5 in `meson.build` and `python/pi0/__init__.py`.

## v0.2.4 (2026-11-13), patch `v0.2.4.patch`

Applies on top of v0.2.3. Nothing in the production changes: without the new flag Stage B behaves
exactly as in v0.2.3. Stage B must be rebuilt (`meson compile -C build`) for the HERMES chain.

### Added: RG-D R_A in the HERMES bins, with the HERMES cuts (RG-D only)
- Reference: Airapetian et al. (HERMES), Nucl. Phys. B 780 (2007) 1, arXiv:0704.3270. Edges from their
  Table A.1, cuts from Sec. 5 (nu > 6 GeV, z > 0.2, x_F > 0). Each panel cuts on the other variables,
  not its own: R vs nu (z > 0.2), R vs z (nu > 6), R vs Q^2 and R vs p_T^2 (nu > 6, z > 0.2); x_F > 0 in all.
- The cuts apply to this plot only. Production grids, cuts.json and binned files are not touched.
- **Grids** `config/binning/grid_A_q2_nu_hermes.json` (Q^2 = 1, 1.5, 2, 3, 4, 5, 6, 8, 25 GeV^2;
  nu = 4, 6, 8, 10 GeV) and `config/binning/grid_B_z_pt2_hermes.json` (z = 0.1, ..., 0.9, 1.3;
  p_T^2 = 0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.4, 1.95, 5 GeV^2). nu > 6 and z > 0.2 are bin edges,
  applied by bin selection. RG-D reaches nu < 8.74 GeV (E - E' with p_e > 1.8 GeV; corrected in v1.0.0, was 8.96 from y < 0.85): only the first three nu bins.
- **`stageB_bin --xf-min X`**: keeps only pi0 pairs with x_F > X, same- and mixed-event, right after the
  z window. Default: no x_F cut. X must be in (-1, 1). A pair with x_F not computed (NaN) is dropped.
  - New provenance key `pairs.xf_min`, always written (`none (the default: no x_F cut)` without the flag).
  - New line `x_F cut (--xf-min)` in the cutflow printout.
  - `pi0.io.check_provenance`: a file binned with an x_F cut gets a WARNING (not fatal),
    `provenance/pairs.xf_min`. A file without the key (made before v0.2.4) has no x_F cut.
- **`analysis/scripts/hermes_ra.py`**: for each HERMES bin, sums the m_gg spectra of the RG-D 4D bins
  inside it and runs the production extraction on the sum (the method of `clas6_compare.py`); N_DIS summed
  over the same (Q^2, nu) cells. Writes `hermes_ra.csv` and `hermes_ra.pdf/.png` (2x2: R vs nu, z, Q^2,
  p_T^2), standard target colours, label `Preliminary (...)`, the x_F cut read from the files' provenance.
  Points at the target's own mean; for nu that is <Q^2>/(2M<x_B>) (Stage B does not sum nu).
- **`analysis/scripts/run_hermes_chain.sh`**: Stage B on the v0.2.1 slims with the HERMES grids and
  `--xf-min 0`, then `hermes_ra.py`. Output in `analysis/hermes/` (in `.gitignore`).
- `pi0.plots_kinematics --nu-grid config/binning/grid_A_q2_nu_hermes.json` draws the HERMES (nu, Q^2)
  grid in the third panel of `binning_grid_factorized.pdf`.
- Tests: `python/tests/test_hermes.py` (grid edges, the cuts as bin selections, the x_F warning).
- **Not compiled here**: the C++ change was written without a ROOT build. Build and run `meson test -C build`
  on ifarm before `run_hermes_chain.sh`.
- Version 0.2.4 in `meson.build` and `python/pi0/__init__.py`.

## v0.2.3 (2026-11-13), patch `v0.2.3.patch`

Applies on top of v0.2.2. **Stage B must be rerun** for anything binned with the v0.2.2 Grid B
(see "Grid B restored"). Stage A does not need to be rerun.

### Fixed: Grid B restored to the hand-chosen CLAS6 edges
- `config/binning/grid_B_z_pt2.json` is again the v0.2.1 file: z = 0, 0.1, ..., 0.8, 1.0 and
  p_T^2 = 0, 0.1, 0.25, 0.4, 0.55, 0.75, 1.5, 10 GeV^2 (the CLAS6 edges inside 0.3 < z < 1 and p_T^2 < 1.5).
  v0.2.2 had committed by mistake the equal-statistics Grid B that `make_grid` wrote in `run_v021_chain.sh`.
- `run_v021_chain.sh` and `run_full_chain.sh` no longer overwrite Grid B: `make_grid --out-b` now writes
  `$OUT/grid_B_z_pt2_equal_stat.json`, which nothing reads. Grid A is still made by `make_grid` and pinned at
  Q^2 = 2.0 by `pin_q2_edge.py`.
- `clas6_compare.py` needs the CLAS6 edges in Grid B and stops with the v0.2.2 Grid B.

### Changed: "DIAGNOSTIC" replaced by "Preliminary (<version>, RG-A photon ID)"
- New `pi0.io.preliminary_label(...)`: `Preliminary (v0.2.1, RG-A photon ID)`, from the Stage A version(s)
  of the inputs (`stageA_skim.code_version`, collected over every input into `stageA_skim.code_versions`)
  and `RG-A photon ID` if any input used the RG-A GBT fallback.
- Used in: the first header line of the `pi0.ratio`, `pi0.broadening` and `pi0.bsa` CSVs; the `pi0.qa`
  figures (small label in the bottom-right corner instead of the diagonal watermark); the
  `pi0.plots_results` figures (label read from the CSV header; no label if the inputs have no blockers);
  the `pi0.io` provenance printout and error message; the `--allow-unpublishable` help texts.
- Analysis scripts `ra.py`, `dpt2.py`, `bsa.py`, `bsa_grid.py`: `Preliminary (<sample label>, RG-A photon ID)`.
  New sample `v021` (`analysis/v021`, the output of `run_v021_chain.sh`) with label `v0.2.1`.
- Only the wording changed. The provenance checks, the blocker list and `--allow-unpublishable` are unchanged.

### Added: farm configs for the new SIDIS train, `config/farm/v022/`
- Inputs `/cache/clas12/rg-d/production/pass1/recon/<LD2|CxC|CuSn>/dst/train/sidis/` (one file per run,
  prescaled 1 in 5 for Q^2 < 2.0), output `/volatile/clas12/mtenorio/pi0/rgd-pi0/farm_v022/<T>`,
  workflows `mt_rgd_pi0_stageA_v022_<T>`. Resources as in v021.
- New key `farm.skip_foreign_run_runs` (array of runs, default empty). Those runs are skimmed with
  `stageA_skim --skip-foreign-run --run <run>`: events of any other run in the file are dropped, every event
  of the run itself is kept. Every other run gets the same command line as before.
  `config/farm/v022/LD2.farm.json` lists 18540 (its file holds events of run 18539).
- `pi0.farm`: `FarmConfig.skip_foreign_run_runs`; a run in both `exclude_runs` and `skip_foreign_run_runs`
  is refused; `stagea_wrapper(..., run=)`. `pi0.batch` passes the run of each chunk and prints the runs
  that get the flag. Tests in `python/tests/test_farm.py`.

### Added: (nu, Q^2) panel in `binning_grid_factorized.pdf`
- `pi0.plots_kinematics` draws a third panel, nu vs Q^2 with the edges of `--nu-grid`
  (default `config/binning/grid_A_q2_nu.json`, the CLAS6 comparison grid). Same events as Grid A.
  nu is not a fourth variable: it is fixed by (Q^2, x_B), nu = Q^2 / (2 M x_B).

### Changed: CLAS6 colours in `clas6_compare.py`
- CLAS6 C is now navy `#0b2f5a` (was the same blue as RG-D C). Fe maroon and Pb dark green as before.
  RG-D keeps the standard colours (C `#1f77b4`, Cu `#d62728`, Sn `#2ca02c`).

- Version 0.2.3 in `meson.build` and `python/pi0/__init__.py`.

## v0.2.2 (2026-11-05), patch `v0.2.2.patch`

Applies on top of v0.2.1. Nothing needs to be rerun: without the new flag Stage A behaves exactly as in v0.2.1.

### Added: `stageA_skim --skip-foreign-run`
- **Why:** the LD2 SIDIS-train file of run 18540 holds 10 events of run 18539 (event numbers 6309641 to 6309650).
  Stage A refuses a file with two runs (exit 5), so the v0.2.1 production has no slim for 18540.
- **What it does:** events whose `RUN::config.run` is non-zero and differs from `--run` are counted and dropped
  before the electron search. They reach neither the slim nor the QA trees, and they do not enter the single-run check.
- **Requires `--run <N>`:** the run to keep is stated on the command line, never taken from whichever event comes first.
  Without `--run` the program refuses to start.
- **Default off:** without the flag a second run still exits 5, as before. The farm production is unchanged.
- **Recorded:**
  - provenance key `events.foreign_run_skipped`, always present: `run: events` for each skipped run, or `none` with the reason;
  - a cutflow row `events from run N (--skip-foreign-run; skipped ...)` after the first row, only when the flag is given;
  - a `foreign run` line in the startup printout.
- `events.input_tag0` still counts every event read, skipped ones included.
- Version 0.2.2 in `meson.build` and `python/pi0/__init__.py`.

### Not listed when v0.2.2 was released (added in v0.2.3)
- **Grid A** (`config/binning/grid_A_q2_xb.json`) remade by `make_grid` from the 95 v0.2.1 slims
  (`/work/clas12/mtenorio/pi0/stageA_v021/`), equal statistics, then the Q^2 edge nearest 2.0 (2.050) moved
  to 2.0 by `pin_q2_edge.py`, so no Q^2 bin straddles the SIDIS-train prescale threshold.
  Edges: Q^2 = 1.0, 1.088, 1.194, 1.323, 1.489, 1.713, 2.0, 2.674, 10.595 GeV^2;
  x_B = 0.061, 0.100, 0.128, 0.158, 0.192, 0.232, 0.290, 0.768. Recorded in the file (`_q2_edge_pinned`).
- **Grid B** (`config/binning/grid_B_z_pt2.json`) was ALSO overwritten by `make_grid`, with equal-statistics
  edges (z = 0.046, ..., 1.0; p_T^2 = 0, 0.018, ..., 8.153). This was not intended; v0.2.3 restores the
  CLAS6 grid of v0.2.1.
- **Scripts** (`analysis/scripts/`):
  - `run_v021_chain.sh`: stages B to D on the v0.2.1 slims (make_grid, pin_q2_edge, stageB_bin,
    pi0.ratio/broadening/bsa/qa/plots_results), output in `analysis/v021/`.
  - `run_wagon_compare.sh`, `compare_wagon_ra.py`: R_A with the LD2 of the SIDIS train vs the LD2 of the
    RGDSidisWagon test (run 18431), on a test grid with a Q^2 edge at 2.0 (`make_grid_q2cut2.py`).
  - `compare_wagon_train.py`, `compare_wagon_files.py`, `wagon_blocks.py`, `wagon_gaps.py`: wagon vs train
    event by event (prescale rule, missing events per block of event numbers).
  - `make_grid_q2cut2.py`: test copy of Grid A with the edge 1.836 moved to 2.0.
  - `pin_q2_edge.py`: moves the Q^2 edge of a Grid A file nearest a value onto that value.
- **Farm configs** `config/farm/v021/<T>.farm.json`: the v0.2.1 Stage A production
  (`train/SIDIS/`, output `farm_v021`).
- **`.gitignore`**: `/analysis/cmp_wagon/` and `/analysis/v021/`.

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
