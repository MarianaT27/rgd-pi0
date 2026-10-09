# analysis/scripts

All commands below are run from the repo root:
/volatile/clas12/mtenorio/pi0/rgd-pi0

## Bookkeeping scripts

### check_runs.py
What it does: compares the run list in config/runs.json with the runs that
actually exist on /cache for each target (LD2, CxC, CuSn).
Reads: config/runs.json (relative path, so it must be run from the repo root)
       /cache/clas12/rg-d/production/pass1/recon/<target>/dst/recon
Writes: nothing, prints two lists per target:
        runs on disk but not in runs.json
        runs in runs.json but not on disk
Run:
    python3 analysis/scripts/check_runs.py

### beam_energy.py
What it does: queries RCDB for runs 18300 to 19140 and groups them by beam energy.
Reads: RCDB (mysql://rcdb@clasdb-farm.jlab.org/rcdb)
Writes: nothing, prints the grouping
Needs: the python rcdb package (clas12 environment)
Run:
    python3 analysis/scripts/beam_energy.py

### skim_runs.py
What it does: for every run in the SIDIS train skim on /cache
(/cache/clas12/rg-d/production/pass1/recon/<target>/dst/train/SIDIS/),
pulls the run conditions from RCDB and writes one row per run.
Reads: the SIDIS train directories on /cache, RCDB
Writes: skim_runs.csv in the CURRENT directory
Needs: the python rcdb package (clas12 environment)
Run (from the repo root, so the csv lands in analysis/runinfo):
    cd analysis/runinfo
    python3 ../scripts/skim_runs.py
    cd ../..

### findRuns.C
What it does: reads a HIPO file event by event and prints every event whose
run number differs from the expected one. Used to find the 9 events of run
18539 inside SIDIS_018540.hipo (entries 66,918,450 to 66,918,458).
Arguments: file name, expected run number (default 18540)
Writes: nothing, prints "entry N : run R, event E" for each foreign event
        and a total at the end
Run:
    clas12root -b -q 'analysis/scripts/findRuns.C("/path/to/file.hipo", 18540)'
