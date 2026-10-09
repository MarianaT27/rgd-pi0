#!/bin/bash
# =====================================================================================
# run_full_chain.sh : stages B, C and D of the rgd-pi0 chain on the FULL production
#
# Chain stages (numbering as in Table 5 of the note):
#   A. Skim      stageA_skim on the farm (swif2, one job per run). NOT run here.
#                Input for this script: 95 slims in farm_full/<target>/<run>/<job>/slim_0.root
#                (run 18540 was skimmed with --max-events, so its slim is marked truncated)
#   B. Grid      make_grid: equal-statistics binning from all 95 slims
#                -> config/binning/grid_A_q2_xb.json  (Q2, xB)
#                -> config/binning/grid_B_z_pt2.json  (z, pT2)
#   C. Bin       stageB_bin (the executable keeps its old name): fills the 4D/3D bins per target
#                -> analysis/full/binned_full/<target>.root
#   D. Extract   python modules in python/pi0, one block per observable:
#                D.1 ratio       R_A per 4D bin, CxC/Cu/Sn vs LD2  -> analysis/full/results_full/ratio_<T>.csv
#                D.2 broadening  Delta<pT2> per 3D bin, CxC/Cu/Sn  -> analysis/full/results_full/broadening_<T>.csv
#                D.3 bsa         A_LU per 4D bin, LD2/CxC/Cu/Sn    -> analysis/full/results_full/bsa_<T>.csv
#                D.4 qa          QA figures from LD2               -> analysis/full/figs_qa_full/
#                D.5 plots       note-style result figures         -> analysis/full/figs_results_full/
#
# Run (from anywhere, it moves to the repo root itself):
#   nohup bash analysis/scripts/run_full_chain.sh >& analysis/logs/run_full_chain.log &
#
# Stops at the first failing step (set -e). Every step is timestamped in the log.
# Re-running overwrites the grid files and everything in analysis/full/binned_full,
# results_full, figs_qa_full and figs_results_full.
# =====================================================================================
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0

FULL=/volatile/clas12/mtenorio/pi0/rgd-pi0/farm_full    # stage A output (input to B and C)
OUT=analysis/full                                       # everything this script writes, except the grid
GA=config/binning/grid_A_q2_xb.json
GB=config/binning/grid_B_z_pt2.json
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

# -------------------------------------------------------------------------------------
# B. Grid
# -------------------------------------------------------------------------------------
stamp "B. Grid: make_grid"
args=()
for f in $FULL/*/*/*/slim_0.root; do args+=(--input "$f"); done
n=$(( ${#args[@]} / 2 ))
echo "input files: $n"
if [ "$n" -ne 95 ]; then echo "expected 95 files, found $n -- stopping"; exit 1; fi
./build/src/tools/make_grid/make_grid "${args[@]}" --config config/cuts.json --out-a $GA --out-b $GB --threads 4

# -------------------------------------------------------------------------------------
# C. Bin
# -------------------------------------------------------------------------------------
mkdir -p $OUT/binned_full
for t in LD2 CxC Cu Sn; do
  args=()
  for f in $FULL/$t/*/*/slim_0.root; do args+=(--input "$f"); done
  stamp "C. Bin: $t ($(( ${#args[@]} / 2 )) files)"
  ./build/src/stageB_bin/stageB_bin "${args[@]}" --output $OUT/binned_full/$t.root --config config/cuts.json --grid-a $GA --grid-b $GB --allow-truncated-inputs --threads 4
done

# -------------------------------------------------------------------------------------
# D. Extract
# -------------------------------------------------------------------------------------
export PYTHONPATH=python
mkdir -p $OUT/results_full

# D.1 ratio
for t in CxC Cu Sn; do
  stamp "D.1 Extract, ratio: $t"
  python3 -m pi0.ratio --target $OUT/binned_full/$t.root --ld2 $OUT/binned_full/LD2.root --config config --out $OUT/results_full/ratio_$t.csv --allow-unpublishable
done

# D.2 broadening
for t in CxC Cu Sn; do
  stamp "D.2 Extract, broadening: $t"
  python3 -m pi0.broadening --target $OUT/binned_full/$t.root --ld2 $OUT/binned_full/LD2.root --config config --out $OUT/results_full/broadening_$t.csv --allow-unpublishable
done

# D.3 bsa
for t in LD2 CxC Cu Sn; do
  stamp "D.3 Extract, bsa: $t"
  python3 -m pi0.bsa --file $OUT/binned_full/$t.root --config config --out $OUT/results_full/bsa_$t.csv --allow-unpublishable
done

# D.4 qa
stamp "D.4 Extract, qa: LD2"
python3 -m pi0.qa --file $OUT/binned_full/LD2.root --config config --outdir $OUT/figs_qa_full --allow-unpublishable

# D.5 plots
stamp "D.5 Extract, plots_results"
python3 -m pi0.plots_results --results $OUT/results_full --outdir $OUT/figs_results_full

stamp "ALL DONE"