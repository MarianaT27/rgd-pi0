#!/bin/bash
# =====================================================================================
# run_wagon_compare.sh : R_A with the LD2 from the SIDIS train vs the LD2 from the RGDSidisWagon test
#
# Purpose: check the wagon at the level of the published quantity. The wagon file of run 18431 is
# incomplete (70.8% of the run) and prescaled (1 in 5 below Q2 = 2.0), but both scale numerator and
# denominator of Y_D/N_DIS_D equally inside any Q2 bin that lies entirely on one side of 2.0. The
# test grid analysis/cmp_wagon/grid_A_q2_xb_q2cut2.json (make_grid_q2cut2.py) puts an edge at 2.0,
# so every bin should give the same R_A with either LD2, within statistics.
#
# Inputs (v0.2.1 Stage A slims, nothing re-skimmed):
#   LD2 train : /work/clas12/mtenorio/pi0/stageA_v021/LD2/018431/*/slim_*.root
#   LD2 wagon : /work/clas12/mtenorio/pi0/wagon_test/slim_wagon_018431.root
#   CxC 18465, Cu 18564, Sn 18564 : /work/clas12/mtenorio/pi0/stageA_v021/<T>/<run>/*/slim_*.root
#
# Steps:
#   C. Bin      stageB_bin for the 5 inputs, test grid A, committed grid B
#               -> analysis/cmp_wagon/binned/{LD2_train,LD2_wagon,CxC,Cu,Sn}.root
#   D. Extract  pi0.ratio for CxC, Cu, Sn, once with each LD2
#               -> analysis/cmp_wagon/results/ratio_<T>_{train,wagon}.csv
#   E. Compare  compare_wagon_ra.py -> printout and analysis/cmp_wagon/ra_wagon_vs_train.png
#
# BEFORE running: python3 analysis/scripts/make_grid_q2cut2.py
# Run (from anywhere):
#   nohup bash analysis/scripts/run_wagon_compare.sh >& analysis/cmp_wagon/run_wagon_compare.log &
# =====================================================================================
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0

S=/work/clas12/mtenorio/pi0/stageA_v021
GA=analysis/cmp_wagon/grid_A_q2_xb_q2cut2.json
GB=config/binning/grid_B_z_pt2.json
OUT=analysis/cmp_wagon
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

if [ ! -f $GA ]; then
  echo "$GA missing -- run first: python3 analysis/scripts/make_grid_q2cut2.py"
  exit 1
fi
mkdir -p $OUT/binned $OUT/results
export PYTHONPATH=python

declare -A IN=(
  [LD2_train]="$S/LD2/018431/*/slim_*.root"
  [LD2_wagon]="/work/clas12/mtenorio/pi0/wagon_test/slim_wagon_018431.root"
  [CxC]="$S/CxC/018465/*/slim_*.root"
  [Cu]="$S/Cu/018564/*/slim_*.root"
  [Sn]="$S/Sn/018564/*/slim_*.root"
)

# -----------------------------------------------------------------------------------
# C. Bin
# -----------------------------------------------------------------------------------
for k in LD2_train LD2_wagon CxC Cu Sn; do
  args=()
  for f in ${IN[$k]}; do args+=(--input "$f"); done
  stamp "C. Bin: $k ($(( ${#args[@]} / 2 )) files)"
  ./build/src/stageB_bin/stageB_bin "${args[@]}" --output $OUT/binned/$k.root --config config/cuts.json --grid-a $GA --grid-b $GB --threads 4
done

# -----------------------------------------------------------------------------------
# D. Extract
# -----------------------------------------------------------------------------------
for t in CxC Cu Sn; do
  for d in train wagon; do
    stamp "D. Extract, ratio: $t with LD2_$d"
    python3 -m pi0.ratio --target $OUT/binned/$t.root --ld2 $OUT/binned/LD2_$d.root --config config --grid-a $GA --grid-b $GB --out $OUT/results/ratio_${t}_$d.csv --allow-unpublishable
  done
done

# -----------------------------------------------------------------------------------
# E. Compare
# -----------------------------------------------------------------------------------
stamp "E. Compare"
python3 analysis/scripts/compare_wagon_ra.py

stamp "ALL DONE"
