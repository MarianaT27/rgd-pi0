#!/bin/bash
# =====================================================================================
# run_dst_compare.sh : DST vs SIDIS-train skim, same 3 runs, same grid and cuts
#
# Purpose: check whether the SIDIS train keeps every DIS event. Both inputs are the
# COMPLETE runs, so if the skim is inclusive, N_DIS must agree cell by cell.
#   LD2 run 18431, CxC run 18465, CuSn run 18564 (gives both Cu and Sn)
#
# Inputs (Stage A already done, nothing is re-skimmed here):
#   dst  : farm_dst/<T>/<run>/<job>/slim_0.root   (DST recon files, all files of the run)
#   skim : farm_full/<T>/<run>/<job>/slim_0.root  (SIDIS train, the full-production slims)
#
# Steps (same stage labels as run_full_chain.sh):
#   B. Grid      NOT run: uses the committed 5x5 grid in config/binning, so every bin
#                matches analysis/full one to one
#   C. Bin       stageB_bin per target, for dst and for skim
#                -> analysis/cmp_dst/binned_full/<T>.root, analysis/cmp_skim/binned_full/<T>.root
#   D. Extract   D.1 ratio and D.2 broadening for CxC, Cu, Sn (no BSA, QA or note plots)
#                -> analysis/cmp_<dst|skim>/results_full/{ratio,broadening}_<T>.csv
#                D.5 compare_dst_skim.py: N_DIS and pi0-window pair counts, dst vs skim
#
# BEFORE running: config/cuts.json and config/binning must be the committed
# full-production versions (the script stops otherwise), and all four farm_dst
# workflows must be finished.
#
# Run (from anywhere):
#   nohup bash analysis/scripts/run_dst_compare.sh >& analysis/logs/run_dst_compare.log &
# =====================================================================================
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0

GA=config/binning/grid_A_q2_xb.json
GB=config/binning/grid_B_z_pt2.json
declare -A RUN=([LD2]=018431 [CxC]=018465 [Cu]=018564 [Sn]=018564)
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

if ! git diff --quiet -- config/cuts.json config/binning/; then
  echo "config/cuts.json or config/binning differs from the committed version -- restore it first:"
  echo "  git checkout -- config/cuts.json config/binning/"
  exit 1
fi

export PYTHONPATH=python
for S in dst skim; do
  if [ "$S" = dst ]; then SRC=farm_dst; else SRC=farm_full; fi
  OUT=analysis/cmp_$S
  mkdir -p $OUT/binned_full $OUT/results_full

  # -----------------------------------------------------------------------------------
  # C. Bin
  # -----------------------------------------------------------------------------------
  for t in LD2 CxC Cu Sn; do
    args=()
    for f in $SRC/$t/${RUN[$t]}/*/slim_0.root; do args+=(--input "$f"); done
    stamp "C. Bin ($S): $t, run ${RUN[$t]} ($(( ${#args[@]} / 2 )) files)"
    ./build/src/stageB_bin/stageB_bin "${args[@]}" --output $OUT/binned_full/$t.root --config config/cuts.json --grid-a $GA --grid-b $GB --allow-truncated-inputs --threads 4
  done

  # -----------------------------------------------------------------------------------
  # D. Extract
  # -----------------------------------------------------------------------------------
  for t in CxC Cu Sn; do
    stamp "D.1 Extract, ratio ($S): $t"
    python3 -m pi0.ratio --target $OUT/binned_full/$t.root --ld2 $OUT/binned_full/LD2.root --config config --out $OUT/results_full/ratio_$t.csv --allow-unpublishable
  done
  for t in CxC Cu Sn; do
    stamp "D.2 Extract, broadening ($S): $t"
    python3 -m pi0.broadening --target $OUT/binned_full/$t.root --ld2 $OUT/binned_full/LD2.root --config config --out $OUT/results_full/broadening_$t.csv --allow-unpublishable
  done
done

stamp "D.5 compare dst vs skim"
python3 analysis/scripts/compare_dst_skim.py

stamp "ALL DONE"
