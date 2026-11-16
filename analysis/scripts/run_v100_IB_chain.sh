#!/bin/bash
# =====================================================================================
# run_v100_IB_chain.sh : the INBENDING chain (v1.0.0), kept apart from outbending
#
# Input: the v1.0.0 inbending Stage A production, made with
#   python -m pi0.batch config/farm/v100_IB/<T>.farm.json --submit      (T = LD2, CxC, Cu, Sn)
# i.e. slims in farm_v100_IB/<T>/<run>/<job>/slim_*.root, every one with provenance polarity = inbending.
# Expected runs (config/runs.json, inbending): LD2 24, CxC 13, Cu 27, Sn 27. The 10 special-DC CuSn
# runs (Cu_specialDC, Sn_specialDC) are NOT used here: they are checked on their own first.
#
#   B. Grid      default: NONE made. The chain uses the shipped Grid A and Grid B
#                (config/binning/grid_A_q2_xb.json, grid_B_z_pt2.json), the bins of the outbending
#                results, so the two polarities can be compared bin by bin.
#                With MAKE_IB_GRID=1: make_grid on the IB slims -> $OUT/grid_A_q2_xb_IB.json
#                (equal statistics, Q2 edge pinned at 2.0 GeV2), used for C and D instead.
#                config/binning/ is never written by this script.
#   C. Bin       stageB_bin per target -> $OUT/binned/<T>.root. Stage B reads the polarity from the
#                slims and applies the /polarities/inbending block of config/cuts.json.
#   D. Extract   D.1 ratio, D.2 broadening (CxC/Cu/Sn vs LD2), D.3 bsa (all four), D.4 qa (LD2),
#                D.5 plots_results -> $OUT/results, $OUT/figs_qa, $OUT/figs_results
#
# Run (from anywhere, it moves to the repo root itself):
#   nohup bash analysis/scripts/run_v100_IB_chain.sh >& analysis/logs/run_v100_IB_chain.log &
#   MAKE_IB_GRID=1 nohup bash analysis/scripts/run_v100_IB_chain.sh >& analysis/logs/run_v100_IB_chain.log &
#
# One log for the whole run (the nohup line above). Stops at the first failing step (set -e);
# every step is timestamped. Re-running overwrites everything in $OUT.
# =====================================================================================
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0

IB=/volatile/clas12/mtenorio/pi0/rgd-pi0/farm_v100_IB   # Stage A output, inbending (input to B and C)
OUT=analysis/v100_IB                                     # everything this script writes
GA=config/binning/grid_A_q2_xb.json
GB=config/binning/grid_B_z_pt2.json
declare -A NRUNS=([LD2]=24 [CxC]=13 [Cu]=27 [Sn]=27)     # config/runs.json, inbending; Cu/Sn without special DC
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

# -------------------------------------------------------------------------------------
# Inputs: every target complete, and nothing outside the inbending production
# -------------------------------------------------------------------------------------
stamp "Inputs"
for t in LD2 CxC Cu Sn; do
  n=$(find $IB/$t -mindepth 1 -maxdepth 1 -type d -regex '.*/[0-9]+' | wc -l)
  nf=$(find $IB/$t -name 'slim_*.root' | wc -l)
  echo "$t: $n runs, $nf slims (expected ${NRUNS[$t]} runs)"
  if [ "$n" -ne "${NRUNS[$t]}" ]; then echo "expected ${NRUNS[$t]} runs for $t, found $n -- stopping"; exit 1; fi
done
mkdir -p $OUT

# -------------------------------------------------------------------------------------
# B. Grid (only with MAKE_IB_GRID=1)
# -------------------------------------------------------------------------------------
if [ "${MAKE_IB_GRID:-0}" = "1" ]; then
  stamp "B. Grid: make_grid on the IB slims"
  args=()
  for t in LD2 CxC Cu Sn; do
    while IFS= read -r f; do args+=(--input "$f"); done < <(find $IB/$t -name 'slim_*.root' | sort)
  done
  echo "input files: $(( ${#args[@]} / 2 ))"
  GA=$OUT/grid_A_q2_xb_IB.json
  ./build/src/tools/make_grid/make_grid "${args[@]}" --config config/cuts.json --out-a $GA --out-b $OUT/grid_B_z_pt2_equal_stat_IB.json --threads 4
  python3 analysis/scripts/pin_q2_edge.py $GA 2.0
else
  stamp "B. Grid: none made, using $GA and $GB (the outbending bins)"
fi

# -------------------------------------------------------------------------------------
# C. Bin
# -------------------------------------------------------------------------------------
mkdir -p $OUT/binned
for t in LD2 CxC Cu Sn; do
  args=()
  while IFS= read -r f; do args+=(--input "$f"); done < <(find $IB/$t -name 'slim_*.root' | sort)
  stamp "C. Bin: $t ($(( ${#args[@]} / 2 )) files)"
  ./build/src/stageB_bin/stageB_bin "${args[@]}" --output $OUT/binned/$t.root --config config/cuts.json --grid-a $GA --grid-b $GB --threads 4
done

# -------------------------------------------------------------------------------------
# D. Extract
# -------------------------------------------------------------------------------------
export PYTHONPATH=python
G="--grid-a $GA --grid-b $GB"
mkdir -p $OUT/results

for t in CxC Cu Sn; do
  stamp "D.1 Extract, ratio: $t"
  python3 -m pi0.ratio --target $OUT/binned/$t.root --ld2 $OUT/binned/LD2.root --config config $G --out $OUT/results/ratio_$t.csv --allow-unpublishable
done

for t in CxC Cu Sn; do
  stamp "D.2 Extract, broadening: $t"
  python3 -m pi0.broadening --target $OUT/binned/$t.root --ld2 $OUT/binned/LD2.root --config config $G --out $OUT/results/broadening_$t.csv --allow-unpublishable
done

for t in LD2 CxC Cu Sn; do
  stamp "D.3 Extract, bsa: $t"
  python3 -m pi0.bsa --file $OUT/binned/$t.root --config config $G --out $OUT/results/bsa_$t.csv --allow-unpublishable
done

stamp "D.4 Extract, qa: LD2"
python3 -m pi0.qa --file $OUT/binned/LD2.root --config config $G --outdir $OUT/figs_qa --allow-unpublishable

stamp "D.5 Extract, plots_results"
python3 -m pi0.plots_results --results $OUT/results --outdir $OUT/figs_results

stamp "ALL DONE"
