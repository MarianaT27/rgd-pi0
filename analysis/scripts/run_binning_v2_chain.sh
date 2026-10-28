#!/bin/bash
# =====================================================================================
# run_binning_v2_chain.sh : stages C and D with the hand-chosen Grid B (binning_clas6.patch)
#
# Input: the inclusive-N_DIS slims in farm_ndis (same as run_ndis_chain.sh). Stage A is
# not rerun: the grid does not enter Stage A.
#
# Two passes, one after the other:
#   1. MAIN       Grid A (Q2, xB)  x  Grid B (z, pT2), committed grids
#                 -> analysis/binning_v2/   binned_full, results_full, figs, xf
#   2. CLAS6 (Q2, nu) comparison binning: Grid A = config/binning/grid_A_q2_nu.json
#                 -> analysis/clas6_q2nu/   binned_full, results_full (ratio, broadening)
#                 then clas6_compare.py     -> analysis/clas6_q2nu/compare/
#
# make_grid is NOT run: it would overwrite the hand-chosen Grid B with equal-statistics
# edges. The script refuses to start if config/binning differs from the committed version.
#
# Run (from anywhere, it moves to the repo root itself):
#   nohup bash analysis/scripts/run_binning_v2_chain.sh >& analysis/logs/run_binning_v2_chain.log &
# =====================================================================================
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0

FULL=/volatile/clas12/mtenorio/pi0/rgd-pi0/farm_ndis    # stage A output (input to C)
GA=config/binning/grid_A_q2_xb.json
GA_NU=config/binning/grid_A_q2_nu.json
GB=config/binning/grid_B_z_pt2.json
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

if ! git diff --quiet -- config/binning/ || [ -n "$(git status --porcelain -- config/binning/)" ]; then
  echo "config/binning has uncommitted changes -- commit the binning patch first (git status config/binning)"
  exit 1
fi
export PYTHONPATH=python

# stageB over all four targets: bin_all <grid_A> <out_dir>
bin_all() {
  mkdir -p $2/binned_full
  for t in LD2 CxC Cu Sn; do
    args=()
    for f in $FULL/$t/*/*/slim_0.root; do args+=(--input "$f"); done
    stamp "C. Bin ($2): $t ($(( ${#args[@]} / 2 )) files)"
    ./build/src/stageB_bin/stageB_bin "${args[@]}" --output $2/binned_full/$t.root --config config/cuts.json --grid-a $1 --grid-b $GB --allow-truncated-inputs --threads 4
  done
}

# ratio and broadening for C, Cu, Sn: ratio_broadening <grid_A> <out_dir>
ratio_broadening() {
  mkdir -p $2/results_full
  for t in CxC Cu Sn; do
    stamp "D.1 ratio ($2): $t"
    python3 -m pi0.ratio --target $2/binned_full/$t.root --ld2 $2/binned_full/LD2.root --config config --grid-a $1 --grid-b $GB --out $2/results_full/ratio_$t.csv --allow-unpublishable
    stamp "D.2 broadening ($2): $t"
    python3 -m pi0.broadening --target $2/binned_full/$t.root --ld2 $2/binned_full/LD2.root --config config --grid-a $1 --grid-b $GB --out $2/results_full/broadening_$t.csv --allow-unpublishable
  done
}

# -------------------------------------------------------------------------------------
# 1. MAIN binning: (Q2, xB) x new Grid B
# -------------------------------------------------------------------------------------
OUT=analysis/binning_v2
bin_all $GA $OUT
ratio_broadening $GA $OUT
for t in LD2 CxC Cu Sn; do
  stamp "D.3 bsa: $t"
  python3 -m pi0.bsa --file $OUT/binned_full/$t.root --config config --out $OUT/results_full/bsa_$t.csv --allow-unpublishable
done
stamp "D.4 qa: LD2"
python3 -m pi0.qa --file $OUT/binned_full/LD2.root --config config --outdir $OUT/figs_qa_full --allow-unpublishable
stamp "D.5 plots_results"
python3 -m pi0.plots_results --results $OUT/results_full --outdir $OUT/figs_results_full
stamp "D.6 x_F vs z"
python3 analysis/scripts/xf_z.py $OUT/binned_full $OUT/xf

# -------------------------------------------------------------------------------------
# 2. CLAS6 comparison binning: (Q2, nu) x new Grid B
# -------------------------------------------------------------------------------------
OUT=analysis/clas6_q2nu
bin_all $GA_NU $OUT
ratio_broadening $GA_NU $OUT
stamp "E. comparison with Mineeva et al. (CLAS6 pi0)"
python3 analysis/scripts/clas6_compare.py $OUT/binned_full $OUT/compare

stamp "ALL DONE"
