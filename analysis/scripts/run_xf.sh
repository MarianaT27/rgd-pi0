#!/bin/bash
# run_xf.sh : step C (binning) on the full sample with the x_F trees, then xf_z.py.
# Output: analysis/xf/binned_full/<T>.root and the x_F plots in analysis/xf/
# Run (from the repo root):
#   nohup bash analysis/scripts/run_xf.sh >& analysis/logs/run_xf.log &
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0

FULL=/volatile/clas12/mtenorio/pi0/rgd-pi0/farm_full
GA=config/binning/grid_A_q2_xb.json
GB=config/binning/grid_B_z_pt2.json
OUT=analysis/xf
mkdir -p $OUT/binned_full
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

for t in LD2 CxC Cu Sn; do
  args=()
  for f in $FULL/$t/*/*/slim_0.root; do args+=(--input "$f"); done
  stamp "C. Bin: $t ($(( ${#args[@]} / 2 )) files)"
  ./build/src/stageB_bin/stageB_bin "${args[@]}" --output $OUT/binned_full/$t.root --config config/cuts.json --grid-a $GA --grid-b $GB --allow-truncated-inputs --threads 4
done

stamp "x_F vs z"
python3 analysis/scripts/xf_z.py $OUT/binned_full $OUT

stamp "ALL DONE"