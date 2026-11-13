#!/bin/bash
# =====================================================================================
# run_hermes_chain.sh : RG-D pi0 R_A in the HERMES 1D bins, with the HERMES cuts (v0.2.4)
#
# Reference: Airapetian et al. (HERMES), Nucl. Phys. B 780 (2007) 1, arXiv:0704.3270.
# The HERMES cuts are applied to THIS plot only:
#   nu > 6 GeV and z > 0.2 : bin edges of the HERMES grids, applied in hermes_ra.py by bin selection
#   x_F > 0                : stageB_bin --xf-min 0, per pi0 pair (same and mixed events)
# Nothing here touches the production grids or binned files.
#
#   C. Bin       stageB_bin with the HERMES grids and --xf-min 0, per target
#                -> analysis/hermes/binned/<target>.root
#   D. Plot      hermes_ra.py: summed spectra per HERMES bin, production extraction, R_A
#                -> analysis/hermes/hermes_ra.csv, hermes_ra.pdf/.png
#
# Input: the v0.2.1 Stage A slims (no Stage A rerun). To use another production, change FULL.
# Run (from anywhere, it moves to the repo root itself):
#   nohup bash analysis/scripts/run_hermes_chain.sh >& analysis/logs/run_hermes_chain.log &
#
# Stops at the first failing step (set -e). Re-running overwrites analysis/hermes/.
# =====================================================================================
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0

FULL=/work/clas12/mtenorio/pi0/stageA_v021    # Stage A output v0.2.1 (input to C)
OUT=analysis/hermes
GA=config/binning/grid_A_q2_nu_hermes.json
GB=config/binning/grid_B_z_pt2_hermes.json
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

# -------------------------------------------------------------------------------------
# C. Bin, HERMES grids, x_F > 0
# -------------------------------------------------------------------------------------
mkdir -p $OUT/binned
for t in LD2 CxC Cu Sn; do
  args=()
  for f in $FULL/$t/*/*/slim_0.root; do args+=(--input "$f"); done
  stamp "C. Bin (HERMES grids, x_F > 0): $t ($(( ${#args[@]} / 2 )) files)"
  ./build/src/stageB_bin/stageB_bin "${args[@]}" --output $OUT/binned/$t.root --config config/cuts.json \
      --grid-a $GA --grid-b $GB --xf-min 0 --allow-truncated-inputs --threads 4
done

# -------------------------------------------------------------------------------------
# D. R_A in the HERMES bins
# -------------------------------------------------------------------------------------
stamp "D. hermes_ra.py"
python3 analysis/scripts/hermes_ra.py $OUT/binned $OUT

stamp "ALL DONE"
