#!/bin/bash
# =====================================================================================
# run_ndis_chain.sh : stages C and D of the rgd-pi0 chain on the INCLUSIVE-N_DIS production
#
# Same as run_full_chain.sh, except:
#   A. Skim      input is farm_ndis/<target>/<run>/<job>/slim_0.root: Stage A with NO photon
#                requirement, so N_DIS counts every DIS event (94 slims: run 18540 is missing,
#                its file holds 10 events of run 18539 and Stage A refuses multi-run files)
#   B. Grid      NOT run: the committed grid in config/binning is reused, so every bin
#                matches analysis/full and the effect of the N_DIS fix can be compared bin by bin
#   C. Bin       -> analysis/ndis/binned_full/<target>.root  (also holds the x_F trees)
#   D. Extract   as in run_full_chain.sh, into analysis/ndis/, plus
#                D.6 x_F vs z    xf_z.py                     -> analysis/ndis/xf/
#
# Run (from anywhere, it moves to the repo root itself):
#   nohup bash analysis/scripts/run_ndis_chain.sh >& analysis/logs/run_ndis_chain.log &
# =====================================================================================
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0

FULL=/volatile/clas12/mtenorio/pi0/rgd-pi0/farm_ndis    # stage A output (input to C)
OUT=analysis/ndis                                       # everything this script writes
GA=config/binning/grid_A_q2_xb.json
GB=config/binning/grid_B_z_pt2.json
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

# B. Grid: NOT run, the committed grid is reused (see header)
if ! git diff --quiet -- config/binning/; then
  echo "config/binning differs from the committed grid -- restore it first: git checkout -- config/binning/"
  exit 1
fi

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

# D.6 x_F vs z
stamp "D.6 x_F vs z"
python3 analysis/scripts/xf_z.py $OUT/binned_full $OUT/xf

stamp "ALL DONE"
