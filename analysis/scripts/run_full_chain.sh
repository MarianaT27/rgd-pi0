#!/bin/bash
# Full chain on the full production: make_grid -> Stage B -> QA, ratio, broadening, BSA -> plots.
# Stops at the first failing step (set -e). Every step is timestamped in the log.
set -euo pipefail
cd /volatile/clas12/mtenorio/pi0/rgd-pi0
FULL=/volatile/clas12/mtenorio/pi0/rgd-pi0/farm_full
GA=config/binning/grid_A_q2_xb.json
GB=config/binning/grid_B_z_pt2.json
stamp() { echo; echo "=== $(date '+%F %T')  $* ==="; }

stamp "make_grid"
args=()
for f in $FULL/*/*/*/slim_0.root; do args+=(--input "$f"); done
n=$(( ${#args[@]} / 2 ))
echo "input files: $n"
if [ "$n" -ne 95 ]; then echo "expected 95 files, found $n -- stopping"; exit 1; fi
./build/src/tools/make_grid/make_grid "${args[@]}" --config config/cuts.json --out-a $GA --out-b $GB --threads 4

mkdir -p binned_full
for t in LD2 CxC Cu Sn; do
  args=()
  for f in $FULL/$t/*/*/slim_0.root; do args+=(--input "$f"); done
  stamp "Stage B $t ($(( ${#args[@]} / 2 )) files)"
  ./build/src/stageB_bin/stageB_bin "${args[@]}" --output binned_full/$t.root --config config/cuts.json --grid-a $GA --grid-b $GB --allow-truncated-inputs --threads 4
done

export PYTHONPATH=python
stamp "QA (LD2)"
python3 -m pi0.qa --file binned_full/LD2.root --config config --outdir figs_qa_full --allow-unpublishable

mkdir -p results_full
for t in CxC Cu Sn; do
  stamp "ratio + broadening $t"
  python3 -m pi0.ratio --target binned_full/$t.root --ld2 binned_full/LD2.root --config config --out results_full/ratio_$t.csv --allow-unpublishable
  python3 -m pi0.broadening --target binned_full/$t.root --ld2 binned_full/LD2.root --config config --out results_full/broadening_$t.csv --allow-unpublishable
done

for t in LD2 CxC Cu Sn; do
  stamp "BSA $t"
  python3 -m pi0.bsa --file binned_full/$t.root --config config --out results_full/bsa_$t.csv --allow-unpublishable
done

stamp "plots_results"
python3 -m pi0.plots_results --results results_full --outdir figs_results_full

stamp "ALL DONE"
