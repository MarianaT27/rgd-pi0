import rcdb
from collections import defaultdict

db = rcdb.RCDBProvider("mysql://rcdb@clasdb-farm.jlab.org/rcdb")
rows = db.select_runs("", 18300, 19140).get_values(["beam_energy", "target"], insert_run_number=True)

groups = defaultdict(list)
for run, energy, target in rows:
    if energy is None:
        continue
    groups[energy].append(run)

for e in sorted(groups):
    runs = groups[e]
    print(f"{e} MeV: {len(runs)} runs, {min(runs)} to {max(runs)}")