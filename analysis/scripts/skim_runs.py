import os, re, csv, rcdb

base = "/cache/clas12/rg-d/production/pass1/recon/%s/dst/train/SIDIS/"
db = rcdb.RCDBProvider("mysql://rcdb@clasdb.jlab.org/rcdb")
conds = ["run_start_time", "event_count", "beam_energy", "beam_current_request",
         "half_wave_plate", "torus_scale", "solenoid_scale", "target", "run_type"]

w = csv.writer(open("skim_runs.csv", "w"))
w.writerow(["skim_target", "run", "file_size_GB"] + conds)

for t in ["LD2", "CxC", "CuSn"]:
    for f in sorted(os.listdir(base % t)):
        m = re.search(r"(\d{6})", f)
        if not m:
            continue
        run = int(m.group(1))
        size = round(os.path.getsize(base % t + f) / 1e9, 1)
        row = [t, run, size]
        for c in conds:
            try:
                row.append(db.get_condition(run, c).value)
            except Exception:
                row.append("")
        w.writerow(row)

print("wrote skim_runs.csv")
