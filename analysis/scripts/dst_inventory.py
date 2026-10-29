"""
dst_inventory.py : how many DST recon files are on /cache, per run and per target,
for the runs of the full production (the run folders in farm_full/<target>/).

Usage (from the repo root):
    python3 analysis/scripts/dst_inventory.py

Reads:   farm_full/<target>/<run>/                       (which runs the full production used)
         /cache/clas12/rg-d/production/pass1/recon/<LD2|CxC|CuSn>/dst/recon/<run>/*.hipo
Writes:  analysis/runinfo/dst_inventory.csv              (target, run, n_files, size_GB)
Prints:  one line per run, then a summary per target.

Cu and Sn come from the same CuSn files, so they show the same counts.
Only looks at /cache: files that are only on tape (/mss) are not counted.
"""
import csv
import os

BASE = "/cache/clas12/rg-d/production/pass1/recon"
SRC = {"LD2": "LD2", "CxC": "CxC", "Cu": "CuSn", "Sn": "CuSn"}
rows = []
for t, src in SRC.items():
    runs = sorted(int(r) for r in os.listdir(f"farm_full/{t}") if r.isdigit())
    for run in runs:
        d = f"{BASE}/{src}/dst/recon/{run:06d}"
        files = [f for f in os.listdir(d) if f.endswith(".hipo")] if os.path.isdir(d) else []
        gb = sum(os.path.getsize(os.path.join(d, f)) for f in files) / 1e9
        rows.append((t, run, len(files), gb))
        print(f"{t:4s} {run}  {len(files):5d} files  {gb:8.1f} GB")

os.makedirs("analysis/runinfo", exist_ok=True)
with open("analysis/runinfo/dst_inventory.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["target", "run", "n_files", "size_GB"])
    w.writerows((t, r, n, f"{g:.2f}") for t, r, n, g in rows)

print("\ntarget  runs  runs_with_files  total_files  total_GB  min_files  max_files")
for t in SRC:
    s = [r for r in rows if r[0] == t]
    n = [r[2] for r in s]
    print(f"{t:6s}  {len(s):4d}  {sum(1 for x in n if x > 0):15d}  {sum(n):11d}  {sum(r[3] for r in s):8.1f}"
          f"  {min(n):9d}  {max(n):9d}")
print("\nwrote analysis/runinfo/dst_inventory.csv")
