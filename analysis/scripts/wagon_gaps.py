import numpy as np
import uproot

t = uproot.open("/work/clas12/mtenorio/pi0/stageA_v021/LD2/018431/mt_rgd_pi0_stageA_v021_LD2_job_00000/slim_0.root")["events"].arrays(["event", "q2"], library="np")
w = uproot.open("/work/clas12/mtenorio/pi0/wagon_test/slim_wagon_018431.root")["events"].arrays(["event"], library="np")["event"]

# only the train events the wagon should have kept
exp = (t["q2"] > 2.0) | (t["event"] % 5 == 0)
ev = np.sort(t["event"][exp])
miss = ~np.isin(ev, w)
print(f"expected {len(ev)}, missing {miss.sum()} ({miss.mean():.3f})")

# stretches of consecutive missing expected events
d = np.diff(np.concatenate([[0], miss.astype(int), [0]]))
start, stop = np.where(d == 1)[0], np.where(d == -1)[0]
n = stop - start
print(f"missing stretches: {len(n)}, median length {np.median(n):.0f} events")
print("longest 20 (first event, last event, events missing):")
for i in np.argsort(n)[::-1][:20]:
    print(f"   {ev[start[i]]:12d}  {ev[stop[i] - 1]:12d}  {n[i]:8d}")
