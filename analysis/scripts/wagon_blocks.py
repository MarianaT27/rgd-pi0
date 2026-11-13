import numpy as np
import uproot

t = uproot.open("/work/clas12/mtenorio/pi0/stageA_v021/LD2/018431/mt_rgd_pi0_stageA_v021_LD2_job_00000/slim_0.root")["events"].arrays(["event"], library="np")["event"]
w = uproot.open("/work/clas12/mtenorio/pi0/wagon_test/slim_wagon_018431.root")["events"].arrays(["event"], library="np")["event"]
in_w = np.isin(t, w)
blk = t // 1_000_000
print("block (M events)   train   found in wagon   fraction")
for b in np.unique(blk):
    m = blk == b
    print(f"{b:6d}   {m.sum():10d}   {in_w[m].sum():10d}   {in_w[m].mean():.3f}")
