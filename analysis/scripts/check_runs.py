import json, os
d = json.load(open('config/runs.json'))
base = '/cache/clas12/rg-d/production/pass1/recon'
for t in ['LD2', 'CxC', 'CuSn']:
    have = set(int(r) for r in os.listdir(f'{base}/{t}/dst/recon') if r.isdigit())
    lst = set(d['outbending'][t]) | set(d['inbending'][t])
    print(t, 'on disk but not in runs.json:', sorted(have - lst))
    print(t, 'in runs.json but not on disk:', sorted(lst - have))