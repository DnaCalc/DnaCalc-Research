import pickle, sys, math, collections
sys.path.insert(0,"/home/ubuntu/codex-runs/r43")
from qr2 import ls2
cases=pickle.load(open("/home/ubuntu/codex-runs/r43/dfcases.pkl","rb"))
rows=[]
for c in cases:
    o=dict(metric="runclip",tie="first",zero="compete",move="swap",onstop="stop",oneslast=True,ratios=[])
    k=ls2([list(x) for x in c["cols"]],c["y"],c["const"],o)
    zero=any(all(v==0 for v in col) for col in c["cols"])
    for name,ratio,j in o["ratios"]:
        if 0<ratio<2**-36:
            rows.append((math.log2(ratio),len(c["y"]),c["const"],zero,len(k)==c["count"],len(k),c["count"],j,len(c["cols"])+c["const"]))
rows.sort()
for r in rows:
    if r[0]>-46: print("%.3f n=%d const=%d zero=%d ok=%d emu=%d excel=%d j=%d p=%d"%r)
