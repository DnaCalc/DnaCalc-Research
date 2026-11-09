import sys, json, struct, time, numpy as np
sys.path.insert(0, "/home/ubuntu/projects/DnaCalc/DnaCalc/tools/excel-service")
from client import Client, number
def bits(x): return struct.unpack("<Q",struct.pack("<d",x))[0]
c=Client(timeout=600)
ts=np.round(np.arange(0.40,0.7068,0.002),6)
out=[]; t0=time.time()
reqs=[]
for t in ts:
    p=0.5+t/2; b=bits(p)
    reqs.append(("NORM.S.INV",[[number("%016x"%(b+j))] for j in range(512)]))
replies=c.ask_many(reqs) if hasattr(c,"ask_many") else [c.ask(*r) for r in reqs]
for t,(fn,rows),r in zip(ts,reqs,replies):
    out.append(dict(t=float(t),start="%016x"%bits(0.5+t/2),rows=[row.get("bits") if row.get("type")=="number" else None for row in r["rows"]]))
json.dump(dict(windows=out,build=replies[0]["excel_build"],env=replies[0]["env_id"],seconds=time.time()-t0),open("/home/ubuntu/codex-runs/r43/fresh/NORM.S.INV.lattice.replies.json","w"))
print(len(out),"windows in %.1fs"%(time.time()-t0))
