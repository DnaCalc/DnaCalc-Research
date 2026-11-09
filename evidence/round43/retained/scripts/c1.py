from nl import *
from racer import invoke
d=load(); t=d["t"]; x=d["x"]
m=np.nonzero((t>0.8505)&(t<=0.875))[0][::40]
P=preimages(x[m])
rows=[];meta=[]
for i,k in enumerate(m):
    ys=[v for v in P[i] if not np.isnan(v)]
    y0=ys[0]
    cand=[np.nextafter(y0,-1) if j<0 else y0 for j in [0]]
    cs=[y0]
    for j in range(1,4): cs.append(float(np.nextafter(cs[-1],2))); 
    c2=[y0]
    for j in range(1,4): c2.append(float(np.nextafter(c2[-1],-2)))
    allc=sorted(set(c2+cs))
    for c in allc: rows.append(("ERF",[c])); meta.append((i,c,c in ys))
out=invoke(rows)
best=collections=None
import collections
res=collections.defaultdict(list)
for (i,c,isx),e in zip(meta,out): res[i].append((abs(LD(e)-LD(t[m[i]])),c,isx,e))
hit=0;tot=0
for i,v in res.items():
    v.sort(key=lambda a:a[0]); tot+=1; hit+=v[0][2]
print("excel y minimises |ERF(y)-t| in",hit,"of",tot)
