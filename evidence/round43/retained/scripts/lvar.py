import pickle, itertools, sys, collections
sys.path.insert(0,"/home/ubuntu/codex-runs/r43")
from qr2 import ls2
from multiprocessing import Pool
cases=pickle.load(open("/home/ubuntu/codex-runs/r43/dfcases.pkl","rb"))
opts=[dict(metric=m,tie=t,zero=z,move=mv,onstop=o) for m in ["run","runclip","fresh","abs"] for t in ["first","last","const","noconst"] for z in ["compete","exclude"] for mv in ["swap","rot"] for o in ["stop","continue"]]
def ev(o):
    ok=0; bad=collections.Counter()
    for c in cases:
        k=len(ls2([list(x) for x in c["cols"]],c["y"],c["const"],o))
        if k==c["count"]: ok+=1
        else: bad[c["src"][:12]]+=1
    return ok,o,dict(bad)
with Pool(16) as p: res=p.map(ev,opts)
res.sort(key=lambda r:-r[0])
print(len(cases))
for r in res[:20]: print(r)
base=[r for r in res if r[1]==dict(metric="run",tie="const",zero="compete",move="swap",onstop="stop")]; print("current",base)
