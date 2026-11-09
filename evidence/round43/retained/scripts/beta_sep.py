import json, struct, os, numpy as np
import racer
racer.RACER=os.path.expanduser("~/codex-runs/r43/racer-beta")
rng=np.random.default_rng(4330)
Q=[]
for k in range(20000):
    x=float(np.exp(rng.uniform(np.log(1e-300),np.log(0.0999))))
    b=float(np.exp(rng.uniform(np.log(1e-3),np.log(1e6))))
    Q.append((x,1.0,b))
res={}
for v in [6,7]:
    os.environ["R43_BETA_A1"]=str(v)
    res[v]=racer.invoke([("BETA.DIST",[x,a,b,"T"]) for x,a,b in Q])
d=[q for q,a,b in zip(Q,res[6],res[7]) if a!=b]
print("differ",len(d),"of",len(Q)); print(d[:5])
json.dump(d,open("beta_sep.json","w"))
