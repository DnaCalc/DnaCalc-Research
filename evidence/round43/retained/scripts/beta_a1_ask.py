import sys, json, numpy as np, time
sys.path.insert(0, "/home/ubuntu/projects/DnaCalc/DnaCalc/tools/excel-service")
from client import Client, number, boolean
rng=np.random.default_rng(4321)
Q=[]
for b in [0.3,0.7,1.5,2.0,3.0,5.5,10.0,37.0,100.0,1000.0]:
    for x in np.exp(rng.uniform(np.log(1e-9),np.log(0.0999),300)): Q.append((float(x),1.0,b))
c=Client(timeout=300); t0=time.time()
r=c.ask("BETA.DIST",[[number(x),number(a),number(b),boolean(True)] for x,a,b in Q],priority="bulk")
json.dump(dict(questions=Q,answers=r["rows"],build=r["excel_build"],seconds=time.time()-t0),open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.a1.json","w"))
print(len(Q), time.time()-t0)
