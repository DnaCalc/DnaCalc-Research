import sys, json, struct, time, numpy as np
sys.path.insert(0, "/home/ubuntu/projects/DnaCalc/DnaCalc/tools/excel-service")
from client import Client, number, boolean
from racer import invoke
def bits(x): return struct.unpack("<Q",struct.pack("<d",float(x)))[0]
AB=[0.3,0.5,0.9,1.0,1.5,2.0,2.5,3.0,4.7,7.0,10.0,20.0,50.0,100.0,500.0]
xs=sorted(set([float(v) for v in np.concatenate([np.logspace(-6,-0.31,20),np.linspace(0.05,0.95,24),1-np.logspace(-6,-0.31,20)])]))
Q=[(x,a,b) for a in AB for b in AB for x in xs]
c=Client(timeout=600)
rows=[[number(x),number(a),number(b),boolean(True)] for x,a,b in Q]
t0=time.time(); ans=[]
for k in range(0,len(rows),5000):
    r=c.ask("BETA.DIST",rows[k:k+5000],priority="bulk"); ans+=r["rows"]
dt=time.time()-t0
ours=invoke([("BETA.DIST",[x,a,b,"T"]) for x,a,b in Q])
json.dump(dict(questions=Q,answers=ans,ours=[o if isinstance(o,str) else "%016x"%bits(o) for o in ours],seconds=dt,build=r["excel_build"]),open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.map.json","w"))
print("asked",len(Q),"in %.1fs"%dt)
