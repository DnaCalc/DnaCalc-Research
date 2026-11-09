import sys, json, struct, time, math, numpy as np
sys.path.insert(0, "/home/ubuntu/projects/DnaCalc/DnaCalc/tools/excel-service")
from client import Client, number
def bits(x): return struct.unpack("<Q",struct.pack("<d",x))[0]
def f(b): return struct.unpack("<d",struct.pack("<Q",b))[0]
c=Client(timeout=120)
centers={"t0.7071+":0.5+0.35355,"tsqrt+":0.5+1/(2*math.sqrt(2)),"t0.7071-":0.5-0.35355,"tsqrt-":0.5-1/(2*math.sqrt(2))}
out={}
for k,p in centers.items():
    b=bits(p); rows=[[number("%016x"%(b-1024+j))] for j in range(2048)]
    t0=time.time(); r=c.ask("NORM.S.INV",rows,priority="bulk"); dt=time.time()-t0
    out[k]=dict(start="%016x"%(b-1024),rows=r["rows"],seconds=dt,build=r["excel_build"],env=r["env_id"])
    print(k,dt,len(r["rows"]))
json.dump(out,open("/home/ubuntu/codex-runs/r43/fresh/NORM.S.INV.split.replies.json","w"))
