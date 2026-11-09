from ask import *
import random, json
random.seed(21)
rows=[]
while len(rows)<400:
    r=random.choice([random.uniform(0.001,0.03),random.uniform(0.03,0.3)])
    n=random.uniform(2,400)
    if random.random()<0.3: n=float(int(n))
    s=random.uniform(1,n); e=random.uniform(s,n)
    if random.random()<0.3: s=float(int(s)) if int(s)>=1 else 1.0
    pv=round(random.uniform(100,1e6),2); t=random.choice([0,1])
    rows.append([r,n,pv,s,e,t])
# edge designs
rows+= [[0.05,10.5,1000,10,10.5,0],[0.05,10.5,1000,1,1,1],[0.05,10.5,1000,1.2,1.8,0],[0.05,10.5,1000,0.99,2,0],[0.05,10.5,1000,2,10.6,0],[0.05,1,1000,1,1,0],[0.05,1.5,1000,1,1.5,1]]
out={}
for fn in ["CUMPRINC","CUMIPMT"]:
    out[fn]=ask(fn,rows)["rows"]
json.dump(dict(rows=[[h(float(v)) for v in r] for r in rows],answers=out),open("frac1.json","w"))
print(len(rows))
