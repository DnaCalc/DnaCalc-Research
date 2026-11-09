from ask import *
import random, json
random.seed(11)
rows=[];meta=[]
for i in range(30):
    r=random.choice([random.uniform(0.0005,0.02),random.uniform(0.02,0.3)])
    n=random.randint(5,400); pv=round(random.uniform(100,1e6),2)
    ks=sorted(random.sample(range(1,n+1),min(n,40)))
    for k in ks: rows.append([r,n,pv,k,k,0])
a=ask("CUMPRINC",rows)["rows"]
with open("d1.txt","w") as o:
    for row,x in zip(rows,a):
        if "bits" in x: o.write(" ".join(h(float(v)) for v in row)+" "+x["bits"]+"\n")
print(len(rows))
