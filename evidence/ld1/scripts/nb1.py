from ask import *
import json
rows=[];exp=[]
for l in open("tmp/NEGBINOM.DIST.trace"):
    r=json.loads(l)
    if not r["expected"].startswith("0x"): continue
    a=[f(x[2:]) for x in r["args"]]
    rows.append(a); exp.append((r["expected"][2:], r["actual"]==r["expected"], r["actual"]))
cum=[a for a in rows if a[3]==1]
print(len(rows), len(cum))
q=[[a[0],a[1],a[2],boolean(a[3]==1)] for a in rows]
nb=ask("NEGBINOM.DIST",q)["rows"]
be=ask("BETA.DIST",[[a[2],int(a[1]),int(a[0])+1,boolean(True)] for a in rows])["rows"]
same_live=sum(1 for x,e in zip(nb,exp) if x.get("bits")==e[0])
rel=0; relc=0; ourok=0; ourbad_rel=0
for a,x,b,e in zip(rows,nb,be,exp):
    if a[3]==1:
        relc+=1
        if x.get("bits")==b.get("bits"): rel+=1
print("live==bank",same_live,"of",len(rows)," cumulative: NEGBINOM==BETA.DIST(p,r,k+1)",rel,"of",relc)
