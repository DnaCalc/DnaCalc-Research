from ask import *
import json
cases=[]
for l in open("tmp/BINOM.INV.trace"):
    r=json.loads(l)
    if not r["expected"].startswith("0x"): continue
    cases.append(([f(a[2:]) for a in r["args"]], f(r["expected"][2:]), r["actual"]==r["expected"]))
rows=[];idx=[]
for i,(a,e,ok) in enumerate(cases):
    n=int(a[0])
    for k in range(0,n+1):
        rows.append([k,a[0],a[1],boolean(True)]); idx.append((i,k))
res=ask("BINOM.DIST",rows)["rows"]
cdf={}
for (i,k),x in zip(idx,res): cdf[(i,k)]=f(x["bits"])
good=0; bad=[]
for i,(a,e,ok) in enumerate(cases):
    n=int(a[0]); alpha=a[2]
    k=next((k for k in range(n+1) if cdf[(i,k)]>=alpha), n)
    if k==e: good+=1
    else: bad.append((a,e,k,[cdf[(i,j)] for j in range(n+1)]))
print(len(cases),"smallest k with BINOM.DIST>=alpha matches",good)
for b in bad[:10]: print(b)
