from ask import *
import json, sys
sys.path.insert(0,"/home/koderbot/projects/DnaCalc/DnaCalc-Research/racer")
from fractions import Fraction as Q
from intervals import B64,E64
cases=[]
for l in open("tmp/BINOM.INV.trace"):
    r=json.loads(l)
    if not r["expected"].startswith("0x"): continue
    a=[f(x[2:]) for x in r["args"]]
    if a[0]>200: continue
    cases.append((a, f(r["expected"][2:])))
rows=[];idx=[]
for i,(a,e) in enumerate(cases):
    n=int(a[0])
    for k in range(0,n+1):
        rows.append([k,a[0],a[1],boolean(False)]); idx.append((i,k))
res=ask("BINOM.DIST",rows)["rows"]
pmf={}
for (i,k),x in zip(idx,res): pmf[(i,k)]=f(x["bits"])
tally={}
for i,(a,e) in enumerate(cases):
    n=int(a[0]); alpha=a[1+1]
    accs={"d":0.0}
    cd=[];cx=[];s=0.0;sx=Q(0)
    for k in range(n+1):
        s+=pmf[(i,k)]; cd.append(s)
        sx=E64.rn(sx+Q(pmf[(i,k)])); cx.append(float(B64.rn(sx)))
    for name,c in [("dbl",cd),("x87",cx)]:
        for cmp in (">=",">"):
            k=next((k for k in range(n+1) if (c[k]>=alpha if cmp==">=" else c[k]>alpha)), n)
            key=name+cmp; tally[key]=tally.get(key,0)+(k==e)
    # reverse: 1 - upper tail
print(len(cases),tally)
