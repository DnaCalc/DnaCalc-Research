from ask import *
import sys, random
sys.path.insert(0,"/home/"+__import__("getpass").getuser()+"/projects/DnaCalc/DnaCalc-Research/racer")
from intervals import B64,E64
from fractions import Fraction as Q
random.seed(1)
loans=[(0.05/12,360,200000.0),(0.09/12,360,125000.0),(0.1,30,1000.0),(0.013,100,5e5),(0.2,50,777.0)]
def sums(xs):
    out={}
    s=0.0
    for x in xs: s+=x
    out["l2r"]=s
    s=0.0
    for x in reversed(xs): s+=x
    out["r2l"]=s
    s=Q(0)
    for x in xs: s=E64.rn(s+Q(x))
    out["x87acc"]=float(B64.rn(s))
    s=Q(0)
    for x in reversed(xs): s=E64.rn(s+Q(x))
    out["x87acc_r"]=float(B64.rn(s))
    out["exact"]=float(B64.rn(sum(Q(x) for x in xs)))
    return out
tally={}
tot=0
for (r,n,pv) in loans:
    sing=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv,k,k,0] for k in range(1,n+1)])["rows"]]
    ranges=[]
    for _ in range(60):
        s=random.randint(1,n); e=random.randint(s,n); ranges.append((s,e))
    got=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv,s,e,0] for s,e in ranges])["rows"]]
    for (s,e),g in zip(ranges,got):
        if e==s: continue
        tot+=1
        for k,v in sums(sing[s-1:e]).items():
            tally[k]=tally.get(k,0)+(v==g)
print(tot,tally)
