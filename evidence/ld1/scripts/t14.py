from ask import *
import sys, random
sys.path.insert(0,"/home/koderbot/projects/DnaCalc/DnaCalc-Research/racer")
from intervals import B64,E64
from fractions import Fraction as Q
random.seed(3)
tally={};tot=0
for (r,n,pv) in [(1.0,40,12345.67),(1.0,60,987.65),(3.0,25,777.77),(1.0,50,0.3)]:
    sing=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv,k,k,0] for k in range(1,n+1)])["rows"]]
    ranges=set()
    while len(ranges)<120:
        s=random.randint(1,n); e=random.randint(s,n)
        if e>s: ranges.add((s,e))
    ranges=sorted(ranges)
    got=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv,s,e,0] for s,e in ranges])["rows"]]
    for (s,e),g in zip(ranges,got):
        xs=sing[s-1:e]; tot+=1
        c={}
        a=0.0
        for x in xs: a+=x
        c["l2r"]=a
        a=0.0
        for x in reversed(xs): a+=x
        c["r2l"]=a
        a=Q(0)
        for x in xs: a=E64.rn(a+Q(x))
        c["x87"]=float(B64.rn(a))
        a=Q(0)
        for x in reversed(xs): a=E64.rn(a+Q(x))
        c["x87r"]=float(B64.rn(a))
        c["exact"]=float(B64.rn(sum(map(Q,xs))))
        for k,v in c.items(): tally[k]=tally.get(k,0)+(v==g)
print(tot,tally)
