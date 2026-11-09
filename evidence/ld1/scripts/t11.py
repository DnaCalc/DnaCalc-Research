from ask import *
import sys
sys.path.insert(0,"/home/koderbot/projects/DnaCalc/DnaCalc-Research/racer")
from fractions import Fraction as Q
from intervals import B64,E64
def x87sub(a,b): return float(B64.rn(E64.rn(Q(a)-Q(b))))
for fn in ["d2.txt","d4.txt","pvscan-big1.txt"]:
    rows=[[f(x) for x in l.split()] for l in open(fn)]
    if fn.startswith("pvscan"):
        rows=[[0.75,1000,p[0],1000,1000,0,p[1],p[2]] for p in rows[:500]]
    ip=[f(x["bits"]) for x in ask("IPMT",[[p[0],p[3],p[1],p[2],0,0] for p in rows])["rows"]]
    tot=0; ok=0; byk={}
    for p,i in zip(rows,ip):
        pred=x87sub(p[6],i); good= pred==p[7]
        m=int(p[1]-p[3]+1)
        b=min(m,6)
        byk.setdefault(b,[0,0]); byk[b][0]+=good; byk[b][1]+=1
        ok+=good; tot+=1
    print(fn, ok, tot, byk)
