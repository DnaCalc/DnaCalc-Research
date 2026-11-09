import sys, struct
sys.path.insert(0,"/home/koderbot/projects/DnaCalc/DnaCalc-Research/racer")
from fractions import Fraction as Q
from intervals import B64, E64
def f(hx): return struct.unpack(">d", bytes.fromhex(hx))[0]
def pre(P): return (P+B64.prev(P))/2, (P+B64.next(P))/2
def cover(xs,Ps):
    ev=[]
    for x,P in zip(xs,Ps):
        a,b=pre(P); c1,c2=sorted((a/x,b/x)); ev.append((c1,0));ev.append((c2,1))
    ev.sort(); cur=best=0; at=None
    for x,t in ev:
        if t==0:
            cur+=1
            if cur>best: best=cur; at=x
        else: cur-=1
    return best, at
rows=[l.split() for l in open(sys.argv[1])][:int(sys.argv[2])]
def parse(v):
    if v.startswith("e64:"): return E64.rn(Q(v[4:]))
    if "/" in v: return Q(v)
    return Q(f(v))
C={k:parse(v) for k,v in (x.split("=") for x in sys.argv[3:])}
pvs=[Q(f(r[0])) for r in rows]; Ps=[Q(f(r[2])) for r in rows]
res=[]
for a,A in C.items():
    for op in "*/\\":
        for prec in (53,64):
            Rr=B64.rn if prec==53 else E64.rn
            ys=[Rr(pv*A) if op=="*" else Rr(pv/A) if op=="/" else Rr(A/pv) for pv in pvs]
            b,at=cover(ys,Ps); res.append((b,a,op,prec,float(at)))
res.sort(reverse=True)
for x in res[:8]: print(x)
