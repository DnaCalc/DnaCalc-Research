import sys, struct, math
sys.path.insert(0,"/home/koderbot/projects/DnaCalc/DnaCalc-Research/racer")
from fractions import Fraction as Q
from decimal import Decimal as D, getcontext
getcontext().prec=50
from intervals import B64
def f(hx): return struct.unpack(">d", bytes.fromhex(hx))[0]
rows=[[f(x) for x in l.split()] for l in open(sys.argv[1])]
L=dict()
for l in open(sys.argv[2]):
    r,lh=l.split(); L[float(r)]=f(lh)
sel=float(sys.argv[3])
for p in rows:
    r,n,pv,k,pm,P=p[0],p[1],p[2],p[3],p[6],p[7]
    if r!=sel: continue
    m=n-k+1
    lo=(Q(P)+B64.prev(Q(P)))/2; hi=(Q(P)+B64.next(Q(P)))/2
    ulo=lo/Q(pm); uhi=hi/Q(pm)
    tl=(D(ulo.numerator)/D(ulo.denominator)).ln(); th=(D(uhi.numerator)/D(uhi.denominator)).ln()
    t0=-(m*L[r])
    ul=D(math.ulp(t0))
    a=(tl-D(t0))/ul; b=(th-D(t0))/ul
    a,b=min(a,b),max(a,b)
    # exact t = -m*ln(1+r)
    R=D(r); te=-(D(m)*(1+R).ln())
    print(int(m), "%.2f %.2f"%(a,b), "exact %.2f"%((te-D(t0))/ul))
