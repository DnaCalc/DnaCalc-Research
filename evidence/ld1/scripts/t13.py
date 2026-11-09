from ask import *
import struct, math
from fractions import Fraction as Q
def nxt(x,i): return struct.unpack(">d",struct.pack(">Q",struct.unpack(">Q",struct.pack(">d",x))[0]+i))[0]
pvs=[nxt(3.3,i) for i in range(400)]
for r in [0.25,0.5,0.625,0.75,0.875,1.0,1.5,3.0,0.1,0.3,2.0]:
    n=3000
    P=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv,n,n,0] for pv in pvs])["rows"]]
    K=-Q(r)/(1+Q(r))
    res=[float((Q(p)-K*Q(pv))/(K*Q(pv))) for p,pv in zip(P,pvs)]
    m=sum(res)/len(res)
    print(r, "mean rel %.3e = %.3f*2^-52"%(m, m/2**-52), "min %.2f max %.2f"%(min(res)/2**-52, max(res)/2**-52))
