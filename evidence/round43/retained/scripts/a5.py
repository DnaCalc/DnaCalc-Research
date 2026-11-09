from nl import *
import sys
d=load(); t=d["t"]; x=d["x"]; src=d["src"]
def xtrue(m):  # S2true*erfinv in units: compute true X = sqrt2*(hi+lo) as float + residual
    hi=d["hi"][m]; lo=d["lo"][m]
    X=np.sqrt(LD(2))*(LD(hi)+LD(lo))
    return X
for s in [int(a) for a in sys.argv[1:]]:
    m=src==s
    X=xtrue(m); u=ulp(x[m])
    r=np.asarray((LD(x[m])-X)/LD(u),dtype=F); r=np.diff(r)
    qs=np.arange(2.5,16,0.01); S=np.array([abs(np.mean(np.exp(2j*np.pi*r/q))) for q in qs])
    o=np.argsort(-S)[:6]
    print(s, "t=%.4f"%t[m][0], "u=%g"%u[0], " ".join("%.2f:%.2f"%(qs[i],S[i]) for i in o))
