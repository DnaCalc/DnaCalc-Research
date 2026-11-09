from nl import *
d=load(); t=d["t"]; x=d["x"]; src=d["src"]
hi=d["hi"]; lo=d["lo"]
yc=hi.copy()  # hi is float(true) = RN(true) approximately (mp float conversion rounds nearest)
for name,m in [("B",(t>0.85)&(t<=0.875)),("inner0",np.abs(t)<0.01),("outer",(np.abs(t)>0.7072)&(np.abs(t)<=0.85))]:
    xm=x[m]; y=yc[m]
    for k in [-2,-1,0,1,2]:
        yy=y.copy()
        for _ in range(abs(k)): yy=np.nextafter(yy, np.inf if k>0 else -np.inf)
        print(name, m.sum(), "y=RN(erfinv)%+d"%k, "match %.4f"%(mul(yy,S2)==xm).mean())
    # direct x=RN(sqrt2*erfinv)
    X=np.sqrt(LD(2))*(LD(hi[m])+LD(lo[m]))
    print(name, "x=RN(true X)", (st(X)==xm).mean())
