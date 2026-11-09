import sys, json, struct, numpy as np
sys.path.insert(0,"/home/ubuntu/codex-runs/r42/pylib")
import mpmath as mp
mp.mp.dps=40
def f(b): return struct.unpack("<d",struct.pack("<Q",b))[0]
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/NORM.S.INV.lattice.replies.json"))
qs=np.arange(1.5,40,0.01)
for w in D["windows"]:
    s0=int(w["start"],16); xb=np.array([int(b,16) for b in w["rows"]],dtype=np.int64)
    p0=f(s0); p1=f(s0+511)
    X0=mp.sqrt(2)*mp.erfinv(2*mp.mpf(p0)-1); X1=mp.sqrt(2)*mp.erfinv(2*mp.mpf(p1)-1)
    u=2.0**-52 if float(X0)>=1 else 2.0**-53
    Dsl=float((X1-X0)/511/u)
    # x in units of u
    xv=np.array([f(int(b)) for b in xb])/u
    r=np.diff(xv)-Dsl
    sd=r.std()
    S=np.abs(np.exp(2j*np.pi*np.outer(1/qs,r)).mean(axis=1))
    ok=(qs<2.5*sd)
    if not ok.any(): print("t=%.3f sd=%.2f none"%(w["t"],sd)); continue
    Sm=S[ok].max(); cand=np.nonzero(ok&(S>=0.6*Sm))[0]; i=cand[-1]
    j=max(0,i-30)+np.argmax(S[max(0,i-30):i+30])
    print("t=%.3f sd=%.2f Q=%.2f (units 2^-53: %.2f) S=%.2f"%(w["t"],sd,qs[j],qs[j]*u/2**-53,S[j]))
