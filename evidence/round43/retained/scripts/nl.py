import numpy as np, struct, math
LD=np.longdouble; F=np.float64
def f(b): return struct.unpack("<d", struct.pack("<Q", int(b)))[0]
def bits(x): return struct.unpack("<Q", struct.pack("<d", float(x)))[0]
S2=F(math.sqrt(2.0)); K=F(0.88622692545275801365)  # RN(sqrt(pi)/2)
def st(v): return np.asarray(v).astype(F)
def mul(a,b): return st(LD(a)*LD(b)) if np.isscalar(a) and np.isscalar(b) else st(np.asarray(a,dtype=LD)*np.asarray(b,dtype=LD))
def add(a,b): return st(np.asarray(a,dtype=LD)+np.asarray(b,dtype=LD))
def div(a,b): return st(np.asarray(a,dtype=LD)/np.asarray(b,dtype=LD))
def ulp(x): x=np.abs(np.asarray(x,dtype=F)); return np.spacing(x)
def load():
    D=np.load("/home/ubuntu/codex-runs/r43/normal-all.npz"); T=np.load("/home/ubuntu/codex-runs/r43/central-truth.npz")
    idx=T["idx"]; p=D["p"].view(F)[idx]; x=D["x"].view(F)[idx]; src=D["src"][idx]
    return dict(p=p,x=x,src=src,t=T["t"],hi=T["hi"],lo=T["lo"])
def preimages(x):
    """all doubles y with st(LD(y)*LD(S2))==x, as list per point (vectorised: returns array (n,5) with nan where not)"""
    y0=st(LD(x)/LD(S2)); out=np.full((len(x),5),np.nan)
    for k,d in enumerate(range(-2,3)):
        y=y0.copy()
        if d>0:
            for _ in range(d): y=np.nextafter(y,np.inf)
        elif d<0:
            for _ in range(-d): y=np.nextafter(y,-np.inf)
        ok=mul(y,S2)==x
        out[ok,k]=y[ok]
    return out
