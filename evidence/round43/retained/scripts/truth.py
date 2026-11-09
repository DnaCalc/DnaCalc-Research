import sys, numpy as np
sys.path.insert(0, "/home/ubuntu/codex-runs/r42/pylib")
import mpmath as mp
from multiprocessing import Pool
from scipy.special import erfinv
mp.mp.dps = 34
def tv(t):
    # true erfinv(t) for double t via Newton in mpmath from scipy start
    x = mp.mpf(float(erfinv(t))); tt = mp.mpf(t)
    for _ in range(3):
        x = x - (mp.erf(x) - tt) / (2/mp.sqrt(mp.pi)*mp.exp(-x*x))
    return x
def work(ts):
    out=[]
    for t in ts:
        x=tv(float(t)); y=float(x); r=float(x-mp.mpf(y))
        out.append((y,r))
    return out
if __name__=="__main__":
    D=np.load("/home/ubuntu/codex-runs/r43/normal-all.npz")
    p=D["p"].view(np.float64)
    m=(p>=0.075)&(p<=0.9375)
    q=p[m]-0.5; t=2*q
    chunks=np.array_split(t,64)
    with Pool(16) as pool: res=pool.map(work,chunks)
    res=[a for c in res for a in c]
    hi=np.array([a for a,b in res]); lo=np.array([b for a,b in res])
    np.savez_compressed("/home/ubuntu/codex-runs/r43/central-truth.npz", idx=np.nonzero(m)[0], t=t, hi=hi, lo=lo)
    print(len(t))
