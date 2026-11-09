import sys, numpy as np
sys.path.insert(0, "/home/ubuntu/codex-runs/r42/pylib")
import mpmath as mp
mp.mp.dps=60
rows=np.load("/home/ubuntu/codex-runs/r43/bias.npy")
sel=(rows[:,0]>0.3)&(rows[:,0]<0.7071)
T=rows[sel,0]; B=rows[sel,1]; E=np.maximum(rows[sel,2],0.02)
# erfinv maclaurin: c_k recursion; erfinv(z)=sum c_k/(2k+1) (sqrt(pi)/2 z)^(2k+1)
N=60; c=[mp.mpf(1)]
for k in range(1,N):
    c.append(sum(c[m]*c[k-1-m]/((m+1)*(2*m+1)) for m in range(k)))
K=mp.sqrt(mp.pi)/2
a=[c[k]/(2*k+1)*K**(2*k+1) for k in range(N)]   # erfinv(t)=sum a_k t^(2k+1) => F(s)=sum a_k s^k
tv={tt:mp.erfinv(mp.mpf(tt)) for tt in T}
def score(m,n):
    if n==0:
        p=a[:m+1]; q=[1]
    else:
        p,q=mp.pade(a[:m+n+1],m,n)
    err=[]
    for tt in T:
        s=mp.mpf(tt)**2
        v=mp.mpf(tt)*mp.polyval(p[::-1],s)/mp.polyval(q[::-1],s)
        err.append(float((v-tv[tt])/mp.mpf(2)**-53))
    err=np.array(err)
    # fit observed B = lin*t + err  (lin free)
    A=np.vstack([T]).T
    r=B-err; lin=np.sum(r*T/E**2)/np.sum(T*T/E**2)
    chi=np.sum(((r-lin*T)/E)**2)
    return chi,lin,err
res=[]
for tot in []:
    for n in range(0,tot+1):
        m=tot-n
        try: chi,lin,err=score(m,n)
        except Exception as ex: continue
        res.append((chi,m,n,lin,err[-1]))
res.sort()
for r in res[:25]: print("chi %.1f [%d/%d] lin %.3f err_end %.2f"%r)
print("dof", len(T))
print("t obs-0.89t | taylor N=6..12 errs")
errs={}
for m in [5,6,7,8,9,10,11,12,13,14]:
    errs[m]=score(m,0)[2]
for i,tt in enumerate(T):
    print("%.4f %+6.2f | "%(tt,B[i]-0.89*tt)+" ".join("%+7.2f"%errs[m][i] for m in errs))
