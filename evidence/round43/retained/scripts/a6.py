from nl import *
import sys
d=load(); t=d["t"]; x=d["x"]; src=d["src"]
s=int(sys.argv[1]); Q=float(sys.argv[2])
m=src==s
X=np.sqrt(LD(2))*(LD(d["hi"][m])+LD(d["lo"][m])); u=ulp(x[m])
r=np.asarray((LD(x[m])-X)/LD(u),dtype=F)
# phase: find phase maximizing
ph=np.angle(np.mean(np.exp(2j*np.pi*r/Q)))*Q/(2*np.pi)
k=np.round((r-ph)/Q); res=r-ph-k*Q
print("phase %.3f"%ph, "res sd %.3f"%res.std(), "k hist", np.unique(k,return_counts=True))
# does phase drift across window?
for a in range(0,len(r),1024):
    rr=r[a:a+1024]; print(a, "%.3f"%(np.angle(np.mean(np.exp(2j*np.pi*rr/Q)))*Q/(2*np.pi)), "%.3f"%abs(np.mean(np.exp(2j*np.pi*rr/Q))))
