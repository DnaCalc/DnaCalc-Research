from nl import *
d=load(); t=d["t"]; x=d["x"]; src=d["src"]
rows=[]
for s in np.unique(src):
    m=src==s; tt=t[m][0]
    if not (0.3<abs(tt)<0.875): continue
    X=np.sqrt(LD(2))*(LD(d["hi"][m])+LD(d["lo"][m])); u=ulp(x[m])
    r=np.diff(np.asarray((LD(x[m])-X)/LD(2.0**-53),dtype=F))   # units 2^-53 always
    qs=np.arange(0.8,40,0.002); S=np.abs(np.exp(2j*np.pi*np.outer(1/qs,r)).mean(axis=1))
    # pick the largest q with S > 0.5*max(S[qs<30]) but avoid trivial large q: use q<=r.std()*0.8
    lim=max(2.0, r.std()*0.9); sel=qs<lim
    i=np.argmax(np.where(sel,S,0)); rows.append((abs(tt),qs[i],S[i],r.std(),m.sum()))
for a in sorted(rows): print("t=%.5f Q=%.3f S=%.2f dsd=%.2f n=%d"%a)
