from nl import *
d=load(); t=d["t"]; x=d["x"]; src=d["src"]
rows=[]
for s in np.unique(src):
    m=src==s; tt=t[m]
    X=np.sqrt(LD(2))*(LD(d["hi"][m])+LD(d["lo"][m]))
    ea=np.asarray((LD(x[m])-X)/np.sqrt(LD(2)),dtype=F)/2.0**-53   # abs y error in units 2^-53
    rows.append((tt.mean(),ea.mean(),ea.std()/np.sqrt(m.sum()),m.sum()))
rows=np.array(sorted(rows))
np.save("/home/ubuntu/codex-runs/r43/bias.npy",rows)
for r in rows: print("%.5f %+.3f %.3f %d"%tuple(r))
