from nl import *
d=load(); t=d["t"]; x=d["x"]; p=d["p"]; src=d["src"]
P=preimages(x)
y=np.where(np.isnan(P[:,2]), np.nanmean(P,axis=1), P[:,2])
e=((y-d["hi"])-d["lo"])/ulp(d["hi"])
# x error
X=np.sqrt(2)*0  # skip
for s in np.unique(src):
    m=src==s; tq=t[m]
    print("%5d n=%5d t=[%.6f,%.6f] ey mean %+.3f sd %.3f"%(s,m.sum(),tq.min(),tq.max(),e[m].mean(),e[m].std()))
