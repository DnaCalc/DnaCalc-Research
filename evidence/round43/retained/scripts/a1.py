from nl import *
d=load(); t=d["t"]; x=d["x"]; p=d["p"]
assert np.all(t==2*(p-0.5))
P=preimages(x); n=np.sum(~np.isnan(P),axis=1)
print("preimage counts", np.bincount(n))
y=np.nanmin(np.abs(P),axis=1)*np.sign(x)  # choose one
# simple formula agreement by |q| band
yk=mul(t,K); xk=mul(yk,S2)
aq=np.abs(p-0.5)
bands=np.arange(0,0.4251,0.025)
for a,b in zip(bands[:-1],bands[1:]):
    m=(aq>=a)&(aq<b)
    if m.sum(): print("%.3f-%.3f n=%7d simple-match %.4f"%(a,b,m.sum(),(xk[m]==x[m]).mean()))
