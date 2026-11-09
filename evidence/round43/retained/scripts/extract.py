import json, numpy as np
W=json.load(open("/home/ubuntu/codex-runs/r42/work/normal-dense-frozen.json"))["windows"]
A=np.load("/home/ubuntu/codex-runs/r42/work/normal-dense-outputs.npz")["answers"]
V=np.load("/home/ubuntu/codex-runs/r42/work/normal-dense-validity.npz")["numeric"]
P=[];X=[];S=[]
for i,w in enumerate(W):
    s0=int(w["start"],16); n=w["count"]
    pb=np.arange(s0,s0+n,dtype=np.uint64)
    m=V[i][:n].astype(bool)
    P.append(pb[m]); X.append(A[i][:n][m]); S.append(np.full(m.sum(),i))
d=json.load(open("/home/ubuntu/codex-runs/claude-r42/fresh/NORM.S.INV.mapping.replies.json"))
for k,r in enumerate(d["replies"]):
    s0=int(r["window"]["start"],16); rows=r["reply"]["rows"]
    pb=[];xb=[]
    for j,row in enumerate(rows):
        if row.get("type")=="number": pb.append(s0+j); xb.append(int(row["bits"],16))
    P.append(np.array(pb,dtype=np.uint64)); X.append(np.array(xb,dtype=np.uint64)); S.append(np.full(len(pb),1000+k))
P=np.concatenate(P); X=np.concatenate(X); S=np.concatenate(S)
p=P.view(np.float64)
print(len(P), (p<0.075).sum(), ((p>=0.075)&(p<=0.925)).sum(), ((p>0.925)&(p<=0.9375)).sum(), (p>0.9375).sum())
np.savez_compressed("/home/ubuntu/codex-runs/r43/normal-all.npz", p=P, x=X, src=S)
