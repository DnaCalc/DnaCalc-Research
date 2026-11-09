import json, random, math, struct, sys, collections
sys.path.insert(0,"/home/ubuntu/codex-runs/r42/pylib")
import mpmath as mp
mp.mp.dps=80
f=lambda h: struct.unpack("<d", struct.pack("<Q", int(h,16)))[0]
Q=json.load(open("/home/ubuntu/codex-runs/r43/q/LINEST.thr3-questions.json"))["witnesses"]
ans={}
for ch in json.load(open("/home/ubuntu/codex-runs/r43/fresh/LINEST.thr3-questions.replies.json"))["replies"]:
    for i,row in zip(ch["ids"],(ch["reply"] or {}).get("rows") or []):
        ans[i]=f(row["bits"])
def house_steps(cols,const,order):
    # exact (mp) Householder with pivot order given (list of column indices), ones first if const; returns for each step (x0, s, vv, before, run_rem)
    n=len(cols[0]); M=[]
    if const: M.append([mp.mpf(1)]*n)
    for c in cols:
        c=[mp.mpf(v) for v in c]
        if const: m=sum(c)/n; c=[v-m for v in c]
        M.append(c)
    seq=([0] if const else [])+[k+int(const) for k in order]
    before={k:sum(v*v for v in M[k]) for k in range(len(M))}
    out=[]
    for j,k in enumerate(seq):
        a=M[k]; x0=a[j]; s=sum(v*v for v in a[j:]); rest=s-x0*x0
        alpha=-mp.sqrt(s) if x0>0 else mp.sqrt(s); v0=x0-alpha; vv=rest+v0*v0
        out.append(dict(col=k,x0=x0,s=s,vv=vv,before=before[k]))
        v=[v0]+a[j+1:]
        for kk in seq[j+1:]:
            b=M[kk]; t=sum(vi*bi for vi,bi in zip(v,b[j:])); fct=2*t/vv
            for i in range(j,n): b[i]-=fct*v[i-j]
    return out
by=collections.defaultdict(list)
for w in Q:
    a=ans[w["id"]]; n=len(w["args"][0]); const=w["args"][2]
    X=[[f(v) for v in r] for r in w["args"][1]]; cols=[[X[i][j] for i in range(n)] for j in range(len(X[0]))]
    st=house_steps(cols[:2],const,[1,0])   # z then x
    sx=st[-1]
    by[w["tag"]].append((w["log2ratio"], n-int(const)-int(a), float(mp.log(sx["s"]/sx["before"],2)), float(mp.log(sx["vv"]/sx["before"],2))))
for t,v in by.items():
    if not t.endswith("xz0") and not t.endswith("xz"): continue
    v.sort(); lo=max(r for r in v if r[1]==min(q[1] for q in v)); hi=min(r for r in v if r[1]!=min(q[1] for q in v))
    print("%-12s drop: s %.4f vv %.4f | keep: s %.4f vv %.4f"%(t,lo[2],lo[3],hi[2],hi[3]))
