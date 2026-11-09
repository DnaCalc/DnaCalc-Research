import json, struct, collections, math
from fractions import Fraction as Fr
f=lambda h: struct.unpack("<d", struct.pack("<Q", int(h,16)))[0]
Q=json.load(open("/home/ubuntu/codex-runs/r43/q/LINEST.thr2-questions.json"))["witnesses"]
ans={}
for ch in json.load(open("/home/ubuntu/codex-runs/r43/fresh/LINEST.thr2-questions.replies.json"))["replies"]:
    for i,row in zip(ch["ids"],(ch["reply"] or {}).get("rows") or []):
        ans[i]=f(row["bits"]) if row.get("type")=="number" else None
def resid_ratio(cols,xi,const):
    # remainder of column xi after projecting out the other nonzero columns (and ones), / its own centred ss
    n=len(cols[0]); C=[[Fr(v) for v in c] for c in cols]
    if const: C=[[v-sum(c)/n for v in c] for c in C]
    basis=[]
    for k,c in enumerate(C):
        if k==xi or all(v==0 for v in c): continue
        v=c[:]
        for b in basis:
            pr=sum(a*bb for a,bb in zip(v,b))/sum(bb*bb for bb in b); v=[a-pr*bb for a,bb in zip(v,b)]
        basis.append(v)
    x=C[xi][:]; sx=sum(v*v for v in x)
    for b in basis:
        pr=sum(a*bb for a,bb in zip(x,b))/sum(bb*bb for bb in b); x=[a-pr*bb for a,bb in zip(x,b)]
    return sum(v*v for v in x)/sx
by=collections.defaultdict(list)
for w in Q:
    if not w["tag"].split("-")[1].startswith("xzw"): continue
    a=ans.get(w["id"]); const=w["args"][2]
    X=[[f(v) for v in r] for r in w["args"][1]]; cols=[[X[i][j] for i in range(6)] for j in range(len(X[0]))]
    r=math.log2(resid_ratio(cols,0,const))
    by[w["tag"]].append((r, 6-int(const)-int(a)))
for t,v in by.items():
    v.sort(); cs=[c for r,c in v]
    lo=max([r for r,c in v if c==min(cs)]); hi=min([r for r,c in v if c!=min(cs)])
    print("%-12s drop<=%.3f keep>=%.3f mono=%s"%(t,lo,hi,all(cs[i]<=cs[i+1] for i in range(len(cs)-1))))
