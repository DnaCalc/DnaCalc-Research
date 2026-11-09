import json, sys, os, collections, struct, math
from lnorm import norm
f=lambda h: struct.unpack("<d", struct.pack("<Q", int(h,16)))[0]
D="/home/ubuntu/codex-runs/r43/lin"; A=sys.argv[1]
for fn in sys.argv[2:]:
    W=json.load(open(os.path.join(D,fn+".inputs.json")))["witnesses"]
    pa={json.loads(l)["id"]:json.loads(l)["calc_func_prediction"] for l in open(os.path.join(D,fn+".%s.jsonl"%A))}
    c=collections.Counter()
    for w in W:
        e=norm(w["expected"]); p=pa[w["id"]]
        if e==p: continue
        a=w["args"]
        try:
            X=[[f(v) for v in r] for r in a[1]]; n=len(a[0]); k=len(X[0]) if len(X)==n else -1
            cols=[[X[i][j] for i in range(n)] for j in range(k)] if k>0 else []
            z=any(all(v==0 for v in cc) for cc in cols)
        except Exception: n=k=-1; z=None
        if e=="missing": kind="missing"
        elif e.startswith("0x") and p.startswith("0x"):
            ve,vp=f(e),f(p)
            if ve==0 or vp==0: kind="zero-vs-nonzero"
            else:
                u=abs(int(e,16)-int(p,16)); kind="ulp<=4" if u<=4 else ("ulp<=1e3" if u<1000 else "big")
        else: kind="type:"+e[:12]+"/"+p[:12]
        c[(kind,str(w["result_index"]),"n=%d"%n,"k=%d"%k,"const=%s"%a[2] if len(a)>2 else "","zero=%s"%z)]+=1
    print(fn)
    for k,v in c.most_common(40): print("  ",v,k)
