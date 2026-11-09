import json, struct, os, pickle
f=lambda h: struct.unpack("<d", struct.pack("<Q", int(h,16)))[0]
D="/home/ubuntu/codex-runs/r43/lin"
cases=[]
for fn in sorted(os.listdir(D)):
    if not fn.endswith(".inputs.json") or "LOGEST" in fn or "R38G" in fn: continue
    Q=json.load(open(os.path.join(D,fn)))
    for w in Q["witnesses"]:
        if w["result_index"]!=[4,2] or not w["expected"].startswith("0x"): continue
        a=w["args"]
        if len(a)<4: continue
        X=[[f(v) for v in r] for r in a[1]]; y=[f(r[0]) for r in a[0]]
        n=len(y); k=len(X[0]); const=a[2]
        if len(X)!=n: continue
        cols=[[X[i][j] for i in range(n)] for j in range(k)]
        df=f(w["expected"])
        cases.append(dict(src=fn,id=w["id"],cols=cols,y=y,const=const,count=n-int(const)-int(df)))
pickle.dump(cases,open("/home/ubuntu/codex-runs/r43/dfcases.pkl","wb"))
import collections
print(len(cases), collections.Counter(c["src"] for c in cases))
