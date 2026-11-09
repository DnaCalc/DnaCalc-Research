import pickle, sys, collections, json
sys.path.insert(0,"/home/ubuntu/codex-runs/r43")
from qr2 import ls2
cases=pickle.load(open("/home/ubuntu/codex-runs/r43/dfcases.pkl","rb"))
o=json.loads(sys.argv[1])
def desc(cols):
    out=[]
    for j,c in enumerate(cols):
        if all(v==0 for v in c): out.append("0"); continue
        lab="c%d"%j
        for i in range(j):
            d=cols[i]
            if any(v!=0 for v in d):
                ratios=[a/b for a,b in zip(c,d) if b!=0]
                if ratios and max(ratios)-min(ratios) < 1e-5*abs(ratios[0]) and all((b==0)==(a==0) for a,b in zip(c,d)):
                    r=ratios[0]; exact=max(ratios)==min(ratios)
                    lab="%s%g*c%d"%("" if exact else "~",round(r,4),i); break
        out.append(lab)
    return out
cnt=collections.Counter(); ex={}
for c in cases:
    k=ls2([list(x) for x in c["cols"]],c["y"],c["const"],o)
    if len(k)!=c["count"]:
        key=(c["const"],len(c["y"]),tuple(desc(c["cols"])),"excel",c["count"],"emu",len(k))
        cnt[key]+=1; ex.setdefault(key,c["src"][:14]+":"+c["id"])
for k,v in sorted(cnt.items(),key=lambda t:-t[1]): print(v,k,ex[k])
