import json, struct, collections
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.map.json"))
c=collections.Counter(); n=collections.Counter(); ex=[]
for (x,a,b),ans,o in zip(D["questions"],D["answers"],D["ours"]):
    if a!=1.0 or b==1.0 or ans.get("type")!="number": continue
    e=int(ans["bits"],16); u=e-int(o,16)
    band="x<0.1" if x<0.1 else ("x>0.9" if x>0.9 else "mid")
    n[(b,band)]+=1; c[(b,band)]+=u!=0
    if u!=0 and len(ex)<12 and b in (2.0,500.0): ex.append((x,b,u))
for b in sorted(set(k[0] for k in n)): print(b, {band:"%d/%d"%(c[(b,band)],n[(b,band)]) for band in ["x<0.1","mid","x>0.9"]})
print(ex)
