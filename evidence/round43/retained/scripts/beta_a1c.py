import json, struct, math, numpy as np
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.map.json"))
f=lambda h: struct.unpack("<d",struct.pack("<Q",int(h,16)))[0]
LD=np.longdouble
rows=[(x,b,f(ans["bits"]),f(o)) for (x,a,b),ans,o in zip(D["questions"],D["answers"],D["ours"]) if a==1.0 and b!=1.0 and x<0.1 and ans.get("type")=="number"]
def st(v): return float(v)
C={
 "-expm1(b*log1p(-x))": lambda x,b: -math.expm1(b*math.log1p(-x)),
 "-expm1(b*log(1-x))": lambda x,b: -math.expm1(b*math.log(1-x)),
 "1-exp(b*log1p(-x))": lambda x,b: 1-math.exp(b*math.log1p(-x)),
 "1-(1-x)**b": lambda x,b: 1-(1-x)**b,
 "ext -expm1(b*log1p(-x))": lambda x,b: -math.expm1(st(LD(b)*LD(math.log1p(-x)))),
 "series bpser-like": None,
}
for k,fn in C.items():
    if fn is None: continue
    print("%-28s %d/%d"%(k,sum(fn(x,b)==e for x,b,e,o in rows),len(rows)))
print("ours", sum(o==e for x,b,e,o in rows))
