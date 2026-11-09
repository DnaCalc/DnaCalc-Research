import json, struct, math, numpy as np
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.map.json"))
f=lambda h: struct.unpack("<d",struct.pack("<Q",int(h,16)))[0]
LD=np.longdouble
cands={
 "1-pow(1-x,b)": lambda x,a,b: 1-(1-x)**b,
 "-expm1(b*log1p(-x))": lambda x,a,b: -math.expm1(b*math.log1p(-x)),
 "1-exp(b*log(1-x))": lambda x,a,b: 1-math.exp(b*math.log(1-x)),
 "-expm1(b*log(1-x))": lambda x,a,b: -math.expm1(b*math.log(1-x)),
}
for side,(ia) in [("a=1",1)]:
    rows=[(x,a,b,f(ans["bits"]),o) for (x,a,b),ans,o in zip(D["questions"],D["answers"],D["ours"]) if a==1.0 and b!=1.0 and ans.get("type")=="number"]
    print(len(rows),"ours exact",sum(int(o,16)==struct.unpack("<Q",struct.pack("<d",e))[0] for x,a,b,e,o in rows))
    for k,fn in cands.items():
        print(k, sum(fn(x,a,b)==e for x,a,b,e,o in rows))
rows=[(x,a,b,f(ans["bits"])) for (x,a,b),ans,o in zip(D["questions"],D["answers"],D["ours"]) if b==1.0 and a!=1.0 and ans.get("type")=="number"]
print("b=1",len(rows),"pow",sum(x**a==e for x,a,b,e in rows),"exp(a log x)",sum(math.exp(a*math.log(x))==e for x,a,b,e in rows))
