import json, struct, math, collections, numpy as np
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.a1.json"))
f=lambda h: struct.unpack("<d",struct.pack("<Q",int(h,16)))[0]
Q=D["questions"]; E=[f(r["bits"]) for r in D["answers"]]
c=collections.Counter()
def bits(v): return struct.unpack("<q",struct.pack("<d",v))[0]
for (x,a,b),e in zip(Q,E):
    dec=int(math.floor(math.log10(x)))
    w=-math.expm1(b*math.log1p(-x))
    c[(dec,"n")]+=1; c[(dec,"ok")]+=w==e; c[(dec,"d")]+=abs(bits(w)-bits(e))
for d in range(-9,-1): print(d, "%d/%d"%(c[(d,"ok")],c[(d,"n")]), "mean|ulp| %.2f"%(c[(d,"d")]/max(1,c[(d,"n")])))
