import json, struct, math, collections
from racer import invoke
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.a1.json"))
f=lambda h: struct.unpack("<d",struct.pack("<Q",int(h,16)))[0]
Q=D["questions"]; E=[f(r["bits"]) for r in D["answers"]]
O=invoke([("BETA.DIST",[x,a,b,"T"]) for x,a,b in Q])
c=collections.Counter()
for (x,a,b),e,o in zip(Q,E,O):
    reg="bpser" if x*b<=0.7 else "other"
    c[(b,reg,"n")]+=1; c[(b,reg,"ours")]+=o==e; c[(b,reg,"expm1log1p")]+=(-math.expm1(b*math.log1p(-x)))==e
for b in sorted(set(q[2] for q in Q)):
    print(b, " ".join("%s: n=%d ours=%d em=%d"%(r,c[(b,r,"n")],c[(b,r,"ours")],c[(b,r,"expm1log1p")]) for r in ["bpser","other"]))
