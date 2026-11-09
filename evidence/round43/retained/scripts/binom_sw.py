import json, os, struct, collections
import racer
racer.RACER=os.path.expanduser("~/codex-runs/r43/racer-binom")
d=json.load(open("fresh/BINOM.vs-BETA.map.json")); Q=d["questions"]; A=[r.get("bits") for r in d["binom"]]
h=lambda v: "%016x"%struct.unpack("<Q",struct.pack("<d",v))[0] if not isinstance(v,str) else v
for v in range(0,5):
    os.environ["R43_BINOM"]=str(v)
    O=racer.invoke([("BINOM.DIST",[k,n,p,"T"]) for k,n,p in Q])
    c=collections.Counter()
    for (k,n,p),a,o in zip(Q,A,O):
        side="p<0.5" if p<0.5 else "p>=0.5"; c[side+":n"]+=1; c[side+":ok"]+=h(o)==a
    print(v, dict(c))
