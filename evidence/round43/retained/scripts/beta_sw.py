import json, struct, os, subprocess, collections
import racer
racer.RACER=os.path.expanduser("~/codex-runs/r43/racer-beta")
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.a1.json"))
f=lambda h: struct.unpack("<d",struct.pack("<Q",int(h,16)))[0]
Q=D["questions"]; E=[f(r["bits"]) for r in D["answers"]]
for v in range(0,8):
    os.environ["R43_BETA_A1"]=str(v)
    O=racer.invoke([("BETA.DIST",[x,a,b,"T"]) for x,a,b in Q])
    print(v, sum(o==e for o,e in zip(O,E)), len(Q))
