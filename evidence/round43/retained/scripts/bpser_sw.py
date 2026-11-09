import json, os, struct, collections, sys
sys.argv=["x"]
exec(open("beta_branch.py").read().split("D=json.load")[0])
import racer
racer.RACER=os.path.expanduser("~/codex-runs/r43/racer-bpser")
D=json.load(open("fresh/BETA.DIST.map.json"))
cells=[(q,a["bits"]) for q,a in zip(D["questions"],D["answers"]) if a.get("type")=="number" and branch(*q).endswith("bpser")]
h=lambda v: "%016x"%struct.unpack("<Q",struct.pack("<d",v))[0] if not isinstance(v,str) else v
for sw in range(16):
    os.environ["R43_BPSER"]=str(sw)
    O=racer.invoke([("BETA.DIST",[x,a,b,"T"]) for (x,a,b),e in cells])
    c=collections.Counter()
    for ((x,a,b),e),o in zip(cells,O): c[branch(x,a,b)]+=h(o)==e
    print(sw, dict(c), len(cells))
