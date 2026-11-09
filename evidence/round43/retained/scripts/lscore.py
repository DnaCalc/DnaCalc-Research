import json, sys, subprocess, os, collections
from lnorm import norm
from concurrent.futures import ThreadPoolExecutor
racer, label = sys.argv[1], sys.argv[2]
D="/home/ubuntu/codex-runs/r43/lin"
files=sorted(f for f in os.listdir(D) if f.endswith(".inputs.json"))
def run(f):
    out=os.path.join(D,f.replace(".inputs.json",".%s.jsonl"%label))
    if not os.path.exists(out):
        subprocess.run([racer,"invoke-inputs",os.path.join(D,f),out],check=True,capture_output=True)
    return f,out
with ThreadPoolExecutor(8) as ex: res=list(ex.map(run,files))
tot=collections.Counter()
for f,out in res:
    P={json.loads(l)["id"]:json.loads(l)["calc_func_prediction"] for l in open(out)}
    W=json.load(open(os.path.join(D,f)))["witnesses"]
    ok=sum(P.get(w["id"])==norm(w["expected"]) for w in W)
    print("%-28s %6d / %6d"%(f.replace(".inputs.json",""),ok,len(W)))
    tot["ok"]+=ok; tot["n"]+=len(W)
print("TOTAL",tot["ok"],tot["n"])
