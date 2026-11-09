import json, sys, os, collections
from lnorm import norm

D="/home/ubuntu/codex-runs/r43/lin"; A,B=sys.argv[1:3]
for f in sorted(os.listdir(D)):
    if not f.endswith(".inputs.json"): continue
    W=json.load(open(os.path.join(D,f)))["witnesses"]
    pa={json.loads(l)["id"]:json.loads(l)["calc_func_prediction"] for l in open(os.path.join(D,f.replace(".inputs.json",".%s.jsonl"%A)))}
    pb={json.loads(l)["id"]:json.loads(l)["calc_func_prediction"] for l in open(os.path.join(D,f.replace(".inputs.json",".%s.jsonl"%B)))}
    g=sum(pa[w["id"]]!=norm(w["expected"]) and pb[w["id"]]==norm(w["expected"]) for w in W)
    l=[w for w in W if pa[w["id"]]==norm(w["expected"]) and pb[w["id"]]!=norm(w["expected"])]
    print("%-24s +%d -%d"%(f[:-12],g,len(l)), [ (w["id"],w["result_index"]) for w in l[:4]])
