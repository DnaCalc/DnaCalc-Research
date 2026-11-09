import json, struct, collections
f=lambda h: struct.unpack("<d", struct.pack("<Q", int(h,16)))[0]
Q=json.load(open("/home/ubuntu/codex-runs/r43/q/LINEST.thr-questions.json"))["witnesses"]
ans={}
for ch in json.load(open("/home/ubuntu/codex-runs/r43/fresh/LINEST.thr-questions.replies.json"))["replies"]:
    for i,row in zip(ch["ids"],(ch["reply"] or {}).get("rows") or []):
        ans[i]=f(row["bits"]) if row.get("type")=="number" else row
by=collections.defaultdict(list)
for w in Q:
    a=ans.get(w["id"]); n=len(w["args"][0]); const=w["args"][2]
    if not isinstance(a,float): by[w["tag"]].append((w["log2ratio"],"?")); continue
    cnt=n-int(const)-int(a); by[w["tag"]].append((w["log2ratio"],cnt))
for t,v in by.items():
    v.sort(); print("%-14s"%t," ".join("%.1f:%s"%(r,c) for r,c in v))
