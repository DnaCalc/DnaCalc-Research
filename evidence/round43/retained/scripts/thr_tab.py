import json, struct, collections, sys
f=lambda h: struct.unpack("<d", struct.pack("<Q", int(h,16)))[0]
name=sys.argv[1]
Q=json.load(open("/home/ubuntu/codex-runs/r43/q/LINEST.%s-questions.json"%name))["witnesses"]
ans={}
for ch in json.load(open("/home/ubuntu/codex-runs/r43/fresh/LINEST.%s-questions.replies.json"%name))["replies"]:
    for i,row in zip(ch["ids"],(ch["reply"] or {}).get("rows") or []):
        ans[i]=f(row["bits"]) if row.get("type")=="number" else None
by=collections.defaultdict(list)
for w in Q:
    a=ans.get(w["id"]); n=len(w["args"][0]); const=w["args"][2]
    by[w["tag"]].append((w["log2ratio"], None if a is None else n-int(const)-int(a)))
for t,v in by.items():
    v.sort(); cs=[c for r,c in v if c is not None]
    lo=max([r for r,c in v if c==min(cs)],default=None); hi=min([r for r,c in v if c!=min(cs)],default=None)
    mono=all(cs[i]<=cs[i+1] for i in range(len(cs)-1))
    print("%-16s drop<=%s keep>=%s mono=%s counts %s"%(t,"%.4f"%lo if lo is not None else None,"%.4f"%hi if hi is not None else None,mono,sorted(set(cs))))
