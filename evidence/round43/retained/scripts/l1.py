import json, struct, collections
QD="/home/ubuntu/codex-runs/r42/questions"; RD="/home/ubuntu/codex-runs/claude-r42/fresh"
f=lambda h: struct.unpack("<d", struct.pack("<Q", int(h,16)))[0]
for fn in ["LINEST","LOGEST"]:
    ans={}
    for ch in json.load(open(RD+"/%s.pivot3-questions.replies.json"%fn))["replies"]:
        for i,row in zip(ch["ids"],(ch["reply"] or {}).get("rows") or []):
            ans[i]="0x"+row["bits"] if row.get("type")=="number" else "error:"+str(row.get("code") or row.get("value"))
    W=json.load(open(QD+"/%s.pivot3-questions.json"%fn))["witnesses"]
    for b in ["base","downc","sign"]:
        P={json.loads(l)["id"]:json.loads(l)["calc_func_prediction"] for l in open(QD+"/%s.pivot3-pred-%s.jsonl"%(fn,b))}
        bad=[w for w in W if ans.get(w["id"])!=P.get(w["id"])]
        print(fn,b,len(W)-len(bad),len(W))
        if b=="downc":
            tags=collections.Counter((w["tag"],tuple(w["result_index"])) for w in bad); print(tags)
            json.dump([dict(w,excel=ans.get(w["id"]),pred=P.get(w["id"])) for w in bad],open("/home/ubuntu/codex-runs/r43/%s-bad.json"%fn,"w"))
