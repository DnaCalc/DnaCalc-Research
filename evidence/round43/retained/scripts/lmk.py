# build inputs files with expected for fresh R42 sets
import json
QD="/home/ubuntu/codex-runs/r42/questions"; RD="/home/ubuntu/codex-runs/claude-r42/fresh"
for qf in ["LINEST.questions.json","LOGEST.questions.json","LINEST.sign-questions.json","LINEST.zero-questions.json","LINEST.pivot-questions.json","LINEST.pivot3-questions.json","LOGEST.pivot3-questions.json"]:
    ans={}
    for ch in json.load(open(RD+"/"+qf.replace(".json",".replies.json")))["replies"]:
        for i,row in zip(ch["ids"],(ch["reply"] or {}).get("rows") or []):
            if row.get("type")=="number": ans[i]="0x"+row["bits"]
            elif row.get("type")=="error": ans[i]="error:"+str(row.get("code") or row.get("value"))
    Q=json.load(open(QD+"/"+qf))
    ws=[]
    for w in Q["witnesses"]:
        w=dict(w); w["expected"]=ans.get(w["id"],"missing"); ws.append(w)
    Q["witnesses"]=ws
    name="F42-"+qf.replace("-questions.json","").replace(".questions.json","").replace(".","-")
    json.dump(Q,open("/home/ubuntu/codex-runs/r43/lin/%s.inputs.json"%name,"w"))
    print(name,len(ws),sum(w["expected"]!="missing" for w in ws))
