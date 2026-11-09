from ask import *
import json
o=json.load(open("oy1.json"))[:40]
cands={"rate":lambda a:a[4],"0.1":lambda a:0.1,"0.05":lambda a:0.05,"0.01":lambda a:0.01,"0":lambda a:0.0,"rate/2":lambda a:a[4]/2,"1":lambda a:1.0,"0.5":lambda a:0.5}
for name,g in cands.items():
    rows=[]; gs=[]
    for x in o:
        a=x["args"]; y=g(a)
        if y<=0: y=0.0
        gs.append(y); rows.append([a[0],a[1],a[2],a[3],a[4],y,a[6],a[7],a[8]])
    p=ask("ODDFPRICE",rows)["rows"]
    q=[]
    for x,pp,y in zip(o,p,gs):
        a=x["args"]; q.append([a[0],a[1],a[2],a[3],a[4],f(pp["bits"]),a[6],a[7],a[8]])
    yy=ask("ODDFYIELD",q)["rows"]
    hits=sum(1 for z,y in zip(yy,gs) if z.get("bits")==h(y))
    near=sum(1 for z,y in zip(yy,gs) if "bits" in z and abs(f(z["bits"])-y)<1e-9)
    print(name,"exact",hits,"near",near,"of",len(o))
