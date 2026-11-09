import json, numpy as np, collections
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.map.json"))
Q=D["questions"]; A=D["answers"]; O=D["ours"]
AB=sorted(set(a for x,a,b in Q))
tab=collections.defaultdict(lambda:[0,0,0])
for (x,a,b),ans,o in zip(Q,A,O):
    if ans.get("type")!="number": continue
    e=int(ans["bits"],16); u=abs(e-int(o,16)) if len(o)==16 else 10**9
    t=tab[(a,b)]; t[0]+=1; t[1]+=u==0; t[2]=max(t[2],u)
print("rows a, cols b: exact fraction (%)")
print("      "+" ".join("%5g"%b for b in AB))
for a in AB: print("%5g "%a+" ".join("%5d"%(100*tab[(a,b)][1]/max(1,tab[(a,b)][0])) for b in AB))
print("max ulp")
for a in AB: print("%5g "%a+" ".join("%5s"%(tab[(a,b)][2] if tab[(a,b)][2]<10**5 else "big") for b in AB))
