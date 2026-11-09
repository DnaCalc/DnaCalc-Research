import json, struct, sys
sys.path.insert(0,"/home/ubuntu/codex-runs/r43")
from qr import ls
f=lambda h: struct.unpack("<d", struct.pack("<Q", int(h,16)))[0]
B=json.load(open("/home/ubuntu/codex-runs/r43/LINEST-bad.json"))
seen=set()
for w in B:
    key=w["id"].rsplit(":",2)[0]
    if key in seen: continue
    seen.add(key)
    a=w["args"]; X=[[f(v) for v in r] for r in a[1]]; y=[f(r[0]) for r in a[0]]
    cols=[[X[i][j] for i in range(len(X))] for j in range(len(X[0]))]
    kp,log=ls(cols,y,a[2],{},trace=True)
    print(key,"const",a[2],"kept",kp)
    for l in log: print("   ",l)
