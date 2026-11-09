import json, struct
f=lambda h: struct.unpack(">d",bytes.fromhex(h))[0]
d=json.load(open("oz1.json")); T=24.61111111111111
prs=d["prs"]; ys=[f(y) for y in d["ys"]]
def P(Y): return 100*(1+Y/2)**(-T)
def dP(Y): return -T/2*100*(1+Y/2)**(-T-1)
def newton(pr,y0,stop,mode,deriv):
    y=y0
    for it in range(100):
        if not (y>-2): return float("nan"),it
        fy=P(y)-pr
        if isinstance(fy,complex): return float("nan"),it
        if deriv=="an": d=dP(y)
        else:
            h=deriv; d=(P(y+h)-P(y))/h
        step=fy/d
        if mode=="before":
            if abs(step)<stop: return y,it
            y-=step
        else:
            y-=step
            if abs(step)<stop: return y,it
    return y,100
best=[]
for y0 in [0.0,0.01,0.05,0.1,0.2,0.5,1.0]:
  for mode in ["before","after"]:
    for deriv in ["an",1e-5,1e-6,1e-8]:
      for stop in [1e-10,1e-9,1e-11,1e-12]:
        err=0; 
        for pr,Y in zip(prs,ys):
            yy,it=newton(pr,y0,stop,mode,deriv)
            err+=abs(yy-Y) if yy==yy else 1.0
        best.append((err,y0,mode,deriv,stop))
best.sort(key=lambda b:b[0])
for b in best[:8]: print(b)
y0=best[0][1]
for pr,Y in list(zip(prs,ys))[:12]:
    yy,it=newton(pr,best[0][1],best[0][4],best[0][2],best[0][3]); print(pr, Y-yy, it)
print("scan y0")
res=[]
import math
for i in range(0,400):
    y0=i*0.0005
    err=0
    for pr,Y in zip(prs,ys):
        yy,it=newton(pr,y0,1e-10,"before","an")
        err+=abs(yy-Y) if yy==yy else 1.0
    res.append((err,y0))
res.sort()
print(res[:5])
