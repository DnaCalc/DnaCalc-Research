import json, sys, math, struct, numpy as np
sys.path.insert(0,"/home/ubuntu/codex-runs/r42/pylib")
import mpmath as mp
mp.mp.dps=40
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.map.json"))
f=lambda h: struct.unpack("<d",struct.pack("<Q",int(h,16)))[0]
def betacf(a,b,x):
    MAXIT=300; EPS=3e-16; FPMIN=1e-300
    qab=a+b; qap=a+1; qam=a-1; c=1.0; d=1-qab*x/qap
    if abs(d)<FPMIN: d=FPMIN
    d=1/d; h=d
    for m in range(1,MAXIT+1):
        m2=2*m; aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1+aa*d; d=FPMIN if abs(d)<FPMIN else d; c=1+aa/c; c=FPMIN if abs(c)<FPMIN else c
        d=1/d; h*=d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1+aa*d; d=FPMIN if abs(d)<FPMIN else d; c=1+aa/c; c=FPMIN if abs(c)<FPMIN else c
        d=1/d; de=d*c; h*=de
        if abs(de-1)<EPS: break
    return h
def hnr(x,a,b):
    bt=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log(1-x))
    if x<(a+1)/(a+b+2): return bt*betacf(a,b,x)/a
    return 1-bt*betacf(b,a,1-x)/b
ex=[];hh=[];oo=[]
for (x,a,b),ans,o in zip(D["questions"],D["answers"],D["ours"]):
    if a<20 or ans.get("type")!="number": continue
    e=f(ans["bits"]); 
    if e==0 or e==1: continue
    tr=mp.betainc(a,b,0,x,regularized=True); u=math.ulp(e)
    ex.append(float((mp.mpf(e)-tr)/u)); hh.append(float((mp.mpf(hnr(x,a,b))-tr)/u)); oo.append(float((mp.mpf(f(o))-tr)/u))
ex=np.array(ex); hh=np.array(hh); oo=np.array(oo)
print("n",len(ex),"excel err sd %.1f  ours sd %.1f  hnr sd %.1f"%(ex.std(),oo.std(),hh.std()))
print("corr(excel,ours) %.3f corr(excel,hnr) %.3f"%(np.corrcoef(ex,oo)[0,1],np.corrcoef(ex,hh)[0,1]))
m=np.abs(ex)<1e5
print("kept",m.sum(),"excel sd %.1f ours sd %.1f hnr sd %.1f"%(ex[m].std(),oo[m].std(),hh[m].std()))
print("corr(excel,ours) %.3f corr(excel,hnr) %.3f"%(np.corrcoef(ex[m],oo[m])[0,1],np.corrcoef(ex[m],hh[m])[0,1]))
print("median |err| excel %.1f ours %.1f hnr %.1f"%(np.median(np.abs(ex[m])),np.median(np.abs(oo[m])),np.median(np.abs(hh[m]))))
print("excel-ours median |d| %.1f"%np.median(np.abs(ex[m]-oo[m])))
