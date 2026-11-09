import numpy as np, math
LD=np.longdouble
def m87(a,b): return float(LD(a)*LD(b))
def a87(a,b): return float(LD(a)+LD(b))
def s87(a,b): return float(LD(a)-LD(b))
def d87(a,b): return float(LD(a)/LD(b))
EPS=2.0**-41
def fmean(v):
    s=0.0
    for t in v: s=s+t
    return s/len(v)
def sq(v,a):
    s=0.0
    for t in v[a:]: s=s+m87(t,t)
    return s
def ls(cols,y,const,opt,trace=False):
    """returns (kept predictor ids in order, info). opt keys: move swap|rot, tie: first|last|const (later w/ const), zeroskip, remclip"""
    n=len(y)
    means=[fmean(c) for c in cols]
    mats=[]
    if const: mats.append(("1",[1.0]*n))
    for i,c in enumerate(cols): mats.append((i,[v-means[i] for v in c] if const else list(c)))
    scale=[sq(c,0) for _,c in mats]
    j=0; kept=[]; log=[]
    while j<len(mats):
        if j>=n: break
        if mats[j][0]!="1":
            def rem(k):
                r=scale[k]
                for e in mats[k][1][:j]: r=s87(r,m87(e,e))
                return r
            cands=[k for k in range(j,len(mats))]
            if opt.get("zeroskip"): 
                nz=[k for k in cands if scale[k]>0]
                cands=nz if nz else cands
            best=cands[0]; br=rem(best)
            if opt.get("remclip"): br=max(br,0.0)
            tie=opt.get("tie","const")
            later=(tie=="last") or (tie=="const" and const)
            for k in cands[1:]:
                r=rem(k)
                if opt.get("remclip"): r=max(r,0.0)
                if r>br or (later and r==br): br=r; best=k
            if trace: log.append(("j",j,[(mats[k][0],rem(k)) for k in range(j,len(mats))],"pick",mats[best][0]))
            if opt.get("move","swap")=="swap":
                mats[j],mats[best]=mats[best],mats[j]; scale[j],scale[best]=scale[best],scale[j]
            else:
                m=mats.pop(best); s=scale.pop(best); mats.insert(j,m); scale.insert(j,s)
        col=mats[j][1]; x0=col[j]; rest=sq(col,j+1); squares=sq(col,j)
        before=scale[j]; r=before
        for e in col[:j]: r=r-e*e
        coll=before>0 and r<before*EPS
        if before>0 and mats[j][0]!="1" and (squares==0 or coll):
            if trace: log.append(("stop at",mats[j][0],r,before))
            break
        if squares==0:
            kept.append(mats[j][0]); j+=1; continue
        norm=math.sqrt(squares); alpha=-norm if x0>0 else norm
        v=[s87(x0,alpha)]+col[j+1:]
        vv=a87(rest,m87(v[0],v[0]))
        for k in range(j+1,len(mats)):
            a=mats[k][1]; t=0.0
            for vi,ai in zip(v,a[j:]): t=a87(t,m87(vi,ai))
            fct=d87(m87(2.0,t),vv)
            for i in range(j,n): a[i]=s87(a[i],m87(fct,v[i-j]))
        kept.append(mats[j][0]); j+=1
    kp=[k for k in kept if k!="1"]
    return kp,log
