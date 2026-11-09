import numpy as np, math
from qr import m87,a87,s87,d87,sq,fmean,EPS
def ls2(cols,y,const,opt):
    n=len(y)
    means=[fmean(c) for c in cols]
    mats=[]
    if const and not opt.get("oneslast"): mats.append(("1",[1.0]*n))
    for i,c in enumerate(cols): mats.append((i,[v-means[i] for v in c] if const else list(c)))
    if const and opt.get("oneslast"):
        mats.append(("1",[1.0]*n)); mats[0],mats[-1]=mats[-1],mats[0]
    scale=[sq(c,0) for _,c in mats]
    j=0; kept=[]
    metric=opt["metric"]; tie=opt["tie"]; later=(tie=="last") or (tie=="const" and const) or (tie=="noconst" and not const)
    while j<len(mats):
        if j>=n: break
        if mats[j][0]!="1":
            def rem(k):
                if metric=="fresh": return sq(mats[k][1],j)
                r=scale[k]
                for e in mats[k][1][:j]: r=s87(r,m87(e,e))
                if metric=="runclip": r=max(r,0.0)
                if metric=="abs": r=abs(r)
                return r
            cands=list(range(j,len(mats)))
            if opt["zero"]=="exclude":
                nz=[k for k in cands if scale[k]>0]
                if nz: cands=nz
            best=cands[0]; br=rem(best)
            for k in cands[1:]:
                r=rem(k)
                if r>br or (later and r==br): br=r; best=k
            if opt["move"]=="swap":
                mats[j],mats[best]=mats[best],mats[j]; scale[j],scale[best]=scale[best],scale[j]
            else:
                m=mats.pop(best); s=scale.pop(best); mats.insert(j,m); scale.insert(j,s)
        col=mats[j][1]; x0=col[j]; rest=sq(col,j+1); squares=sq(col,j)
        before=scale[j]; r=before
        for e in col[:j]: r=r-e*e
        coll=before>0 and r<before*EPS
        if before>0 and mats[j][0]!="1" and "ratios" in opt: opt["ratios"].append((mats[j][0],r/before,j))
        if before>0 and mats[j][0]!="1" and (squares==0 or coll):
            if opt["onstop"]=="stop": break
            del mats[j]; del scale[j]; continue
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
    return [k for k in kept if k!="1"]
