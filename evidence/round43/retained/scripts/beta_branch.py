import json, collections, struct, math
eps=max(2.220446049250313e-16,1e-15)
def branch(x,a,b):
    y=1-x
    if b==1.0 and x<=0.9: return "unit-b"
    if a==1.0: return "unit-a"
    a0,b0,x0=a,b,x
    if min(a0,b0)<=1.0:
        sw=x>0.5
        if sw: a0,b0,x0=b,a,y
        if b0<min(eps,eps*a0): br="fpser"
        elif a0<min(eps,eps*b0) and b0*x0<=1: br="apser"
        elif max(a0,b0)>1:
            if b0<=1: br="bpser"
            elif x0>=0.3: br="bpser_sym"
            elif x0<0.1 and (x0*b0)**a0<=0.7: br="bpser"
            elif b0>15: br="bgrat_150"
            else: br="bup_140"
        elif a0>=min(0.2,b0): br="bpser"
        elif x0**a0<=0.9: br="bpser"
        elif x0>=0.3: br="bpser_sym"
        else: br="bup_140"
        return "le1:"+br
    lam=(a+b)*y-b if a>b else a-(a+b)*x
    if lam<0: a0,b0,x0=b,a,y; lam=abs(lam)
    if b0<40 and b0*x0<=0.7: br="bpser"
    elif b0<40: br="bup_160"
    elif a0>b0: br="bfrac" if (b0<=100 or lam>0.03*b0) else "basym"
    elif a0<=100: br="bfrac"
    elif lam>0.03*a0: br="bfrac"
    else: br="basym"
    return "gt1:"+br
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.map.json"))
c=collections.Counter()
for (x,a,b),ans,o in zip(D["questions"],D["answers"],D["ours"]):
    if ans.get("type")!="number": continue
    br=branch(x,a,b); c[(br,"n")]+=1; c[(br,"ok")]+=ans["bits"]==o
for br in sorted(set(k[0] for k in c)): print("%-16s %5d/%5d"%(br,c[(br,"ok")],c[(br,"n")]))
