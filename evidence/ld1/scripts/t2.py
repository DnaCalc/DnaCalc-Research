from ask import *
from fractions import Fraction as Q
import json
r,n,pv=0.05/12,360,200000.0
pm = f(ask("PMT",[[r,n,pv,0,0]])["rows"][0]["bits"])
ks=list(range(1,n+1))
cp=[f(x["bits"]) for x in ask("CUMPRINC",[[r,n,pv,k,k,0] for k in ks])["rows"]]
pp=[f(x["bits"]) for x in ask("PPMT",[[r,k,n,pv,0,0] for k in ks])["rows"]]
fvs=[f(x["bits"]) for x in ask("FV",[[r,k,pm,pv,0] for k in range(0,n+1)])["rows"]]
json.dump(dict(pm=pm,cp=cp,pp=pp,fv=fvs),open("t2.json","w"))
R=Q(r);P=Q(pv)
v=1/(1+R)
def ulp(x):
    import math
    return math.ulp(x)
# exact principal k: -pv*r*v^(n-k+1)/(1-v^n)
ex=[ -P*R*v**(n-k+1)/(1-v**n) for k in ks]
same=sum(a==b for a,b in zip(cp,pp))
fvd=sum(cp[k-1]==fvs[k-1]-fvs[k] for k in ks)
print("cp==pp",same,"cp==fv diff",fvd, "of",n)
for k in [1,2,3,10,50,100,180,300,359,360]:
    e=ex[k-1]
    print(k, "cp err ulp %.2f"%(float(Q(cp[k-1])-e)/ulp(cp[k-1])), "pp err %.2f"%(float(Q(pp[k-1])-e)/ulp(pp[k-1])))
import statistics
ce=[abs(float(Q(cp[k-1])-ex[k-1])/ulp(cp[k-1])) for k in ks]
pe=[abs(float(Q(pp[k-1])-ex[k-1])/ulp(pp[k-1])) for k in ks]
print("max cp err",max(ce),"mean",statistics.mean(ce)," pp max",max(pe),statistics.mean(pe))
