from ask import *
import json, struct, math
from fractions import Fraction as Q
def nb(x,i): return struct.unpack(">d",struct.pack(">Q",struct.unpack(">Q",struct.pack(">d",x))[0]+i))[0]
rows=[];meta=[]
for p in [0.125,0.25,0.375,0.5,0.625,0.75,0.875]:
    for n in range(1,9):
        P=Q(p); qq=1-P
        cdf=Q(0)
        for k in range(n):
            cdf+= math.comb(n,k)*P**k*qq**(n-k)
            c=float(cdf)
            assert Q(c)==cdf
            for d in range(-3,4):
                a=nb(c,d)
                if 0<a<1:
                    rows.append([n,p,a]); meta.append((p,n,k,d))
res=ask("BINOM.INV",rows)["rows"]
# for each (p,n,k): find the boundary: internal CDF(k) position
from collections import defaultdict
pos=defaultdict(dict)
for (p,n,k,d),x in zip(meta,res):
    pos[(p,n,k)][d]=f(x["bits"])
out={}
for key,v in sorted(pos.items()):
    p,n,k=key
    # answer <=k means CDF(k) "meets" alpha
    meets={d:(v[d]<=k) for d in v}
    # internal threshold: largest d such that meets
    ds=sorted(meets)
    patt="".join("1" if meets[d] else "0" for d in ds)
    out[str(key)]=patt
    print(key, patt)
json.dump(out,open("bi3.json","w"))
