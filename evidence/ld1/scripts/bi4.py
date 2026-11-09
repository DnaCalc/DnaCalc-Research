from ask import *
import json, struct, math
from fractions import Fraction as Q
pat=json.load(open("bi3.json"))
keys=[eval(k) for k in pat]
cum=ask("BINOM.DIST",[[k,n,p,boolean(True)] for (p,n,k) in keys])["rows"]
mass=ask("BINOM.DIST",[[k,n,p,boolean(False)] for (p,n,k) in keys])["rows"]
def nb(x,i): return struct.unpack(">d",struct.pack(">Q",struct.unpack(">Q",struct.pack(">d",x))[0]+i))[0]
M={}
for key,c,m in zip(keys,cum,mass): M[key]=(f(c["bits"]),f(m["bits"]))
res={"pubcdf>":0,"pubcdf>=":0,"summass>":0,"summass>=":0}
bad=[]
for key in keys:
    p,n,k=key
    P=Q(p); c=float(sum(math.comb(n,j)*P**j*(1-P)**(n-j) for j in range(k+1)))
    s=0.0
    for j in range(k+1): s+=M[(p,n,j)][1] if (p,n,j) in M else float(math.comb(n,j)*P**j*(1-P)**(n-j))
    for name,val in [("pubcdf",M[key][0]),("summass",s)]:
        for cmp in (">",">="):
            pp="".join("1" if (val>nb(c,d) if cmp==">" else val>=nb(c,d)) else "0" for d in range(-3,4))
            ok= pp==pat[str(key)]
            res[name+cmp]+=ok
            if name=="pubcdf" and cmp==">=" and not ok: bad.append((key,pat[str(key)],pp))
print(len(keys),res); [print(b) for b in bad]
