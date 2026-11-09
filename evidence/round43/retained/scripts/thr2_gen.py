import json, random, math, struct
from fractions import Fraction as Fr
random.seed(4302)
def h(x): return "0x%016x"%struct.unpack("<Q",struct.pack("<d",float(x)))[0]
def ratio(x,z,const):
    n=len(x); X=[Fr(v) for v in x]; Z=[Fr(v) for v in z]
    if const:
        mx=sum(X)/n; mz=sum(Z)/n; X=[v-mx for v in X]; Z=[v-mz for v in Z]
    sx=sum(v*v for v in X); sz=sum(v*v for v in Z); xz=sum(a*b for a,b in zip(X,Z))
    return (sx-xz*xz/sz)/sx
n=6
x=[3.0,-7.0,12.0,5.0,-11.0,2.0]; w0=[9.0,4.0,-6.0,13.0,1.0,-8.0]; y=[4.0,-2.0,7.0,1.0,-9.0,3.0]
W=[];k=0
lays=[("xz",None),("xz0",None),("xz00",None),("0xz",None)]+[("xzw",s) for s in [-20,-10,-4,0,4,10,20]]+[("xzww",0)]
for const in [True,False]:
  for lay,s in lays:
    for ti in range(49):
        target=-44+ti*0.125
        d0=2.0**-10; z0=[2*v for v in x]; z0[-1]+=d0; c=float(ratio(x,z0,const))/d0**2
        d=math.sqrt(2.0**target/c); e=math.floor(math.log2(d))-9; d=round(d/2**e)*2**e
        z=[2*v for v in x]; z[-1]+=d
        r=float(ratio(x,z,const))
        w=[v*2.0**s for v in w0] if s is not None else None
        w2=[float(v) for v in [5,-3,8,-12,6,2]]
        cols={"x":x,"z":z,"0":[0.0]*n,"w":w}
        X=[]
        for i in range(n):
            row=[]; wi=0
            for ch in lay:
                if ch=="w":
                    row.append((w if wi==0 else w2)[i]); wi+=1
                else: row.append(cols[ch][i])
            X.append(row)
        k+=1
        W.append(dict(id="U%d:4:2"%k,tag="%s-%s-%s"%("c" if const else "nc",lay,s),log2ratio=math.log2(r),
           args=[[[h(v)] for v in y],[[h(v) for v in row] for row in X],const,True],result_index=[4,2]))
json.dump(dict(function="LINEST",witnesses=W),open("/home/ubuntu/codex-runs/r43/q/LINEST.thr2-questions.json","w"))
print(len(W))
