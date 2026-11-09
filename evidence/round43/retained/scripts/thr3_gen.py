import json, random, math, struct
from fractions import Fraction as Fr
random.seed(4303)
def h(x): return "0x%016x"%struct.unpack("<Q",struct.pack("<d",float(x)))[0]
def ratio(x,z,const):
    n=len(x); X=[Fr(v) for v in x]; Z=[Fr(v) for v in z]
    if const:
        mx=sum(X)/n; mz=sum(Z)/n; X=[v-mx for v in X]; Z=[v-mz for v in Z]
    sx=sum(v*v for v in X); sz=sum(v*v for v in Z); xz=sum(a*b for a,b in zip(X,Z))
    return (sx-xz*xz/sz)/sx
W=[];k=0
for n in [4,6,9,16,32]:
  x=[float(random.randint(-30,30)) for _ in range(n)]; w=[float(random.randint(-30,30)) for _ in range(n)]; y=[float(random.randint(-20,20)) for _ in range(n)]
  for const in [True,False]:
    for lay in ["xz","xz0","xzt600","xzt400"]:
      for ti in range(81):
        target=-43+ti/32
        d0=2.0**-10; z0=[2*v for v in x]; z0[-1]+=d0; c=float(ratio(x,z0,const))/d0**2
        d=math.sqrt(2.0**target/c); e=math.floor(math.log2(d))-12; d=round(d/2**e)*2**e
        z=[2*v for v in x]; z[-1]+=d
        r=float(ratio(x,z,const))
        extra=[] if lay=="xz" else ([0.0]*n if lay=="xz0" else [v*2.0**-int(lay[3:]) for v in w])
        cols=[x,z]+([extra] if extra else [])
        X=[[c[i] for c in cols] for i in range(n)]
        k+=1
        W.append(dict(id="V%d:4:2"%k,tag="n%d-%s-%s"%(n,"c" if const else "nc",lay),log2ratio=math.log2(r),
           args=[[[h(v)] for v in y],[[h(v) for v in row] for row in X],const,True],result_index=[4,2]))
json.dump(dict(function="LINEST",witnesses=W),open("/home/ubuntu/codex-runs/r43/q/LINEST.thr3-questions.json","w"))
print(len(W))
