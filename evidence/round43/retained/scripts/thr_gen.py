import json, random, math, struct
from fractions import Fraction as Fr
random.seed(4301)
def h(x): return "0x%016x"%struct.unpack("<Q",struct.pack("<d",float(x)))[0]
def ratio(x,z,const):
    n=len(x); X=[Fr(v) for v in x]; Z=[Fr(v) for v in z]
    if const:
        mx=sum(X)/n; mz=sum(Z)/n; X=[v-mx for v in X]; Z=[v-mz for v in Z]
    sx=sum(v*v for v in X); sz=sum(v*v for v in Z); xz=sum(a*b for a,b in zip(X,Z))
    return (sx-xz*xz/sz)/sx
W=[]; k=0
for n in [4,5,6,8]:
  for const in [True,False]:
    for lay in ["xz","zx","xz0","0xz","x0z","xzw","wxz","xz00"]:
      p=len(lay)+int(const)
      if p>n: continue
      for ti in range(17):
        target=-47+0.5*ti
        while True:
            x=[random.randint(-30,30) for _ in range(n)]
            w=[random.randint(-30,30) for _ in range(n)]
            if len(set(x))>2 and len(set(w))>2: break
        # c: ratio at delta=1 (small-delta regime ~ c*delta^2): use delta=2^-10 to measure
        d0=2.0**-10; z0=[2*v for v in x]; z0[-1]+=d0
        r0=ratio(x,z0,const)
        if r0==0: continue
        c=float(r0)/d0**2
        d=math.sqrt(2.0**target/c)
        e=math.floor(math.log2(d))-7; d=round(d/2**e)*2**e
        z=[2.0*v for v in x]; z[-1]=z[-1]+d
        assert z[-1]-2*x[-1]==d
        r=float(ratio(x,z,const))
        cols={"x":[float(v) for v in x],"z":z,"0":[0.0]*n,"w":[float(v) for v in w]}
        X=[[cols[c][i] for c in lay] for i in range(n)]
        y=[float(random.randint(-20,20)) for _ in range(n)]
        k+=1
        W.append(dict(id="T%d:4:2"%k,tag="n%d-%s-%s"%(n,"c" if const else "nc",lay),log2ratio=math.log2(r),
             args=[[[h(v)] for v in y],[[h(v) for v in row] for row in X],const,True],result_index=[4,2]))
json.dump(dict(function="LINEST",witnesses=W),open("/home/ubuntu/codex-runs/r43/q/LINEST.thr-questions.json","w"))
print(len(W))
