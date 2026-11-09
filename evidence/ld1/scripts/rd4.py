import json, struct, sys, itertools
sys.path.insert(0,"/home/koderbot/projects/DnaCalc/DnaCalc-Research/racer")
from fractions import Fraction as Q
from intervals import E64
from decimal import Decimal, getcontext
getcontext().prec=60
f=lambda h: struct.unpack(">d",bytes.fromhex(h[2:]))[0]
rows=[];seen=set()
for fn in ["TRUNC","ROUNDDOWN","ROUNDUP"]:
  for l in open(f"tmp/rd-{fn}-.trace"):
    r=json.loads(l)
    if not r["expected"].startswith("0x"): continue
    x=f(r["args"][0]); d=f(r["args"][1]) if len(r["args"])>1 else 0.0
    if x==0 or x!=x: continue
    e=Decimal(abs(x)).adjusted(); shift=14-e
    if abs(shift)<=27 or int(d)<shift: continue
    ev=abs(f(r["expected"]))
    if not (2.3e-308<ev<1e308): continue
    R=int((Decimal(ev).scaleb(shift)).to_integral_value(rounding="ROUND_HALF_EVEN"))
    if x in seen: continue
    seen.add(x); rows.append((abs(x),shift,R))
def chunks(n):
    out=[];sh=0
    while n:
        p=(n&15)<<sh
        if p: out.append(p)
        n>>=4; sh+=4
    return out
base={}
def const(p,neg):
    k=(p,neg)
    if k not in base: base[k]=E64.rn(Q(1,10**p) if neg else Q(10)**p)
    return base[k]
def ulp64(c):
    return Q(2)**(c.numerator.bit_length()-c.denominator.bit_length()-64+ (1 if c>=Q(2)**(c.numerator.bit_length()-c.denominator.bit_length()) else 0))
def nextc(c,j):
    # move c by j units in 64-bit significand
    for _ in range(abs(j)):
        c=E64.next(c) if j>0 else E64.prev(c)
    return c
off={}
def model(x,shift):
    v=Q(x); neg=shift<0
    for p in chunks(abs(shift)):
        c=const(p,neg)
        j=off.get((p,neg),0)
        if j: c=nextc(c,j)
        v=E64.rn(v*c)
    w=v.numerator//v.denominator
    return w+1 if v-w>=Q(1,2) else w
def score(sub=None):
    rs=rows if sub is None else sub
    return sum(1 for x,s,R in rs if model(x,s)==R)
# restrict to near-tie rows for speed
near=[]
for x,s,R in rows:
    ex=Q(x)*Q(10)**s; w=ex.numerator//ex.denominator
    if abs(ex-w-Q(1,2))< Q(1,10**6)*w/10**8*0+ Q(2)**(w.bit_length()-60): near.append((x,s,R))
print("rows",len(rows),"near",len(near),"base score",score(near), score())
keys=sorted(set((p,s<0) for x,s,R in near for p in chunks(abs(s)) if p>=16))
improved=True
while improved:
    improved=False
    for k in keys:
        best=(score(near),off.get(k,0))
        cur=off.get(k,0)
        for j in (-2,-1,1,2):
            off[k]=j; sc=score(near)
            if sc>best[0]: best=(sc,j)
        off[k]=best[1]
        if best[1]!=cur: improved=True; print("set",k,best, flush=True)
print("final near",score(near),"all",score(),"of",len(rows))
print({k:v for k,v in off.items() if v})
