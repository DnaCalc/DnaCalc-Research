import json, struct, sys
sys.path.insert(0,"/home/koderbot/projects/DnaCalc/DnaCalc-Research/racer")
from fractions import Fraction as Q
from intervals import E64, Format
from decimal import Decimal, getcontext
getcontext().prec=60
f=lambda h: struct.unpack(">d",bytes.fromhex(h[2:]))[0]
rows=[]
seen=set()
for fn in ["TRUNC","ROUNDDOWN","ROUNDUP"]:
  for l in open(f"tmp/rd-{fn}-.trace"):
    r=json.loads(l)
    if not r["expected"].startswith("0x"): continue
    x=f(r["args"][0]); d=f(r["args"][1]) if len(r["args"])>1 else 0.0
    if x==0 or x!=x: continue
    e=Decimal(abs(x)).adjusted(); shift=14-e
    if abs(shift)<=27: continue
    if int(d) < shift: continue     # keeps all 15 digits
    ev=abs(f(r["expected"]))
    if ev==0 or ev!=ev or ev==float("inf"): continue
    R=int((Decimal(ev).scaleb(shift)).to_integral_value(rounding="ROUND_HALF_EVEN"))
    key=(x,)
    if key in seen: continue
    seen.add(key)
    rows.append((abs(x),shift,R))
print("rows",len(rows))
def chunks(n,bits=4):
    out=[];sh=0;m=(1<<bits)-1
    while n:
        p=(n&m)<<sh
        if p: out.append(p)
        n>>=bits; sh+=bits
    return out
def rhu(v):  # round half up to integer
    w=v.numerator//v.denominator
    return w+1 if v-w>=Q(1,2) else w
def model(name, x, shift):
    X=Q(x)
    if name=="exact": return rhu(X*Q(10)**shift)
    if name=="old":
        c=E64.rn(Q(10)**abs(shift)); v=E64.rn(X*c) if shift>0 else E64.rn(X/c); return rhu(v)
    if name.startswith("S"):
        prec=int(name[1:3]) if len(name)>1 else 64
        F=Format(prec,-16382,16383)
        v=X
        for p in chunks(abs(shift)):
            c=E64.rn(Q(10)**p if shift>0 else Q(1,10**p))
            v=F.rn(v*c)
        return rhu(E64.rn(v))
for name in sys.argv[1:]:
    ok=sum(1 for x,s,R in rows if model(name,x,s)==R)
    print(name, ok, "of", len(rows))
if "--miss" in sys.argv:
    from collections import Counter
    C=Counter()
    for x,s,R in rows:
        r=model("S64",x,s)
        if r!=R:
            X=Q(x); ex=X*Q(10)**s; w=ex.numerator//ex.denominator
            ulp=Q(2)**(w.bit_length()-64)
            o=model("old",x,s)
            C[(R-r, round(float((ex-w-Q(1,2))/ulp),1) if abs(float((ex-w-Q(1,2))/ulp))<5 else "far", o==R, s>0)]+=1
    for k,v in sorted(C.items(), key=lambda t:-t[1])[:30]: print(k,v)
