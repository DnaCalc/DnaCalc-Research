import json, struct, sys, collections
sys.path.insert(0,"/home/ubuntu/codex-runs/r42/pylib")
import mpmath as mp
mp.mp.dps=60
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/BETA.DIST.a1.json"))
f=lambda h: struct.unpack("<d",struct.pack("<Q",int(h,16)))[0]
Q=D["questions"]; E=[f(r["bits"]) for r in D["answers"]]
def rn(v,p):
    if v==0: return mp.mpf(0)
    m,e=mp.frexp(v); s=mp.mpf(2)**p
    return mp.ldexp(mp.nint(m*s),e-p)
c=collections.Counter()
for (x,a,b),e in zip(Q,E):
    x=mp.mpf(x); b=mp.mpf(b)
    Lx=mp.log1p(-x)
    L53=rn(Lx,53); L64=rn(Lx,64)
    cands={"V1":rn(-mp.expm1(rn(b*L53,53)),53),"V2":rn(-mp.expm1(rn(b*L53,64)),53),"V3":rn(-mp.expm1(rn(b*L64,64)),53),
           "V4":rn(-mp.expm1(b*Lx),53),"V5":rn(rn(-mp.expm1(rn(b*L64,64)),64),53),"V6":rn(rn(-mp.expm1(rn(b*L53,64)),64),53)}
    for k,v in cands.items(): c[k]+=float(v)==e
print(len(Q),dict(c))
