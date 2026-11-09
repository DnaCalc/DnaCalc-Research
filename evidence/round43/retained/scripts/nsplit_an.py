import sys, json, struct, numpy as np
sys.path.insert(0, "/home/ubuntu/codex-runs/r42/pylib")
import mpmath as mp
mp.mp.dps=40
def f(b): return struct.unpack("<d",struct.pack("<Q",b))[0]
D=json.load(open("/home/ubuntu/codex-runs/r43/fresh/NORM.S.INV.split.replies.json"))
for k,v in D.items():
    s0=int(v["start"],16); e=[]
    for j,row in enumerate(v["rows"]):
        p=f(s0+j); x=f(int(row["bits"],16)); t=2*(p-0.5)
        X=mp.sqrt(2)*mp.erfinv(mp.mpf(t))
        e.append(float((mp.mpf(x)-X)/mp.mpf(2)**-52))
    e=np.abs(np.array(e))
    # sliding fraction of |e|>5
    big=(e>5).astype(int)
    cs=[big[i:i+64].mean() for i in range(0,2048,64)]
    print(k, " ".join("%.2f"%c for c in cs))
    # find last index with big in first part, first ... print indices of big errors near middle
    idx=np.nonzero(big)[0]
    print("   big idx range", idx.min() if len(idx) else None, idx.max() if len(idx) else None, "count",len(idx))
